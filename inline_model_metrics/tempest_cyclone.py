# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import glob
import os
import subprocess
import sys
import shutil

from netCDF4 import Dataset

from .tempest_common import (TempestExtremesAbstract, _is_date_after,
                             _is_date_after_or_equal)

from tempest_helper import (
    count_trajectories,
    get_trajectories,
    plot_trajectories_cartopy,
    save_trajectories_netcdf,
    storms_overlap_in_time,
    storm_overlap_in_space,
    write_track_line,
    remove_duplicates_from_track_files
)


class TempestError(Exception):
    """
    Custom TempestExtremes tracking exception
    """
    pass


class TempestExtremesCyclone(TempestExtremesAbstract):
    """
    Run the TempestExtremes cyclone tracker inline with a climate model.
    """

    def __init__(self, arglist=None, **kwargs):
        super().__init__(version="1.0", **kwargs)
        self._parse_args(arglist, desc="Inline TempestExtremes TC tracking")
        self._parse_app_config()
        self._set_message_level()

    @property
    def cli_spec(self):
        """
        Defines the command-line interface specification for the app.
        """
        return [
            {
                "names": ["-c", "--config-file"],
                "help": "Pathname of app configuration file",
            },
        ]

    def run(self, *args, **kwargs):
        """
        Run the app
        """
        self._get_app_options()
        self._get_environment_variables()

        output_dir_native = self.output_directory+'_native'
        if not os.path.exists(output_dir_native):
            try:
                os.makedirs(output_dir_native)
            except PermissionError:
                msg = f"Unable to create output directory {output_dir_native}"
                self.logger.error(msg)
                sys.exit(1)

        # initialise variables that might not get set if no detection step
        for regrid_resol in self.regrid_resolutions:
            self.processed_files_psl[regrid_resol] = ""
        for track_type in self.track_types:
            self.cmd_detect_type[track_type] = ""
            self.cmd_stitch_type[track_type] = ""
        self.outdir = self.output_directory + "_" + "native"

        self.logger.debug(
            f"CYLC_TASK_CYCLE_TIME {self.cylc_task_cycle_time}, "
            f"runid {self.runid}"
        )

        timestamp_current = self.cylc_task_cycle_time[:8]
        timestamp_next = self.next_cycle[:8]
        timestamp_previous = self.previous_cycle[:8]
        timestamp_tm2 = self.tm2_cycle[:8]
        timestamp_tp2 = self.tp2_cycle[:8]
        self.logger.debug(f"timestamp_day at top {timestamp_current}")
        self.logger.debug(f"startdate {self.startdate}")
        self.logger.debug(f"enddate {self.enddate}")
        self.logger.debug(f"lastcycle {self.lastcycle}")
        self.logger.debug(f"previous_cycle {self.previous_cycle}")
        self.logger.debug(f"tm2_cycle {self.tm2_cycle}")
        self.logger.debug(f"next_cycle {self.next_cycle}")
        self.logger.debug(f"is_last_cycle {self.is_last_cycle}")
        self.logger.debug(f"inline_tracking {self.inline_tracking}")

        # Check whether the cylc date is after the most recent post-processed file
        condition = self._get_tracking_date(timestamp_current)
        if condition == "AlreadyComplete":
            return

        # First write a dot_file to document which timestamps are yet to be tracked
        dot_file = "do_tracking"
        candidate_files = []
        self._write_dot_track_file(timestamp_current, timestamp_next, dot_file=dot_file)

        dot_tracking_files = sorted(glob.glob(os.path.join(self.outdir, dot_file+'*')))
        self.logger.debug(f"dot_tracking_files {dot_tracking_files}")

        # loop through all dot files
        # only remove dot file when the detection has run
        # only try to do detection on data with timestamp before the current time
        # (want to make sure that data has been written on the previous step, hence
        # no conflict with writing if postproc on the next step might be running)
        if dot_tracking_files:
            for do_track_file in dot_tracking_files:
                ftimestamp_day = do_track_file.split('.')[1].split("-")[0]
                ftimestamp_endday = do_track_file.split('.')[1].split("-")[1]
                self.logger.debug(f"running dot file {do_track_file} {ftimestamp_day}")
                if self.inline_tracking == "True":
                    self.logger.debug(f"running inline {self.inline_tracking}")
                    # do not want to do calculations on data after or equal to the current
                    # cycle date, unless it is also the last
                    if _is_date_after_or_equal(ftimestamp_day, timestamp_current) \
                            and not self.is_last_cycle == "true":
                        continue

                    # if timestamp_previous is before the start date then no work
                    if _is_date_after(self.startdate, timestamp_previous):
                        continue

                # find the relevant input data using the given file pattern
                #fname = self._file_pattern_processed(ftimestamp_day+"*", "*", "psl",
                #                           frequency=self.data_frequency)
                #file_search = os.path.join(self.input_directory, fname)
                #self.logger.debug(f"file_search {file_search}")

                for regrid_resol in self.regrid_resolutions:
                    self.outdir = self.output_directory+'_'+regrid_resol
                    fname = self._file_pattern_processed(ftimestamp_day + "*", "*",
                                                "psl", frequency=self.data_frequency)
                    file_search = os.path.join(self.outdir, fname)
                    self.logger.debug(f"file_search {file_search}")
                    if glob.glob(file_search):
                        processed_files, variable_units = \
                        self._identify_processed_files(ftimestamp_day,
                                            ftimestamp_endday, grid_resol=regrid_resol)
                        self.processed_files[ftimestamp_day] = processed_files
                        self.processed_files_psl[regrid_resol] = \
                            self.processed_files[ftimestamp_day]["psl"]
                        self.variable_units = variable_units

                        # run TempestExtremes detection
                        candidate_files = self._run_detection(ftimestamp_day)
                        # if this timestep has worked OK, then need to remove the
                        # dot_file (the data is needed later)
                        self._remove_dot_track_file(ftimestamp_day, ftimestamp_endday)

                    else:
                        self.logger.debug(f"no files to process for timestamp "
                                          f"{ftimestamp_day}" f"{file_search}")
        else:
            self.logger.error(f"no dot files to process ")

        # if we run inline, then we need to wait one extra time level, if not then
        # can use time level one later
        if self.inline_tracking == "True":
            timestep_tm2 = timestamp_tm2
            timestep_tm1 = timestamp_previous
            timestep_t = timestamp_current
            timestep_tp1 = timestamp_next
        else:
            timestep_tm2 = timestamp_previous
            timestep_tm1 = timestamp_current
            timestep_t = timestamp_next
            timestep_tp1 = timestamp_tp2

        self._current_year_value = timestep_t[0:4]
        # Run the tracking and node editing for T-2, T-1 files
        self._run_tracking_and_editing(timestep_tm2, timestep_tm1,
                                       timestep_t, do_lastcycle=False)

        # Run the tracking and node editing for T-1, T files if last timestep of run
        if self.is_last_cycle == "true":
            self._run_tracking_and_editing(timestep_tm1, timestep_t,
                                           timestep_tp1, do_lastcycle=True)

    def _run_tracking_and_editing(
        self,
        timestamp_tm2,
        timestamp_previous,
        timestamp_day,
        do_lastcycle
    ):
        """
        Prepare for tracking by setting up required filenames and calling the
        Tempest tracking (stitch), and then call editing of the files to produce
        e.g. profiles if required, and matching tracks across time periods to remove
        duplicates
        Find file for timestep name of T-2 (one before timestamp_previous)
        If this is available, then:
           run stitching for the combined T-2&T-1 timestamps
           run editing (profile), using the T-2 and T-1 data files for TC profiles,
           IKE etc

        :param str timestamp_day: The current timestep of data/tracking to process
        :param str timestamp_previous: The previous timestep of data/tracking
        :                              to process
        :param str timestamp_tm2: The timestep of data/tracking two cycles before now
        :                              to process
        :param bool do_lastcycle: Is this the last cycle in the model run
        """
        fname_psl = self._file_pattern_processed(timestamp_tm2 + "*",
                                                 timestamp_previous + "*", "psl",
                                                 frequency=self.data_frequency)
        file_search = os.path.join(self.outdir, fname_psl)
        files_psl_prev = sorted(glob.glob(file_search))
        self.logger.debug(f"file_search for tm2 {file_search} {files_psl_prev} "
                          f"{len(files_psl_prev)}")

        if len(files_psl_prev) == 1:
            timestamp_tm1 = timestamp_previous

            # run the stitching and edit/profile on native grid data
            for regrid_resol in self.regrid_resolutions:
                outdir = self.output_directory + "_" + regrid_resol
                tracked_files = self._run_stitch_twotimes(outdir, timestamp_tm1,
                                                          timestamp_tm2)

                # run the node editing, and then edit the track files after matching
                # across time periods
                self._run_editing_matching(
                    outdir,
                    timestamp_day,
                    timestamp_tm1,
                    timestamp_tm2,
                    tracked_files,
                    files_psl_prev[0],
                    grid_resol=regrid_resol,
                    do_lastcycle=do_lastcycle
                )

                # archive the track files
                self._archive_track_data(outdir, timestamp_tm2, timestamp_tm1,
                                         do_lastcycle=do_lastcycle)

        # at this point, I can delete the input data for the T-2 timestep
        #if not do_lastcycle:
        #    if self.delete_processed:
        #        self._tidy_data_files(timestamp_tm2,
        #                          timestamp_previous,
        #                          self.variables_rename)

    def _run_editing_matching(
        self,
        outdir,
        timestamp_t,
        timestamp_tm1,
        timestamp_tm2,
        tracked_files,
        processed_files_psl,
        grid_resol="native",
        do_lastcycle=False
    ):
        """
        Run the Tempest node file editing, matching and writing.

        :param str outdir: The output directory
        :param str timestamp_day: The current timestep of data/tracking to process
        :param str timestamp_previous: The previous timestep of data/tracking
        :                              to process
        :param str processed_files_psl: The path name to the psl nc file
        :param str grid_resol: Either 'native', or the grid resolution for
        :                      regridded output (e.g. N216)
        :param bool do_lastcycle: Is this the last cycle in the model run
        """

        if len(tracked_files) > 0:
            edited_files = self._run_node_file_editor(
                            tracked_files,
                            timestamp_t,
                            timestamp_tm1,
                            timestamp_tm2,
                            grid_resol=grid_resol)

            self._match_tracks(
                outdir,
                timestamp_tm1,
                timestamp_tm2,
                processed_files_psl
                )

            self._collect_and_write_tracks(
                outdir,
                timestamp_t,
                timestamp_tm1,
                timestamp_tm2,
                processed_files_psl
            )

        self._collect_and_write_tracks(
            outdir,
            timestamp_t,
            timestamp_tm1,
            timestamp_tm2,
            processed_files_psl,
            test_new_year=True,
            do_lastcycle=False
        )

        if do_lastcycle:
            self._collect_and_write_tracks(
            outdir,
            timestamp_t,
            timestamp_tm1,
            timestamp_tm2,
            processed_files_psl,
            do_lastcycle=do_lastcycle
        )

    def _collect_and_write_tracks(
        self,
        outdir,
        timestamp_next,
        timestamp_current,
        timestamp_previous,
        nc_file_in,
        test_new_year=False,
        do_lastcycle=False
     ):
        """
        Collect together the track files from previous step into updating file,
          with potential to write and plot
        :param str outdir: output directory
        :param str timestamp_day: The current timestep of data/tracking to process
        :param str timestamp_previous: The previous timestep of data/tracking
        :                              to process
        :param str nc_file_in: path to netcdf file for reference calendar information
        :param bool test_new_year: Test whether this is a new year
        :param bool do_lastcycle: Is this the last cycle in the model run
        """
        do_year = None
        for track_type in self.track_types:
            if self._construct_command(track_type)["profile"] is not None:
                trackname = "trackprofile"
                step = "profile"
            else:
                trackname = "track"
                step = "stitch"

            tracked_file_search = os.path.join(
                outdir,
                f"{self.runid}_{trackname}_????????_????????_{track_type}" +
                f"_adjust_b.txt"
            )

            files = sorted(glob.glob(tracked_file_search))
            if len(files) > 0:
                date_start = os.path.basename(files[0]).split('_')[2]
                date_end = os.path.basename(files[-1]).split('_')[3]
                if do_lastcycle:
                    tracked_file_final = os.path.join(
                        outdir,
                        self._archived_files_dir,
                        f"{self.runid}_{trackname}_"
                        f"{date_start}_{timestamp_next}_"
                        f"{track_type}_endrun.txt"
                    )
                    tracked_file_last = os.path.join(
                        outdir,
                        f"{self.runid}_{trackname}_{timestamp_previous}_"
                        f"{timestamp_current}_{track_type}_adjust_f.txt"
                    )
                    if os.path.exists(tracked_file_last):
                        files.extend([tracked_file_last])
                else:
                    tracked_file_final = os.path.join(
                        outdir,
                        f"{self.runid}_{trackname}_{date_start}_{date_end}_"
                        f"{track_type}_ongoing.txt"
                    )

            if test_new_year and self._is_new_year:
                do_year = self._old_year_value
                # all adjust_b files from a given year together make that year of tracks
                tracked_file_search_year = os.path.join(
                    outdir,
                    f"{self.runid}_{trackname}_{do_year}????_????????_"
                    f"{track_type}_adjust_b.txt"
                )
                # files_year = sorted(glob.glob(tracked_file_search_year))
                # if len(files_year) > 0:
                #     for f in files_year:
                #         os.replace(f, os.path.join(outdir, self._archived_files_dir))

                # tracked_file_search_year = os.path.join(
                #     outdir,
                #     self._archived_files_dir,
                #     f"{self.runid}_{trackname}_{do_year}????_????????_{track_type}"
                #     f"_adjust_b.txt"
                # )
                #files = sorted(glob.glob(tracked_file_search_year))
                files = sorted(glob.glob(tracked_file_search_year))
                tracked_file_final = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.runid}_{trackname}_year_{do_year}_{track_type}.txt"
                )
                self.logger.debug(f"This is a new year {do_year} {files}"
                                  f"{test_new_year}")

            self.logger.debug(f"File search for writing tracks {do_year} {files} "
                              f"{test_new_year} {self._is_new_year}")

            if len(files) > 0:
                cmd = "cat "+" ".join(files)+" > "+tracked_file_final
                sts = subprocess.run(
                    cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    check=True,
                )
                self.logger.debug(sts.stdout)
                if sts.stdout:
                    msg = (
                        f"Error found in cat output\n" f"{sts.stdout}"
                    )
                    raise RuntimeError(msg)

                self._read_write_and_plot_tracks(
                    tracked_file_final,
                    track_type,
                    nc_file_in,
                    date_start,
                    date_end,
                    self.variable_units,
                    title_prefix=f"{self.runid} {self.resolution_code} "
                    f"{date_start}_{date_end}",
                    title_suffix=f"{track_type} tracks",
                    write_to_netcdf=True
                )

    def _archive_track_data(self, outdir, timestamp_previous, timestamp_day,
                            do_lastcycle=False):
        """
        Set up archiving of the required data by producing .arch file beside any
        file needing to be archived, and um_postprocess will do archiving
        Want to do year files, last cycle. May want to move old files out of way.
        :param str outdir: output directory
        """
        if not os.path.exists(os.path.join(outdir, self._archived_files_dir)):
            os.makedirs(os.path.join(outdir, self._archived_files_dir))

        for track_type in self.track_types:
            if self._construct_command(track_type)["profile"] is not None:
                trackname = "trackprofile"
                step = "profile"
            else:
                trackname = "track"
                step = "stitch"

            self.logger.debug(f"Archive, is new year? {self._is_new_year}")
            # archive the annual files or remaining files at end of run
            if self._is_new_year or do_lastcycle:
                # first move (if new year) or copy (if end of run) the files to the
                # archive directory, and then add .arch archive indicator
                if self._is_new_year:
                    year = self._old_year_value
                    file_command = "move"
                else:
                    year = self._current_year_value
                    file_command = "copy"
                subdir = outdir.split('/')[-1]
                candidate_files_search = os.path.join(
                    outdir,
                    f"{self.runid}_candidate_{year}????_{track_type}.txt"
                )
                candidate_file_year_tar = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.runid}_candidate_{year}_{track_type}.txt.tar.gz"
                )

                files = sorted(glob.glob(candidate_files_search))
                self.logger.debug(f"Files to archive {candidate_files_search} {files}")
                files_to_tar = []
                if len(files) > 0:
                    for f in files:
                        if file_command == "move":
                            os.replace(f, os.path.join(outdir,
                                        self._archived_files_dir, os.path.basename(f)))
                        else:
                            shutil.copy(f, os.path.join(outdir,
                                        self._archived_files_dir, os.path.basename(f)))
                        files_to_tar.append(os.path.join(subdir,
                                                         self._archived_files_dir,
                                                         os.path.basename(f)))
                    os.chdir(os.path.join(outdir, '../'))
                    if os.path.exists(candidate_file_year_tar):
                        os.remove(candidate_file_year_tar)
                    # tar up the candidate files for this year
                    cmd = "tar -cvzf " + candidate_file_year_tar + " " + " ".\
                        join(files_to_tar)
                    self.logger.debug(f"Tar candidates {cmd}")
                    sts = self._run_cmd(cmd, check=False)
                    with open(candidate_file_year_tar + ".arch", "a"):
                        os.utime(candidate_file_year_tar + ".arch", None)
                    for f in files_to_tar:
                        os.remove(f)

                tracked_file_adjust_b = os.path.join(
                    outdir,
                    f"{self.runid}_{trackname}_{year}????_????????_"
                    f"{track_type}_adjust_b.txt"
                )
                tracked_file_adjust_f = os.path.join(
                    outdir,
                    f"{self.runid}_{trackname}_{timestamp_previous}_{timestamp_day}_"
                    f"{track_type}_adjust_f.txt"
                )
                files_year = sorted(glob.glob(tracked_file_adjust_b))
                if do_lastcycle:
                    if os.path.exists(tracked_file_adjust_f):
                        files_year.extend([tracked_file_adjust_f])

                if len(files_year) > 0:
                    for f in files_year:
                        if file_command == "move":
                            os.replace(f, os.path.join(outdir, self._archived_files_dir, os.path.basename(f)))
                        else:
                            shutil.copy(f, os.path.join(outdir, self._archived_files_dir))

                tracked_files_adjust_b = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.runid}_{trackname}_{year}????_????????_{track_type}"
                    f"_adjust_?.txt"
                )
                tracked_files_adjust_b_tar = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.runid}_{trackname}_{year}_{track_type}"
                    f"_adjust_b.txt.tar.gz"
                )

                files = sorted(glob.glob(tracked_files_adjust_b))
                if len(files) > 0:
                    files_to_tar = []
                    for f in files:
                        files_to_tar.append(os.path.join(subdir,
                                                         self._archived_files_dir,
                                                         os.path.basename(f)))
                    os.chdir(os.path.join(outdir, '../'))
                    if os.path.exists(tracked_files_adjust_b_tar):
                        os.remove(tracked_files_adjust_b_tar)
                    if len(files_to_tar) > 0:
                        # tar up the candidate files for this year
                        cmd = "tar -cvzf " + tracked_files_adjust_b_tar + " " + " ".\
                            join(files_to_tar)
                        self.logger.debug(f"Tar candidates {cmd}")
                        sts = self._run_cmd(cmd, check=False)
                        with open(tracked_files_adjust_b_tar + ".arch", "a"):
                            os.utime(tracked_files_adjust_b_tar + ".arch", None)
                    for f in files_to_tar:
                        os.remove(f)

                tracked_file_final = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.runid}_{trackname}_year_{year}_{track_type}"
                )
                for f_ending in [".txt", ".nc"]:
                    if os.path.exists(tracked_file_final+f_ending):
                        with open(tracked_file_final + f_ending + ".arch", "a"):
                            os.utime(tracked_file_final + f_ending + ".arch", None)

                if self._is_new_year:
                    candidate_files_to_tidy = os.path.join(
                        outdir,
                        f"{self.runid}_candidate_{year}????_????????_"
                        f"{track_type}.txt"
                    )
                    files = sorted(glob.glob(candidate_files_to_tidy))
                    if len(files) > 0:
                        for f in files:
                            os.replace(f, os.path.join(outdir, "tidy",
                                                       os.path.basename(f)))

                    tracked_files_tidy = os.path.join(
                        outdir,
                        f"{self.runid}_{trackname}_{year}????_????????_"
                        f"{track_type}.txt"
                    )
                    files = sorted(glob.glob(tracked_files_tidy))
                    if len(files) > 0:
                        for f in files:
                            os.replace(f, os.path.join(outdir, "tidy",
                                                       os.path.basename(f)))
            print('why not do this bit ',do_lastcycle)
            # archive track the files at the end of the run
            if do_lastcycle:
                do_year = self._current_year_value
                tracked_file_endrun = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.runid}_{trackname}_{do_year}????_????????_{track_type}"
                    f"_endrun*"
                )
                tracked_files = sorted(glob.glob(tracked_file_endrun))
                print('last cycle, endrun ',tracked_files, tracked_file_endrun)
                if len(tracked_files) > 0:
                    for f in tracked_files:
                        f_ending = f[-3:]
                        print('file ending ',f_ending)
                        if f_ending in ["txt", ".nc"]:
                            with open(f + ".arch", "a"):
                                os.utime(f + ".arch", None)

    def _run_detection(self, timestamp):
        """
        Run the Tempest detection.

        :param str timestamp: The timestep of data/tracking to process
        :returns: The path to the candidate files (as a list ordered by track
            type) and details of the processed input files (as a dict).
        :rtype: tuple
        :
        """
        self.logger.debug(f"cwd {os.getcwd()}")

        candidate_files = []

        for track_type in self.track_types:
            self.logger.debug(f"Running {track_type} detection")
            candidatefile = os.path.join(
                self.outdir,
                f"{self.runid}_candidate_{timestamp}_{track_type}.txt",
            )
            self.logger.debug(f"candidatefile {candidatefile}")

            cmd_io = '{} --out {} '.format(
                    self.tc_detect_script,
                    candidatefile)

            # need the first input file not to be orography, since that file has to
            # have a time coordinate
            fnames = []
            for key in self.processed_files[timestamp]:
                if isinstance(self.processed_files[timestamp][key], list):
                    fnames.extend(self.processed_files[timestamp][key])
                else:
                    fnames.append(self.processed_files[timestamp][key])
            fnames.remove(self.processed_files[timestamp]['orog'])

            in_file_list = os.path.join(self.outdir, 'in_file_list_detect.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(fnames) + ';' + \
                           self.processed_files[timestamp]['orog']
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)
            cmd_io += '--in_data_list '+in_file_list+' '

            tracking_phase_commands = self._construct_command(track_type)
            cmd_detect = cmd_io + tracking_phase_commands["detect"]
            self.cmd_detect_type[track_type] = cmd_detect
            self.logger.info(f"Detect command {cmd_detect}")

            sts = subprocess.run(
                cmd_detect,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
            self.logger.debug(sts.stdout)
            if "EXCEPTION" in sts.stdout:
                msg = (
                    f"EXCEPTION found in TempestExtreme detect output\n"
                    f"{sts.stdout}"
                )
                raise RuntimeError(msg)

            candidate_files.append(candidatefile)

        return candidate_files

    def _run_stitch_twotimes(self, outdir, timestamp, timestamp_previous):
        """
        Concatenate the candidate files for the previous year together and then
        stitch this file.

        :param str outdir: directory for the files
        :param str timestamp: The current timestep of data/tracking to process
        :param str timestamp: The previous timestep of data/tracking to process
        :returns: The string paths to the annual stitched file in a list.
        :rtype: list
        """

        tracked_files = {}
        for track_type in self.track_types:
            previous_year = int(timestamp[:4])
            candidate_file_current = os.path.join(
                outdir,
                f"{self.runid}_candidate_{timestamp}_{track_type}.txt",
            )
            candidate_file_previous = os.path.join(
                outdir,
                f"{self.runid}_candidate_{timestamp_previous}_{track_type}.txt",
            )
            self.logger.debug(f"candidate files {candidate_file_previous},"
                              f"{candidate_file_current}")
            if os.path.exists(candidate_file_current) and \
                    os.path.exists(candidate_file_previous):
                candidate_files = [candidate_file_previous]
                candidate_files.append(candidate_file_current)
                self.logger.debug(f"two candidate files {candidate_files}"
                                  f" for stitching")

                candidate_twotimes = os.path.join(
                    outdir,
                    f"{self.runid}_candidate_{timestamp_previous}_{timestamp}"
                    f"_{track_type}.txt",
                )
                with open(candidate_twotimes, "w") as out_file:
                    for candidate_file in candidate_files:
                        with open(candidate_file) as in_file:
                            for line in in_file:
                                out_file.write(line)

                track_twotimes = candidate_twotimes.replace("candidate", "track")
                self._stitch_file(candidate_twotimes, track_twotimes, track_type)
                tracked_files[track_type] = track_twotimes

        return tracked_files

    def _stitch_file(self, candidate_file, tracked_file, track_type):
        """
        Run the Tempest stitching on the specified file.

        :param str candidate_file: The path of the candidate input file.
        :param str tracked_file: The path of the stitched output file.
        :param str track_type: The name of the type of tracking required.
        """

        tracking_phase_commands = self._construct_command(track_type)

        self.logger.debug(f"tracked_file {tracked_file}")

        # stitch candidates together
        cmd_stitch_io = "{} --in {} --out {} ".format(
            self.tc_stitch_script, candidate_file, tracked_file
        )
        cmd_stitch = cmd_stitch_io + tracking_phase_commands["stitch"]
        self.cmd_stitch_type[track_type] = cmd_stitch
        self.logger.info(f"Stitch command {cmd_stitch}")
        sts = subprocess.run(
            cmd_stitch,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=True,
        )
        self.logger.debug(f"sts err {sts.stdout}")

        if "EXCEPTION" in sts.stdout:
            msg = (
                f"EXCEPTION found in TempestExtreme stitch output\n" f"{sts.stdout}"
            )
            raise RuntimeError(msg)

    def _run_node_file_editor(self, tracked_files, timestamp_t, timestamp_tm1,
                              timestamp_tm2, grid_resol='native'):
        """
        Run the Tempest node file editor on all of the tracked files.

        :param list tracked_files: The paths (as strings) of the tracked files
            produced by the stitching process. The files are in order of
            track type.
        :param str timestamp: The timestep of data/tracking to process
        :returns: The path to the candidate file, the path to the tracked file,
            the path of the sea level pressure input netCDF file and , all as
            strings.
        :rtype: tuple
        """

        edited_files = []

        processed_files_prev, variable_units_prev = \
            self._identify_processed_files(timestamp_tm2, timestamp_tm1,
                                           grid_resol=grid_resol)
        processed_files_curr, variable_units_curr = \
            self._identify_processed_files(timestamp_tm1, timestamp_t,
                                           grid_resol=grid_resol)
        self.logger.info(f"processed_files_prev {processed_files_prev}")
        self.logger.info(f"processed_files_curr {processed_files_curr}")

        processed_files_needed = []
        # put the path names of the files with variables needed into list
        for var in self.nodeedit_vars:
            processed_files_needed.append(processed_files_prev[var])
            processed_files_needed.append(processed_files_curr[var])

        self.logger.info(f"processed_files_needed {processed_files_needed}")

        for index, track_type in enumerate(self.track_types):
            tracking_phase_commands = self._construct_command(track_type)["profile"]
            if tracking_phase_commands is None:
                continue

            tracked_file = tracked_files[track_type]
            edited_file = os.path.join(os.path.dirname(tracked_file),
                                       os.path.basename(tracked_file).
                                       replace('track', 'trackprofile'))

            cmd_io = '{} --in_nodefile {} --out_nodefile {} '.format(
                    self.tc_editor_script,
                    tracked_file,
                    edited_file)
            in_file_list = os.path.join(self.outdir, 'in_file_list_edit.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(processed_files_needed)
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)
            cmd_io += '--in_data_list '+in_file_list+' '

            cmd_edit = cmd_io + tracking_phase_commands
            self.cmd_edit_type[track_type] = cmd_edit
            self.logger.info(f"Editor command {cmd_edit}")

            sts = subprocess.run(
                cmd_edit,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
            self.logger.debug(sts.stdout)
            if "EXCEPTION" in sts.stdout:
                msg = (
                    f"EXCEPTION found in TempestExtreme editor output\n"
                    f"{sts.stdout}"
                )
                raise RuntimeError(msg)
            edited_files.append(edited_file)

        return edited_files


    def _read_write_and_plot_tracks(
        self,
        tracked_file,
        track_type,
        nc_file_in,
        timestart,
        timeend,
        variable_units,
        title_prefix="",
        title_suffix="TempestExtremes TCs",
        include_num_tracks=True,
        write_to_netcdf=False
    ):
        """
        Read the tracks, potentially write them to netcdf, and plot them.
        The title is the specified prefix, the
        number of tracks (if requested) and then the suffix.

        :param str tracked_file: The path to the track file.
        :param str nc_file_in: The path to an input data netCDF file, which is
            used to gather additional information about the dates and calendars
            used in the data.
        :param str timestart: The timestep of the start of the data period to process
        :param str timeend: The timestep of the end of the data period to process
        :param dict variable_units: The units of each output variable, as derived from
        :                           the input data files
        :param str title_prefix: The title for the plot.
        :param str title_suffix: The title for the plot.
        :param bool include_num_tracks: Include the number of tracks within the plot
        :                               title
        :param bool write_to_netcdf: Whether to write tracks to netcdf file or not
        """

        # storms returned as dictionary with keys: length, lon, lat, year, month, day,
        # hour, step
        track_step = os.path.basename(tracked_file).split('_')[1]
        if track_step == 'track':
            step = 'stitch'
        elif track_step == 'trackprofile':
            step = 'profile'
        hr_frequency = int(self.data_frequency.strip('h'))
        storms = get_trajectories(tracked_file, nc_file_in, hr_frequency,
                                  self.column_names[track_type+'_'+step])
        self.logger.debug(f"got storms {len(storms)} from {tracked_file}")
        if len(storms) == 0:
            return
        # self.logger.debug(f"got storms {storms[0]}")

        self.logger.debug(f"netcdf file in {nc_file_in}")
        # write the storms to a netcdf file
        if write_to_netcdf:
            with Dataset(nc_file_in, "r") as x:
                time = x.variables["time"]
                calendar = time.calendar
                calendar_units = self.calendar_units
            nc_file_out = self.runid+'_'+tracked_file[:-4]+".nc"
            nc_file_out = os.path.join(os.path.dirname(tracked_file),
                                       os.path.basename(tracked_file)[:-4]+".nc")
            self.logger.debug(f"open netcdf file {nc_file_out}")
            save_trajectories_netcdf(
                self.outdir,
                nc_file_out,
                storms,
                calendar,
                calendar_units,
                variable_units,
                self.data_frequency,
                self.suiteid,
                self.resolution_code,
                self.cmd_detect_type[track_type],
                self.cmd_stitch_type[track_type],
                self.column_names[track_type+'_'+step],
                startperiod=timestart,
                endperiod=timeend)

        title_components = []
        if title_prefix:
            title_components.append(title_prefix)
        if include_num_tracks:
            title_components.append(str(count_trajectories(storms)))
        if title_suffix:
            title_components.append(title_suffix)

        title_full = " ".join(title_components)

        filename = tracked_file[:-4] + ".png"
        self.logger.debug(f"plot storms {filename}")
        if self.plot_tracks:
            try:
                plot_trajectories_cartopy(storms, filename, title=title_full)
            except:
                self.logger.debug(f"Error from plotting trajectories ")

    def _match_tracks(
        self,
        outdir,
        timestamp_current,
        timestamp_previous,
        nc_file_in
    ):
        """
        Read the tracks, potentially write them to netcdf, and plot them.
        The title is the specified prefix, the
        number of tracks (if requested) and then the suffix.

        :param str outdir: output directory
        :param str nc_file_in: The path to an input data netCDF file, which is
            used to gather additional information about the dates and calendars
            used in the data.
        :param str timestamp_current: The current timestamp of the file to process
        :param str timestamp_previous: The previous timestamp of the file to process
        """

        for track_type in self.track_types:
            if self._construct_command(track_type)["profile"] is not None:
                trackname = 'trackprofile'
                step = 'profile'
            else:
                trackname = 'track'
                step = 'stitch'
            column_names = self.column_names[track_type+'_'+step]

            tracked_file_current = os.path.join(
                outdir,
                f"{self.runid}_{trackname}_{timestamp_previous}_{timestamp_current}"
                f"_{track_type}.txt"
            )
            tracked_file_current_adjust = tracked_file_current[:-4]+'_adjust_f.txt'

            # if the _adjust_f file is available, use this, else original track file
            tracked_file_previous_search = os.path.join(
                outdir,
                f"{self.runid}_{trackname}_????????_{timestamp_previous}"
                f"_{track_type}_adjust_f.txt"
            )
            if len(glob.glob(tracked_file_previous_search)) == 1:
                tracked_file_previous = glob.glob(tracked_file_previous_search)[0]
                tracked_file_previous_adjust = \
                    glob.glob(tracked_file_previous_search)[0].\
                    replace('_adjust_f', '_adjust_b')
            elif len(glob.glob(tracked_file_previous_search.replace('_adjust_f', '')))\
                    == 1:
                tracked_file_previous = \
                    glob.glob(tracked_file_previous_search.replace('_adjust_f', ''))[0]
                tracked_file_previous_adjust = tracked_file_previous[:-4] +\
                    '_adjust_b.txt'
            else:
                self.logger.debug("No previous step track file"
                                  "{tracked_file_previous_search.format('_adjust_b')}"
                                  "{tracked_file_previous_search.format('')}")
                continue

            self.logger.debug(f"tracked_files2 {tracked_file_current},"
                              f" {tracked_file_previous}")

            if os.path.exists(tracked_file_previous) and \
                    os.path.exists(tracked_file_current):
                self.logger.debug(f"_match_tracks nc_file_in {nc_file_in}")
                hr_frequency = int(self.data_frequency.strip('h'))
                storms_previous = get_trajectories(tracked_file_previous,
                                                   nc_file_in,
                                                   hr_frequency,
                                                   column_names)
                storms_current = get_trajectories(tracked_file_current,
                                                  nc_file_in,
                                                  hr_frequency,
                                                  column_names)

                storms_time_space_match = []
                for storm_c in storms_current:
                    storms_time = storms_overlap_in_time(storm_c, storms_previous)
                    if len(storms_time) > 0:
                        self.logger.debug(f"storms overlap in time {len(storms_time)}")
                        # self.logger.debug(f"storms_time {storms_time}")
                        # self.logger.debug(f"storm_c {storm_c}")
                        storms_space = storm_overlap_in_space(storm_c, storms_time)
                        # self.logger.debug(f"storms_space {storms_space}")
                        if storms_space is not None:
                            storms_time_space_match.append(storms_space)
                            # self.logger.debug(f"Now need to remove this storm from
                            # storms_previous {storms_space}")

                self.logger.debug(f"storms_time_space_match "
                                  f"{len(storms_time_space_match)}")
                remove_duplicates_from_track_files(
                    tracked_file_previous,
                    tracked_file_current,
                    tracked_file_previous_adjust,
                    tracked_file_current_adjust,
                    storms_time_space_match,
                    column_names
                )

                # if this time period crossing into a new year (wrt the adjust_b file)
                year1 = os.path.basename(tracked_file_previous_adjust).\
                    split('_')[2][0:4]
                year2 = os.path.basename(tracked_file_previous_adjust).\
                    split('_')[3][0:4]
                self._current_year_value = year2
                if year1 != year2:
                    self._is_new_year = True
                    self._old_year_value = year1
                self._tidy_track_files(
                    outdir,
                    tracked_file_previous
                )

    def _get_app_options(self):
        """Get commonly used configuration items from the config file"""

        super()._get_app_options()

        self.tc_detect_script = self.app_config.get_property(
            "common", "tc_detect_script"
        )
        self.tc_stitch_script = self.app_config.get_property(
            "common", "tc_stitch_script"
        )
        self.tc_editor_script = self.app_config.get_property(
            "common", "tc_editor_script"
        )
        self.nodeedit_vars = eval(self.app_config.get_property("common",
                                                               "nodeedit_vars"))
        # track_types is a Python list so eval converts str to list
        self.track_types = eval(
            self.app_config.get_property("common", "track_types"))
        self.variables_input = eval(self.app_config.get_property("common",
                                                            "tc_variables_input"))
        self.variables_rename = eval(self.app_config.get_property("common",
                                                            "tc_variables_rename"))
