# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import glob
import os
import subprocess
import sys


from netCDF4 import Dataset

from .tempest_common import (TempestExtremesAbstract, _is_date_after,
                             _is_date_after_or_equal)

from tempest_helper import (
    count_trajectories,
    get_trajectories,
    plot_trajectories_cartopy,
    save_trajectories_netcdf,
    storms_overlap_in_time,
    storms_overlap_in_space,
    write_track_line,
    rewrite_track_file
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
            self.processed_files_slp[regrid_resol] = ""
        for track_type in self.track_types:
            self.cmd_detect_type[track_type] = ""
            self.cmd_stitch_type[track_type] = ""
        self.outdir = self.output_directory + "_" + "native"

        self.logger.debug(
            f"CYLC_TASK_CYCLE_TIME {self.cylc_task_cycle_time}, "
            f"um_runid {self.um_runid}"
        )

        timestamp_day = self.cylc_task_cycle_time[:8]
        timestamp_endday = self.next_cycle[:8]
        timestamp_previous = self.previous_cycle[:8]
        timestamp_tm2 = self.tm2_cycle[:8]
        self.logger.debug(f"timestamp_day at top {timestamp_day}")
        self.logger.debug(f"startdate {self.startdate}")
        self.logger.debug(f"enddate {self.enddate}")
        self.logger.debug(f"lastcycle {self.lastcycle}")
        self.logger.debug(f"previous_cycle {self.previous_cycle}")
        self.logger.debug(f"tm2_cycle {self.tm2_cycle}")
        self.logger.debug(f"next_cycle {self.next_cycle}")
        self.logger.debug(f"is_last_cycle {self.is_last_cycle}")

        # First write a dot_file to document which timestamps are yet to be tracked
        dot_file = "do_tracking"
        candidate_files = []
        self._write_dot_track_file(timestamp_day, timestamp_endday, dot_file=dot_file)

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
                ftimestamp_endday = do_track_file.split('.')[1].split('-')[1]

                # do not want to do calculations on data after or equal to the current
                # cycle date, unless it is also the last
                if _is_date_after_or_equal(ftimestamp_day, timestamp_day) \
                        and not self.is_last_cycle == "true":
                    continue

                # if timestamp_previous is before the start date then no work
                if _is_date_after(self.startdate, timestamp_previous):
                    continue

                # find the relevant input data using the given file pattern
                fname = self._file_pattern_processed(ftimestamp_day+"*", "*", "slp",
                                           frequency=self.data_frequency)
                file_search = os.path.join(self.input_directory, fname)
                self.logger.debug(f"file_search {file_search}")

                for regrid_resol in self.regrid_resolutions:
                    self.outdir = self.output_directory+'_'+regrid_resol
                    fname = self._file_pattern_processed(ftimestamp_day + "*", "*", "slp",
                                                         frequency=self.data_frequency)
                    file_search = os.path.join(self.outdir, fname)
                    if glob.glob(file_search):
                        processed_files, variable_units = \
                        self._identify_processed_files(ftimestamp_day, ftimestamp_endday,
                                            grid_resol=regrid_resol)
                        #self.source_files[ftimestamp_day] = source_files
                        self.processed_files[ftimestamp_day] = processed_files
                        self.processed_files_slp[regrid_resol] = \
                            self.processed_files[ftimestamp_day]["slp"]
                        self.variable_units = variable_units

                        # run TempestExtremes detection
                        candidate_files = self._run_detection(ftimestamp_day)
                        # if this timestep has worked OK, then need to remove the dot_file
                        # (the data is needed later)
                        self._remove_dot_track_file(timestamp_previous, timestamp_day)

                    else:
                        self.logger.debug(f"no files to process for timestamp "
                                          f"{ftimestamp_day}")
        else:
            self.logger.error(f"no dot files to process ")

        # find file for timestep name of T-2 (one before timestamp_previous)
        # if this is available, then:
        #     run stitching for the combined T-2&T-1 timestamps
        #     run editing (profile), using the T-2 and T-1 data files for TC profiles,
        #     IKE etc
        fname_slp = self._file_pattern_processed(timestamp_tm2+"*", timestamp_previous+"*", "slp",
                                       frequency=self.data_frequency)
        file_search = os.path.join(self.outdir, fname_slp)
        files_slp_prev = sorted(glob.glob(file_search))
        self.logger.debug(f"file_search for tm2 {file_search} {files_slp_prev} "
                          f"{len(files_slp_prev)}")

        self.logger.debug(f"files_slp_prev {files_slp_prev}")

        if len(files_slp_prev) == 1:
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
                    files_slp_prev[0],
                    grid_resol=regrid_resol
                )

                # archive the track files
                self._archive_track_data(outdir)

        # at this point, I can delete the input data for the T-2 timestep
        if self.delete_processed:
            self._tidy_data_files(timestamp_tm2, timestamp_tm1)

        # if self.delete_source:
        #    self._tidy_data_files(timestamp_previous, timestamp_day,
        #    f_remove = 'source')

    def _run_editing_matching(
        self,
        outdir,
        timestamp_t,
        timestamp_tm1,
        timestamp_tm2,
        tracked_files,
        processed_files_slp,
        grid_resol="native"
    ):
        """
        Run the Tempest tracking (stitch), node editing, matching and writing.

        :param str outdir: The output directory
        :param str timestamp_day: The current timestep of data/tracking to process
        :param str timestamp_previous: The previous timestep of data/tracking
        :                              to process
        :param str processed_files_slp: The path name to the slp nc file
        :param str grid_resol: Either 'native', or the grid resolution for
        :                      regridded output (e.g. N216)
        """

        if len(tracked_files) > 0:
            edited_files = self._run_node_file_editor(
                            tracked_files,
                            timestamp_t,
                            timestamp_tm1,
                            timestamp_tm2,
                            grid_resol=grid_resol)

            self.logger.debug(f"run matching {timestamp_t} {timestamp_tm1} {timestamp_tm2}"
                              f"slp file {processed_files_slp}")
            self._match_tracks(
                outdir,
                timestamp_tm1,
                timestamp_tm2,
                processed_files_slp
                )

            self._collect_and_write_tracks(
                outdir,
                timestamp_tm1,
                timestamp_tm2,
                processed_files_slp
            )

        self._collect_and_write_tracks(
            outdir,
            timestamp_tm1,
            timestamp_tm2,
            processed_files_slp,
            test_new_year=True,
            do_lastcycle=False
        )

        if self.is_last_cycle == 'true':
            self._collect_and_write_tracks(
                outdir,
                timestamp_tm1,
                timestamp_tm2,
                processed_files_slp,
                do_lastcycle=True
            )

    def _collect_and_write_tracks(
        self,
        outdir,
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
                f"{self.um_runid}_{trackname}_????????_????????_{track_type}" +
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
                        f"{self.um_runid}_{trackname}_fullrun_{date_start}_{date_end}_"
                        f"{track_type}.txt"
                    )
                    tracked_file_last = os.path.join(
                        outdir,
                        f"{self.um_runid}_{trackname}_{timestamp_previous}_"
                        f"{timestamp_current}_{track_type}_adjust_f.txt"
                    )
                    if os.path.exists(tracked_file_last):
                        files.extend([tracked_file_last])
                else:
                    tracked_file_final = os.path.join(
                        outdir,
                        f"{self.um_runid}_{trackname}_{date_start}_{date_end}_"
                        f"{track_type}_ongoing.txt"
                    )

            if test_new_year and self._is_new_year:
                do_year = self._old_year_value
                tracked_file_search_year = os.path.join(
                    outdir,
                    f"{self.um_runid}_{trackname}_{do_year}????_????????_"
                    f"{track_type}_adjust_b.txt"
                )
                files_year = sorted(glob.glob(tracked_file_search_year))
                if len(files_year) > 0:
                    for f in files_year:
                        cmd = "mv " + f + " " + os.path.join(outdir,
                                                             self._archived_files_dir)
                        sts = self._run_cmd(cmd, check=True)
                        self.logger.debug(f"{cmd} {sts.stdout}")

                tracked_file_search_year = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.um_runid}_{trackname}_{do_year}????_????????_{track_type}"
                    f"_adjust_b.txt"
                )
                files = sorted(glob.glob(tracked_file_search_year))
                tracked_file_final = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.um_runid}_{trackname}_year_{do_year}_{track_type}.txt"
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
                    title_prefix=f"{self.um_runid} {self.resolution_code} "
                    f"{date_start}_{date_end}",
                    title_suffix=f"{track_type} tracks",
                    write_to_netcdf=True
                )

    def _archive_to_mass(
        self,
        archive_directory
    ):
        """
        Archive any files with .arch to MASS system
        :param str archive_directory: The directory to look for any .arch files
        """
        #TODO this is specific to the Met Office and should be moved out

        moosedir = "moose:/crum/{}/{}/"
        archive_files = glob.glob(os.path.join(archive_directory, '*.arch'))
        if len(archive_files) > 0:
            for fname_arch in archive_files:
                fname = fname_arch[:-5]
                if fname[-2:] == "nc":
                    mass_stream = "any.nc.file"
                else:
                    mass_stream = "ady.file"
                if os.stat(fname).st_size == 0:
                    self.logger.debug(f"File is zero length, no archive {fname}")
                return
                cmd = "moo put -F "+fname+" "+moosedir.format(self.um_suiteid,
                                                             mass_stream)
                self.logger.debug(f"Archive cmd {cmd}")

                sts = subprocess.run(
                    cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                )
                if sts.stderr:
                    msg = (
                        f"Error found in cmd {cmd} {sts.stderr} \n"
                    )
                    raise RuntimeError(msg)
                else:
                    os.remove(fname_arch)
                self.logger.debug(sts.stdout)

    def _archive_track_data(self, outdir):
        """
        Archive the required track data to MASS
        Want to do year files, last cycle. May want to move old files out of way.
        Do we want a rolling tidyup, or just at end of year?
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

            self.logger.debug(f"Archive {self._is_new_year}")
            # archive the annual files
            if self._is_new_year:
                # first move the files to the archive directory, and then archive
                year = self._old_year_value
                subdir = outdir.split('/')[-1]
                candidate_files_search = os.path.join(
                    outdir,
                    f"{self.um_runid}_candidate_{year}????_{track_type}.txt"
                )
                candidate_file_year_tar = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.um_runid}_candidate_{year}_{track_type}.txt.tar.gz"
                )
                files = sorted(glob.glob(candidate_files_search))
                self.logger.debug(f"Files to archive {candidate_files_search} {files}")
                files_to_tar = []
                if len(files) > 0:
                    for f in files:
                        cmd = "mv " + f + " " + os.path.join(outdir,
                                                             self._archived_files_dir)
                        sts = self._run_cmd(cmd, check=True)
                        self.logger.debug(f"{cmd} {sts.stdout}")
                        files_to_tar.append(os.path.join(subdir,
                                                         self._archived_files_dir,
                                                         os.path.basename(f)))
                    os.chdir(os.path.join(outdir, '../'))
                    if os.path.exists(candidate_file_year_tar):
                        os.remove(candidate_file_year_tar)
                    # tar up the candidate files for this year
                    cmd = "tar -cvzf "+candidate_file_year_tar+" "+" ".\
                        join(files_to_tar)
                    self.logger.debug(f"Tar candidates {cmd}")
                    sts = self._run_cmd(cmd, check=False)
                    with open(candidate_file_year_tar + ".arch", "a"):
                        os.utime(candidate_file_year_tar + ".arch", None)
                    #self._archive_to_mass(candidate_file_year_tar)

                tracked_files_adjust_b = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.um_runid}_{trackname}_{year}????_????????_{track_type}"
                    f"_adjust_b.txt"
                )
                tracked_files_adjust_b_tar = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.um_runid}_{trackname}_{year}_{track_type}"
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
                        cmd = "tar -cvzf "+tracked_files_adjust_b_tar+" "+" ".\
                            join(files_to_tar)
                        self.logger.debug(f"Tar candidates {cmd}")
                        sts = self._run_cmd(cmd, check=False)
                        with open(tracked_files_adjust_b_tar + ".arch", "a"):
                            os.utime(tracked_files_adjust_b_tar + ".arch", None)
                        #self._archive_to_mass(tracked_files_adjust_b_tar)

                tracked_file_final = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.um_runid}_{trackname}_year_{year}_{track_type}"
                )
                for f_ending in [".txt", ".nc"]:
                    if os.path.exists(tracked_file_final+f_ending):
                        with open(tracked_file_final + f_ending + ".arch", "a"):
                            os.utime(tracked_file_final + f_ending + ".arch", None)
                        #self._archive_to_mass(tracked_file_final+f_ending)

                candidate_files_to_tidy = os.path.join(
                    outdir,
                    f"{self.um_runid}_candidate_{year}????_????????_{track_type}.txt"
                )
                files = sorted(glob.glob(candidate_files_to_tidy))
                if len(files) > 0:
                    for f in files:
                        #TODO change mv command to a Python os.rename()
                        cmd = "mv " + f + " " + os.path.join(outdir, "tidy")
                        sts = self._run_cmd(cmd, check=True)
                        self.logger.debug(f"{cmd} {sts.stdout}")

                tracked_files_tidy = os.path.join(
                    outdir,
                    f"{self.um_runid}_{trackname}_{year}????_????????_{track_type}.txt"
                )
                files = sorted(glob.glob(tracked_files_tidy))
                if len(files) > 0:
                    for f in files:
                        cmd = "mv " + f + " "  + os.path.join(outdir, "tidy")
                        sts = self._run_cmd(cmd, check=True)
                        self.logger.debug(f"{cmd} {sts.stdout}")

            # archive all the files at the end of the run
            if self.is_last_cycle == 'true':
                pass

        # run the archiving on any .arch files that exist
        self._archive_to_mass(os.path.join(outdir, self._archived_files_dir))

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
                f"{self.um_runid}_candidate_{timestamp}_{track_type}.txt",
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

    def _run_stitching(self, candidate_files):
        """
        Run the Tempest stitching on all of the candidate files.

        :param list candidate_files: The paths to the str paths of the
            candidates to stitch.
        :returns: The path to the candidate file, the path to the tracked file,
            the path of the sea level pressure input netCDF file and , all as
            strings.
        :rtype: tuple
        """
        tracked_files = []

        for index, track_type in enumerate(self.track_types):
            self.logger.debug(f"Running {track_type} stitching")

            candidatefile = candidate_files[index]
            trackedfile = os.path.join(
                self.outdir, f"{self.um_runid}_track_{self.time_range}"
                             f"_{track_type}.txt"
            )
            self._stitch_file(candidatefile, trackedfile, track_type)
            tracked_files.append(trackedfile)

        return tracked_files

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
                f"{self.um_runid}_candidate_{timestamp}_{track_type}.txt",
            )
            candidate_file_previous = os.path.join(
                outdir,
                f"{self.um_runid}_candidate_{timestamp_previous}_{track_type}.txt",
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
                    f"{self.um_runid}_candidate_{timestamp_previous}_{timestamp}"
                    f"_{track_type}.txt",
                )
                with open(candidate_twotimes, "w") as out_file:
                    for candidate_file in candidate_files:
                        with open(candidate_file) as in_file:
                            for line in in_file:
                                out_file.write(line)

                track_twotimes = candidate_twotimes.replace('candidate', 'track')
                self._stitch_file(candidate_twotimes, track_twotimes, track_type)
                tracked_files[track_type] = track_twotimes

        return tracked_files

    def _run_annual_stitch(self, timestamp):
        """
        Concatenate the candidate files for the previous year together and then
        stitch this file.

        :param str timestamp: The timestep of data/tracking to process
        :returns: The string paths to the annual stitched file in a list.
        :rtype: list
        """

        tracked_files = []
        for track_type in self.track_types:
            previous_year = int(timestamp[:4])
            candidate_pattern = os.path.join(
                self.outdir,
                f"{self.um_runid}_candidate_{previous_year}*_{track_type}.txt",
            )
            candidate_files = sorted(glob.glob(candidate_pattern))
            annual_candidate = os.path.join(
                self.outdir,
                f"{self.um_runid}_candidate_year_{previous_year}_{track_type}.txt",
            )
            with open(annual_candidate, "w") as out_file:
                for candidate_file in candidate_files:
                    with open(candidate_file) as in_file:
                        for line in in_file:
                            out_file.write(line)

            annual_track = os.path.join(
                self.outdir, f"{self.um_runid}_track_year_{previous_year}"
                             f"_{track_type}.txt"
            )
            self._stitch_file(annual_candidate, annual_track, track_type)
            tracked_files.append(annual_track)

        return tracked_files

    def _run_wholerun_stitch(self):
        """
        Concatenate the candidate files for all periods together and then
        stitch this file.

        :returns: The string paths to the annual stitched file in a list.
        :rtype: list
        """

        tracked_files = []
        for track_type in self.track_types:
            start_period = int(self.startdate[:8])
            end_period = int(self.enddate[:8])
            candidate_pattern = os.path.join(
                self.outdir,
                f"{self.um_runid}_candidate_????????_{track_type}.txt",
            )
            candidate_files = sorted(glob.glob(candidate_pattern))
            wholerun_candidate = os.path.join(
                self.outdir,
                f"{self.um_runid}_candidate_fullrun_{start_period}_{end_period}"
                f"_{track_type}.txt",
            )
            with open(wholerun_candidate, "w") as out_file:
                for candidate_file in candidate_files:
                    with open(candidate_file) as in_file:
                        for line in in_file:
                            out_file.write(line)

            wholerun_track = os.path.join(
                self.outdir, f"{self.um_runid}_track_fullrun_{start_period}"
                             f"_{end_period}_{track_type}.txt"
            )
            self._stitch_file(wholerun_candidate, wholerun_track, track_type)
            tracked_files.append(wholerun_track)

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

    def _run_node_file_editor(self, tracked_files, timestamp_t, timestamp_tm1, timestamp_tm2,
                              grid_resol='native'):
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

        print('timespamps in run_node_file_editor ', timestamp_t, timestamp_tm1, timestamp_tm2)
        processed_files_prev, variable_units_prev = \
            self._identify_processed_files(timestamp_tm2, timestamp_tm1, grid_resol=grid_resol)
        processed_files_curr, variable_units_curr = \
            self._identify_processed_files(timestamp_tm1, timestamp_t, grid_resol=grid_resol)
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
            with Dataset(nc_file_in, 'r') as x:
                time = x.variables['time']
                calendar = time.calendar
                calendar_units = self.calendar_units
            nc_file_out = self.um_runid+'_'+tracked_file[:-4]+'.nc'
            nc_file_out = os.path.join(os.path.dirname(tracked_file),
                                       os.path.basename(tracked_file)[:-4]+'.nc')
            self.logger.debug(f"open netcdf file {nc_file_out}")
            save_trajectories_netcdf(
                self.outdir,
                nc_file_out,
                storms,
                calendar,
                calendar_units,
                variable_units,
                self.data_frequency,
                self.um_suiteid,
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
                f"{self.um_runid}_{trackname}_{timestamp_previous}_{timestamp_current}"
                f"_{track_type}.txt"
            )
            tracked_file_current_adjust = tracked_file_current[:-4]+'_adjust_f.txt'

            # if the _adjust_f file is available, use this, else the original track file
            tracked_file_previous_search = os.path.join(
                outdir,
                f"{self.um_runid}_{trackname}_????????_{timestamp_previous}"
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
                        storms_space = storms_overlap_in_space(storm_c, storms_time)
                        # self.logger.debug(f"storms_space {storms_space}")
                        if storms_space is not None:
                            storms_time_space_match.append(storms_space)
                            # self.logger.debug(f"Now need to remove this storm from
                            # storms_previous {storms_space}")

                self.logger.debug(f"storms_time_space_match "
                                  f"{len(storms_time_space_match)}")
                rewrite_track_file(
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
        # track_types is a Python list so eval converts str to list
        self.track_types = eval(
            self.app_config.get_property("common", "track_types"))
        self.variables_input = eval(self.app_config.get_property("common",
                                                                 "tc_variables_input"))
        self.variables_rename = eval(self.app_config.get_property("common",
                                                                 "tc_variables_rename"))
