# (C) British Crown Copyright 2020, Met Office.
# Please see LICENSE for license details.
import glob
import os
import re
import shutil
import subprocess
import sys
from shutil import copyfile, move

from netCDF4 import Dataset
import numpy as np
from cftime import utime, datetime
import iris

from afterburner.apps import AbstractApp

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


class TempestTracker(AbstractApp):
    """
    Run the TempestExtreme tracker inline with a climate model
    """

    def __init__(self, arglist=None, **kwargs):
        super().__init__(version="1.0", **kwargs)
        self._parse_args(arglist, desc="Inline TempestExtreme tracking")
        self._parse_app_config()
        self._set_message_level()

        # Instance attributes set later
        self.time_range = None
        self.frequency = None
        self.resolution_code = None
        self.cmd_detect_type = {}
        self.cmd_stitch_type = {}
        self.cmd_edit_type = {}
        self.cmd_detect = None
        self.cmd_stitch = None
        self.cmd_edit = None
        self.source_files = {}
        self.processed_files = {}
        self.processed_files_slp = {}
        self.variable_units = {}
        self.calendar_units = "days since 1950-01-01 00:00:00"
        self.regrid_resolutions = None
        self.outdir = None
        self.column_names = {}
        self._is_new_year = False
        self._old_year_value = None
        self._archived_files_dir = "archived_files"
        self.variables_rename = ["slp", "uas", "vas", "viwve", "viwvn", "rvT42_850", "rh_850", "ts", \
                                 "zg_925", "zg_850", "zg_700", "zg_600", "zg_500", "zg_250", \
                                 "rvT63_850", "rvT63_700", "rvT63_600", "rvT63_500", "rvT63_250", \
                                 "ua_925", "ua_850", "va_925", "va_850"]
        self.outputcmd_detect_default = ""
        self.in_fmt_stitch_default = ""
        self.out_fmt_profile1_default = ""
        self.out_fmt_profile2_default = ""

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

#        # this section of code processes data from the current timestep
#        current_time = timestamp_day
#        # find the relevant input data using the given file pattern
#        fname = self._file_pattern(current_time + "*", "*", "slp",
#                                   um_stream="pt", frequency="*")
#        file_search = os.path.join(self.input_directory, fname)
#        self.logger.debug(f"file_search {file_search}")

#        # pre-process the input files and produce standard processed files for later use
#        if glob.glob(file_search):
#            for regrid_resol in self.regrid_resolutions:
#                self.outdir = self.output_directory + '_' + regrid_resol
#                source_files, processed_files, variable_units = \
#                    self._generate_data_files(timestamp_day, timestamp_endday,
#                                              grid_resol=regrid_resol)

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
                if self._is_date_after_or_equal(ftimestamp_day, timestamp_day) \
                        and not self.is_last_cycle == "true":
                    continue

                # if timestamp_previous is before the start date then no work
                if self._is_date_after(self.startdate, timestamp_previous):
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
                self._archive_track_data(
                    outdir,
                    timestamp_tm1,
                    timestamp_tm2
                )

        # at this point, I can delete the input data for the T-2 timestep
        if self.delete_processed:
            self._tidy_data_files(timestamp_tm2, timestamp_tm1)

        # if self.delete_source:
        #    self._tidy_data_files(timestamp_previous, timestamp_day,
        #    f_remove = 'source')

    def _run_cmd(
        self,
        cmd,
        check=True
    ):
        sts = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=check,
        )
        self.logger.debug(sts.stdout)
        if sts.stderr:
            msg = (
                f"Error found in run_cmd output\n" f"{sts.stderr}"
            )
            raise RuntimeError(msg)
        return sts

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

    def _archive_track_data(
        self,
        outdir,
        timestamp_current,
        timestamp_previous
    ):
        """
        Archive the required track data to MASS
        Want to do year files, last cycle. May want to move old files out of way.
        Do we want a rolling tidyup, or just at end of year?
        :param str outdir: output directory
        :param str timestamp_current: the current timestep of data/tracking to process
        :param str timestamp_previous: the previous timestep of data/tracking
        :                              to process
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
                    with open(os.path.join(candidate_file_year_tar, ".arch"), "a"):
                        os.utime(os.path.join(candidate_file_year_tar, ".arch"), None)
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
                        with open(os.path.join(tracked_files_adjust_b_tar, ".arch"), "a"):
                            os.utime(os.path.join(tracked_files_adjust_b_tar, ".arch"), None)
                        #self._archive_to_mass(tracked_files_adjust_b_tar)

                tracked_file_final = os.path.join(
                    outdir,
                    self._archived_files_dir,
                    f"{self.um_runid}_{trackname}_year_{year}_{track_type}"
                )
                for f_ending in [".txt", ".nc"]:
                    if os.path.exists(tracked_file_final+f_ending):
                        with open(os.path.join(tracked_file_final+f_ending, ".arch"), "a"):
                            os.utime(os.path.join(tracked_file_final+f_ending, ".arch"), None)
                        #self._archive_to_mass(tracked_file_final+f_ending)

                candidate_files_to_tidy = os.path.join(
                    outdir,
                    f"{self.um_runid}_candidate_{year}????_????????_{track_type}.txt"
                )
                files = sorted(glob.glob(candidate_files_to_tidy))
                if len(files) > 0:
                    for f in files:
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

    def _is_new_year_old(self, timestamp, timestamp_end):
        """
        Check if this is the first cycle in a new year.

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :returns: True if this is this period crosses a year boundary.
        :rtype: bool
        """

        if timestamp[:4] != timestamp_end[:4]:
            return True
        else:
            return False

    def _is_date_after(self, timetest, timeref):
        """
        Check if timetest is after timeref.

        :returns: True if timetest is after timeref.
        :rtype: bool
        """

        if int(timetest[:8]) > int(timeref[:8]):
            return True
        else:
            return False

    def _is_date_after_or_equal(self, timetest, timeref):
        """
        Check if timetest is equal or after timestamp.

        :returns: True if timetest is equal or after timeref.
        :rtype: bool
        """

        if int(timetest[:8]) >= int(timeref[:8]):
            return True
        else:
            return False

    def _write_dot_track_file(self, timestamp, timestamp_end,
                              dot_file='do_tracking'):
        """
        Write a file indicating that this timestep needs to be tracked

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param str dot_file: The first part of the string of a filename to indicate
        :                    which time periods still need tracking

        """
        do_tracking_file = os.path.join(self.outdir, dot_file+'.' + timestamp +
                                        '-'+timestamp_end)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if not os.path.exists(do_tracking_file):
            os.system('touch '+do_tracking_file)

    def _remove_dot_track_file(self, timestamp, timestamp_end,
                               dot_file='do_tracking'):
        """
        Remove a file indicating that this timestep needs to be tracked

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param str dot_file: The first part of the string of a filename to indicate
        :                    which time periods still need tracking

        """
        do_tracking_file = os.path.join(self.outdir, dot_file+'.' + timestamp +
                                        '-'+timestamp_end)
        if os.path.exists(do_tracking_file):
            os.system('rm '+do_tracking_file)

    def _tidy_data_files(self, timestamp, timestamp_end, dot_file='do_tracking',
                         f_remove='processed'):
        """
        Remove input files and tracking dot file for this timestamp (tidy up)

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param str dot_file: The first part of the string of a filename to indicate
        :                    which time periods still need tracking
        :param str f_remove: An indicator of which files need to be deleted.
                   processed = the (regridded) files read by the tracking code
                   source    = the files produced by the model
        """
        self.logger.info(f"Tidy up input files")
        files_remove = []
        source_files, processed_files, variable_units = \
                    self._generate_data_files(timestamp)

        if f_remove == 'processed':
            files_remove = processed_files
            files_remove.pop('orog', None)
        elif f_remove == 'source':
            if self.delete_source:
                files_remove = source_files
                files_remove.pop('orog', None)

        for v in files_remove:
            self.logger.info(f"Tidy data remove {files_remove[v]}")
            os.remove(files_remove[v])

        # if f_remove == 'processed':
        #    self._remove_dot_track_file(timestamp, timestamp_end, dot_file = dot_file)
        # self.logger.debug(f"removed dot file {timestamp}")

    def _tidy_track_files(
        self,
        outdir,
        track_file
    ):
        """
        Tidy up requested files - initially move to subdir, eventually delete
        """
        tidy_dir = os.path.join(outdir, 'tidy')
        if not os.path.exists(tidy_dir):
            os.makedirs(tidy_dir)
        self.logger.info(f"Tidy track files mv {track_file} to tidy dir")
        # move(track_file, tidy_dir)
        # cmd = 'mv '+ track_file+ ' '+tidy_dir
        # os.system(cmd)

    def _file_pattern(self, timestart, timeend, varname, um_stream='pt',
                      frequency='6h'):
        """
        Derive the input nc filenames from the file pattern, assuming a
        um model filenaming pattern as here, or could be other patterns
        for other models/platforms (which would need to be added)

        :param str timestart: The timestep of the start of the data period to process
        :param str timeend: The timestep of the end of the data period to process
        :param str um_stream: The name of the um output stream (output file
        :                     identification)
        :param str frequency: The frequency of the input data (in hours, needs to
        :                     include "h"), used to determine file naming
        :returns: a filename given the inputs to the pattern
        :rtype: str
        """
        if self.data_frequency is None:
            file_freq = frequency
        else:
            file_freq = str(self.data_frequency)

        if self.um_file_pattern != '':
            if 'atmos' in self.um_file_pattern:
                # file format from postproc
                fname = self.um_file_pattern.format(
                    runid=self.um_runid,
                    frequency=file_freq,
                    date_start=timestart,
                    date_end=timeend,
                    stream=um_stream,
                    variable=varname
                )
            else:
                # file format from direct STASH to netcdf conversion
                fname = self.um_file_pattern.format(
                    runid=self.um_runid,
                    stream=um_stream,
                    date_start=timestart,
                    variable=varname
                )
        self.logger.info(f"fname from pattern {fname} {um_stream} {timestart} {timeend} {varname}")
        return fname.strip('"')

    def _file_pattern_processed(self, timestart, timeend, varname,
                      frequency='6h'):
        """
        For processed files, we know what the filenames look like, so search specifically

        :param str timestart: The timestep of the start of the data period to process
        :param str timeend: The timestep of the end of the data period to process
        :param str um_stream: The name of the um output stream (output file
        :                     identification)
        :param str frequency: The frequency of the input data (in hours, needs to
        :                     include "h"), used to determine file naming
        :returns: a filename given the inputs to the pattern
        :rtype: str
        """
        if self.data_frequency is None:
            file_freq = frequency
        else:
            file_freq = str(self.data_frequency)

        fname = self.file_pattern_processed.format(
            runid=self.um_runid,
            frequency=file_freq,
            date_start=timestart,
            date_end=timeend,
            variable=varname
        )

        self.logger.info(f"fname from pattern {fname} {timestart} {timeend} {varname}")
        return fname.strip('"')

    def _construct_command(self, track_type):
        """
        Read the TempestExtreme command line parameters from the configuration.
        :param str track_type: The name of the type of tracking to run, possible
        :                      values:
               detect, stitch, profile
        :returns: A dictionary with keys of `detect`, `stitch`, `profile`
            and the values for each of these is a string containing the
            command line parameters for each of these TempestExtreme steps.
            The parameters are sorted into alphabetical order in each line.
        :rtype: dict
        """

        # These are the first and last column hesders in the standard
        # Tempest output files
        column_initial = "grid_x,grid_y,"
        column_final = ",year,month,day,hour"

        commands = {}
        fmt_value = {}
        for step in ["detect", "stitch", "profile", "detectblobs", "nodefilefilter"]:
            try:
                step_config = self.app_config.section_to_dict(f"{track_type}_{step}")
                #step_arguments = [
                #    f"--{parameter} {step_config[parameter]}"
                #    for parameter in sorted(list(step_config.keys()))
                #    if step_config[parameter]
                #]
                step_arguments = []
                for parameter in sorted(list(step_config.keys())):
                    if step_config[parameter]:
                        print('param ',step, parameter, step_config[parameter])
                        if '_default' in step_config[parameter]:
                            if "outputcmd" in step_config[parameter]:
                                param_value = self.outputcmd_detect_default
                            elif "in_fmt_stitch" in step_config[parameter]:
                                param_value = self.in_fmt_stitch_default
                                fmt_value['stitch'] = self.in_fmt_stitch_default
                            elif "out_fmt_profile1" in step_config[parameter]:
                                param_value = self.out_fmt_profile1_default
                                fmt_value['profile'] = self.out_fmt_profile1_default
                            elif "out_fmt_profile2" in step_config[parameter]:
                                param_value = self.out_fmt_profile2_default
                                fmt_value['profile'] = self.out_fmt_profile2_default
                            print('_default in ', step_config[parameter], param_value)
                        else:
                            param_value = step_config[parameter]
                            if "in_fmt" in parameter:
                                fmt_value["stitch"] = step_config[parameter]
                            elif "out_fmt" in parameter:
                                fmt_value["profile"] = step_config[parameter]

                        step_arguments.append(f"--{parameter} {param_value}")


                commands[step] = " ".join(step_arguments)
                print('step, commands ',step, commands[step])
            except:
                commands[step] = None

            # set up the column names of the track output file, to be used for
            # naming the storm keys
            if step == "stitch" and commands[step] is not None:
                #col_names = column_initial + step_config["in_fmt"].strip('\"') +\
                #            column_final
                col_names = column_initial + fmt_value["stitch"].strip('\"') +\
                            column_final
                self.column_names[track_type+"_stitch"] = {}
                names = col_names.split(',')
                for im, name in enumerate(names):
                    self.column_names[track_type+"_"+step][name] = im
            if step == "profile" and commands[step] is not None:
                #col_names = column_initial + step_config["out_fmt"].strip('\"') +\
                #            column_final
                col_names = column_initial + fmt_value["profile"].strip('\"') +\
                            column_final
                self.column_names[track_type+"_profile"] = {}
                names = col_names.split(',')
                for im, name in enumerate(names):
                    self.column_names[track_type+"_"+step][name] = im

        return commands

    def _generate_data_files(self, timestamp, timestamp_end, grid_resol='native'):
        """
        Identify and then fix the grids and var_names in the input files.
        The time_range and frequency attributes are set when this method runs.

        :param str timestamp: The timestep of the start of the data period to process
        :param str grid_resol: Either native, or the resolution string if regridding
        :                      is required
        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """
        timestamp_day = timestamp
        self.logger.debug(f"timestamp_day in generate_data_files {timestamp_day}")
        fname = self._file_pattern(timestamp_day+'*', '*', 'slp', um_stream='pt',
                                   frequency='*')
        file_search = os.path.join(
            self.input_directory, fname
        )
        files = sorted(glob.glob(file_search))
        self.logger.debug(f"file_search {file_search}, files {files}")
        if not files:
            msg = f"No input files found for glob pattern {file_search}"
            self.logger.error(msg)
            raise RuntimeError(msg)

        time_ranges = []
        unique_ranges = []
        frequencies = []
        unique_frequencies = []
        for filename in files:
            file_elements = os.path.basename(filename).split("_")
            if len(file_elements) >= 3:
                time_range = os.path.basename(filename).split("_")[3]
                frequency = os.path.basename(filename).split("_")[2]
            else:
                time_range = timestamp_day+'-'+timestamp_day
                frequency = self.data_frequency
            time_ranges.append(time_range)
            frequencies.append(frequency)
            unique_ranges = list(set(time_ranges))
            unique_frequencies = list(set(frequencies))
        if len(unique_ranges) == 0:
            msg = "No tracked_file periods found"
            self.logger.error(msg)
            raise RuntimeError(msg)
        elif len(unique_ranges) != 1:
            msg = "No tracked_file periods found"
            self.logger.error(msg)
            raise RuntimeError(msg)
        else:
            self.time_range = unique_ranges[0]
        if len(unique_frequencies) == 0:
            msg = "No tracked_file frequencies found"
            self.logger.error(msg)
            raise RuntimeError(msg)
        elif len(unique_frequencies) != 1:
            msg = "No tracked_file frequencies found"
            self.logger.error(msg)
            raise RuntimeError(msg)
        else:
            frequency = unique_frequencies[0]
            components = re.match(r"(\d+)", frequency)
            if not components:
                msg = r"No digit found in frequency {frequency}"
                self.logger.error(msg)
                raise ValueError(msg)
            else:
                self.frequency = int(components[1])

        source_files, processed_files, variable_units = \
            self._process_input_files(timestamp, timestamp_end, grid_resol=grid_resol)

        return source_files, processed_files, variable_units

    def _check_lat_lon_coords(self, fname, varname):
        """
        Check that file has latitude and longitude coordinates called
        latitude, longitude, if not change the names

        :param str fname: filename in which to check coord names
        """
        cube = iris.load_cube(fname)
        lat_name = cube.coord('latitude').var_name
        lon_name = cube.coord('longitude').var_name
        time_name = cube.coord('time').var_name
        if lat_name != 'latitude' or lon_name != 'longitude':
            cmd = os.path.join(self.ncodir, "ncks") \
                  + " -O -6 " + fname + " " + fname + ".nc3"
            self.logger.debug(f"cmd {cmd}")
            subprocess.call(cmd, shell=True)
            for old, new in zip([lat_name, lon_name], ["latitude", "longitude"]):
                cmd = os.path.join(self.ncodir, "ncrename") + \
                      " -d " + old + ',' + new + ' -v ' + old + ',' + new + ' ' + fname+'.nc3'
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)
            cmd = "mv " + fname + ".nc3" + " " + fname
            self.logger.debug(f"cmd {cmd}")
            subprocess.call(cmd, shell=True)

        if cube.var_name != varname:
            cmd = os.path.join(self.ncodir, "ncrename") + \
                  " -v " + cube.var_name + "," + varname + " " + fname
            self.logger.debug(f"cmd {cmd}")
            subprocess.call(cmd, shell=True)

    def _check_time_coord(self, fnames):
        """
        Check that file has latitude and longitude coordinates called
        latitude, longitude, if not change the names

        :param list fnames: filenames in which to check coord names
        """
        print('regrid ',fnames)
        for fname in fnames:
            cube = iris.load_cube(fname)
            time_name = cube.coord('time').var_name
            if time_name != 'time':
                cmd = os.path.join(self.ncodir, "ncks") \
                      + " -6 " + fname + " " + fname + ".nc3"
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)
                cmd = os.path.join(self.ncodir, "ncrename") + \
                        " -d " + time_name + ",time -v " + time_name + ",time " + fname+".nc3"
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)

                cmd = "mv " + fname + ".nc3" + " " + fname
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)

    def _process_input_files(self, time_start, time_end, grid_resol='native'):
        """
        Identify and then fix the grids and var_names in the input files.
        The variable names need to have a new var_name, either because the default
        from the UM is confusing or unknown, and these names are needed for the
        variable name inputs for the TempestExtremes scripts

        :param str grid_resol: The resolution string to be used if regridding
        :                      is required
        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """
        source_filenames = {}
        processed_filenames = {}
        variable_units = {}

        variables_required = {}
        # these variables need to have a new var_name, either because the default
        # from the UM is confusing or unknown, and these names are needed for the
        # variable name inputs for the TempestExtremes scripts
        variables_rename = self.variables_rename
        for var in self.variables_input:
            variables_required[var] = {'fname': var}
            if var in variables_rename:
                variables_required[var].update({'varname_new': var})

        reference_name = self._file_pattern(self.time_range.split('-')[0],
                                            self.time_range.split('-')[1],
                                            variables_required["slp"]["fname"],
                                            um_stream='pt')
        reference_path = os.path.join(self.input_directory, reference_name)
        reference = iris.load_cube(reference_path)
        variable_units['slp'] = reference.units

        # Identify the grid and orography file
        if grid_resol == 'native':
            longitude_size = reference.shape[-1]
            resolution = longitude_size // 2
            self.resolution_code = f"N{resolution}"

            processed_filenames["orog"] = os.path.join(
                self.orography_dir, f"orog_HadGEM3-GC31-{self.resolution_code}e.nc"
            )
            shutil.copyfile(processed_filenames["orog"], os.path.join(self.outdir, 'orography.nc'))
        else:
            resolution_code = f"N{grid_resol}"
            processed_filenames["orog"] = os.path.join(
                self.orography_dir, f"orog_HadGEM3-GC31-{grid_resol}e.nc"
            )
        cube_orog = iris.load_cube(processed_filenames["orog"])
        variable_units['orog'] = cube_orog.units
        reference = cube_orog
        iris.save(cube_orog, os.path.join(self.outdir, 'orography.nc'))

        for var in self.variables_input:
            filename = self._file_pattern(self.time_range.split('-')[0],
                                          self.time_range.split('-')[1],
                                          variables_required[var]["fname"],
                                          um_stream='pt')

            input_path = os.path.join(self.input_directory, filename)
            if not os.path.exists(input_path):
                msg = f"Unable to find expected input file {input_path}"
                self.logger.error(msg)
                raise RuntimeError(msg)

            # make the output path filename similar to CMIP6 naming, will be standard
            # regardless of the input filename structure
            # varname_new, freq, time
            var_name = var
            if "varname_new" in variables_required[var]:
                var_name = variables_required[var]["varname_new"]
            output_path = self._file_pattern_processed(time_start,
                                                       time_end,
                                                       var_name,
                                                       self.frequency)

            output_path = os.path.join(self.outdir, output_path)
            # temporary for testing #
            if os.path.exists(output_path):
                source_filenames[var] = input_path
                processed_filenames[var] = output_path
                continue
            if not os.path.exists(os.path.dirname(output_path)):
                os.makedirs(os.path.dirname(output_path))

            # apart from slp, regrid the data to the slp grid (all input data for
            # TempestExtremes needs to be on the same grid)
            # regrid u and v to t grid and rename variable if necessary
            cube = iris.load_cube(input_path)
            cube = iris.util.squeeze(cube)
            ndims = cube.shape
            variable_units[var] = cube.units
            if var == 'uas':
                variable_units['sfcWind'] = cube.units
            self.logger.debug(f"regrid file {input_path}")
            regridded = cube.regrid(reference, iris.analysis.Linear())
            if "varname_new" in variables_required[var]:
                regridded.var_name = variables_required[var]["varname_new"]

            output_paths = []
            if len(ndims) == 3:
                iris.save(regridded, output_path)
                #output_paths.append(output_path)
                output_paths = output_path
            else:
                level_coord = regridded.dim_coords[1]
                for ilevel, level in enumerate(level_coord.points):
                    regridded_level = regridded[:, ilevel, :, :]
                    regridded_level.var_name = regridded.var_name+'_'+str(int(level))
                    iris.save(regridded_level, output_path[:-3]+'_'+str(int(level))+'.nc')
                    output_paths.append(output_path[:-3]+'_'+str(int(level))+'.nc')
            print ('incoming paths ', [output_paths])
            check_paths = [output_paths] if isinstance(output_paths, str) else output_paths
            self._check_time_coord(check_paths)

            source_filenames[var] = input_path
            processed_filenames[var] = output_paths

        self.logger.debug(f"Orography file {processed_filenames['orog']}")

        return source_filenames, processed_filenames, variable_units

    def _identify_processed_files(self, time_start, time_end, grid_resol='native'):
        """
        Identify the processed input files to be used by tracking.
        The files have pseudo-CMIP6 filenames, using the processed variable names

        :param str time_start: The start time of this period of data YYYYMMDD
        :param str time_end: The end time of this period of data YYYYMMDD
        :param str grid_resol: The resolution string to be used if regridding
        :                      is required
        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """

        processed_filenames = {}
        variable_units = {}

        # Identify the grid and orography file
        processed_filenames["orog"] = os.path.join(
            self.outdir, "orography.nc"
        )
        cube = iris.load_cube(processed_filenames["orog"])
        variable_units['orog'] = cube.units
        longitude_size = cube.shape[-1]
        resolution = longitude_size // 2
        self.resolution_code = f"N{resolution}"

        if not os.path.exists(os.path.dirname(self.outdir)):
            raise Exception("Processed file directory does not exist, should come from pre-processing "\
                            +self.outdir)

        for var_name in self.variables_rename:
            # identify the processed path filename similar to CMIP6 naming, will be standard
            # regardless of the input filename structure
            # varname_new, freq, time,
            output_path = self._file_pattern_processed(time_start,
                                                       time_end,
                                                       var_name,
                                                       self.data_frequency)
            output_path = os.path.join(self.outdir, output_path)

            self.logger.debug(f"read file {output_path}")
            if os.path.exists(output_path):
                processed_filenames[var_name] = output_path
                cube = iris.load_cube(output_path)
                cube = iris.util.squeeze(cube)
                variable_units[var_name] = cube.units
                if var_name == 'uas':
                    variable_units['sfcWind'] = cube.units
                    variable_units['wind'] = cube.units
            else:
                raise Exception("Processed file does not exist "+output_path)

            print ('incoming path ', [output_path])
            self._check_time_coord([output_path])

        self.logger.debug(f"Orography file {processed_filenames['orog']}")

        return processed_filenames, variable_units

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

    def _rewrite_track_file(
        self,
        tracked_file_Tm1,
        tracked_file_T,
        tracked_file_Tm1_adjust,
        tracked_file_T_adjust,
        storms_match,
        track_type
    ):
        """
        Rewrite the .txt track files, removing the matching storms from the
        previous timestep which have been found in the current timestep and
        adding them to this current timestep
        :param str tracked_file_Tm1: The path to the track file from the previous
        :                            timestep.
        :param str tracked_file_T: The path to the track file from the current
        :                          timestep.
        :param str tracked_file_Tm1_adjust: The path to the updated track file
            for the previous time for output.
        :param str tracked_file_T_adjust: The path to the updated track file
            for the current timestep for output.
        :param list storms_match: The storms which have been found to match
            with a later time
        :param str track_type: the Tempest step from which the storm contents comes
            either stitch or profile
        """
        header_delim = 'start'

        if len(storms_match) == 0:
            copyfile(tracked_file_Tm1, tracked_file_Tm1_adjust)
            return

        with open(tracked_file_Tm1) as file_input:
            with open(tracked_file_Tm1_adjust, 'w') as file_output:
                for line in file_input:
                    line_array = line.split()
                    if header_delim in line:
                        line_of_traj = 0  # reset trajectory line to zero
                        matching_track = False
                        line_header = line
                        track_length = int(line_array[1])
                        start_date = line_array[2] + line_array[3].zfill(2) + \
                            line_array[4].zfill(2) + line_array[5].zfill(2)
                    else:
                        if line_of_traj <= track_length:
                            lon = float(line_array[2])
                            lat = float(line_array[3])
                            if line_of_traj == 0:
                                for storm in storms_match:
                                    storm_old = storm["early"]
                                    storm_new = storm["late"]
                                    date = self._storm_dates(storm_old)[0]
                                    if date == start_date and \
                                            track_length == storm_old["length"]:
                                        if lon == storm_old["lon"][0] and \
                                                lat == storm_old["lat"][0]:
                                            matching_track = True
                                if not matching_track:
                                    file_output.write(line_header)
                                    file_output.write(line)
                            else:
                                if not matching_track:
                                    file_output.write(line)
                            line_of_traj += 1

        if self._construct_command(track_type)["profile"] is not None:
            column_names = self.column_names[track_type+'_profile']
        else:
            column_names = self.column_names[track_type+'_stitch']

        with open(tracked_file_T) as file_input:
            with open(tracked_file_T_adjust, 'w') as file_output:
                for line in file_input:
                    line_array = line.split()
                    if header_delim in line:
                        line_of_traj = 0  # reset trajectory line to zero
                        matching_track = False
                        line_header = line
                        track_length = int(line_array[1])
                        start_date = line_array[2] + line_array[3].zfill(2) + \
                            line_array[4].zfill(2) + line_array[5].zfill(2)
                    else:
                        if line_of_traj <= track_length:
                            lon = float(line_array[2])
                            lat = float(line_array[3])
                            if line_of_traj == 0:
                                match_type = ''
                                for storm in storms_match:
                                    storm_old = storm["early"]
                                    storm_new = storm["late"]
                                    date = self._storm_dates(storm_new)[0]
                                    if date == start_date and \
                                            track_length == storm_new["length"]:
                                        if lon == storm_new["lon"][0] and \
                                                lat == storm_new["lat"][0]:
                                            matching_track = True
                                            match_type = storm["method"]
                                            storm_old_match = storm_old
                                            storm_new_match = storm_new
                                            match_offset = storm["offset"]
                                if not matching_track:
                                    file_output.write(line_header)
                                    file_output.write(line)
                                else:
                                    if match_type == 'extend':
                                        line_extra = ''
                                        new_length = track_length + match_offset
                                        new_date_line, new_track_lines = \
                                            write_track_line(storm_old_match,
                                                             match_offset,
                                                             new_length,
                                                             column_names)
                                        line_header = new_date_line

                                        for new_line in new_track_lines:
                                            line_extra += new_line

                                    elif match_type == 'remove':
                                        line_extra = ''
                                    else:
                                        line_extra = ''
                                    file_output.write(line_header)
                                    file_output.write(line_extra)
                                    file_output.write(line)

                            else:
                                file_output.write(line)
                            line_of_traj += 1

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

        self.input_directory = self.app_config.get_property("common",
                                                            "input_directory")
        self.output_directory = self.app_config.get_property(
            "common", "output_directory"
        )
        self.tc_detect_script = self.app_config.get_property(
            "common", "tc_detect_script"
        )
        self.tc_stitch_script = self.app_config.get_property(
            "common", "tc_stitch_script"
        )
        self.tc_editor_script = self.app_config.get_property(
            "common", "tc_editor_script"
        )
        self.orography_dir = self.app_config.get_property("common", "orography_dir")
        self.delete_processed = self.app_config.get_bool_property(
            "common", "delete_processed"
        )
        self.delete_source = self.app_config.get_bool_property(
            "common", "delete_source"
        )
        self.plot_tracks = self.app_config.get_bool_property("common", "plot_tracks")
        # track_types is a Python list so eval converts str to list
        self.track_types = eval(self.app_config.get_property("common", "track_types"))
        self.variables_input = eval(self.app_config.get_property("common",
                                                                 "variables_input"))
        self.variables_rename = eval(self.app_config.get_property("common",
                                                                 "variables_rename"))
        self.nodeedit_vars = eval(self.app_config.get_property("common",
                                                               "nodeedit_vars"))
        self.um_file_pattern = self.app_config.get_property("common",
                                                            "um_file_pattern")
        self.file_pattern_processed = self.app_config.get_property("common",
                                                            "file_pattern_processed")
        self.regrid_resolutions = \
            eval(self.app_config.get_property("common", "regrid_resolutions"))
        self.data_frequency = self.app_config.get_property("common",
                                                            "data_frequency")
        self.outputcmd_detect_default = self.app_config.get_property("common", "outputcmd_detect_default")
        self.in_fmt_stitch_default = self.app_config.get_property("common", "in_fmt_stitch_default")
        self.out_fmt_profile1_default = self.app_config.get_property("common", "out_fmt_profile1_default")
        self.out_fmt_profile2_default = self.app_config.get_property("common", "out_fmt_profile2_default")

    def _get_environment_variables(self):
        """
        Get the required environment variables from the suite. A list and
        explanation of the required environment variables is included in the
        documentation.
        """
        try:
            self.um_runid = os.environ["RUNID_OVERRIDE"]
        except:
            self.um_runid = os.environ["RUNID"]
        try:
            self.um_suiteid = os.environ["SUITEID_OVERRIDE"]
        except:
            self.um_suiteid = os.environ["CYLC_SUITE_NAME"]
        self.cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
        self.time_cycle = os.environ["TIME_CYCLE"]
        self.previous_cycle = os.environ["PREVIOUS_CYCLE"]
        self.tm2_cycle = os.environ["TM2_CYCLE"]
        self.next_cycle = os.environ["NEXT_CYCLE"]
        self.startdate = os.environ["STARTDATE"]
        self.enddate = os.environ["ENDDATE"]
        self.lastcycle = os.environ["LASTCYCLE"]
        self.is_last_cycle = os.environ["IS_LAST_CYCLE"]
        self.ncodir = os.environ["NCODIR"]
