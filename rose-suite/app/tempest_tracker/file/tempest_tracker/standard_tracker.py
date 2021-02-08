# (C) British Crown Copyright 2020, Met Office.
# Please see LICENSE for license details.
import glob
import os
import re
import shutil
import subprocess
import sys

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
        self.variable_units = {}
        self.calendar_units = 'days since 1950-01-01 00:00:00'
        self.regrid_resolutions = None
        self.outdir = None

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

        if not os.path.exists(self.output_directory):
            try:
                os.makedirs(self.output_directory)
            except PermissionError:
                msg = f"Unable to create output directory {self.output_directory}"
                self.logger.error(msg)
                sys.exit(1)

        self.logger.debug(
            f"CYLC_TASK_CYCLE_TIME {self.cylc_task_cycle_time}, "
            f"um_runid {self.um_runid}"
        )

        # First write a dot_file to document which timestamps are yet to be tracked
        # Check for input files available
        # for um postproc nc files, these take form 
        timestamp_day = self.cylc_task_cycle_time[:8]
        timestamp_endday = self.next_cycle[:8]
        self.logger.debug(f"time_stamp_day {timestamp_day}")
        self.logger.debug(f"enddate {self.enddate}")
        self.logger.debug(f"lastcycle {self.lastcycle}")
        self.logger.debug(f"previous_cycle {self.previous_cycle}")
        self.logger.debug(f"next_cycle {self.next_cycle}")
        self.logger.debug(f"is_last_cycle {self.is_last_cycle}")

        dot_file = 'do_tracking'
        self._write_dot_track_file(timestamp_day, timestamp_endday, dot_file=dot_file)

        dot_tracking_files = sorted(glob.glob(os.path.join(self.output_directory, dot_file+'*')))
        self.logger.debug(f"dot_tracking_files {dot_tracking_files}")
        if dot_tracking_files:
            for do_track_file in dot_tracking_files:
                ftimestamp_day = do_track_file.split('.')[1].split('-')[0]
                ftimestamp_endday = do_track_file.split('.')[1].split('-')[1]
                # do not want to do calculations on data after the current cycle date
                if self._is_date_after(timestamp_day, ftimestamp_day):
                    continue
                fname = self._file_pattern(ftimestamp_day+'*', '*', 'slp', um_stream = 'pt')
                self.logger.debug(f"fname {fname}")
                file_search = os.path.join(self.input_directory, fname)
                if glob.glob(file_search):
                    self.outdir = self.output_directory
                    source_files, processed_files, variable_units = self._generate_data_files(ftimestamp_day)
                    self.processed_files[ftimestamp_day] = processed_files
                    self.variable_units = variable_units
                    self._run_steps(ftimestamp_day, ftimestamp_endday)
                    if self.delete_processed:
                        self._tidy_data_files(ftimestamp_day, ftimestamp_endday)

                    if self.regrid_resolutions is not None:
                        for regrid_resol in self.regrid_resolutions:
                            self.outdir = self.output_directory+'_'+regrid_resol
                            source_files, processed_files, variable_units = self._generate_data_files(ftimestamp_day, regrid_resol = regrid_resol)
                            self.source_files[ftimestamp_day] = source_files
                            self.processed_files[ftimestamp_day] = processed_files
                            self.variable_units = variable_units
                            self._run_steps(ftimestamp_day, ftimestamp_endday)
                            if self.delete_processed:
                                self._tidy_data_files(ftimestamp_day, ftimestamp_endday)

                    if self.delete_source:
                        self._tidy_data_files(ftimestamp_day, ftimestamp_endday, f_remove = 'source')

                    # if this timestep has worked OK, then need to remove the dot_file and the data
                else:
                    self.logger.error(f"no files to process for timestamp " f"{ftimestamp_day}")
        else:
            self.logger.error(f"no dot files to process ")
 
    def _run_steps(self, timestamp, timestamp_end):
        # Run the detection, stitching and editing for this time period
        candidate_files = self._run_detection(timestamp)
        tracked_files = self._run_stitching(candidate_files)
        edited_files = self._run_node_file_editor(tracked_files, timestamp)

        # Run the plotting for this time period (if required and
        # data available)
        for index, track_type in enumerate(self.track_types):
            candidate_file = candidate_files[index]
            tracked_file = tracked_files[index]
            self.cmd_detect = self.cmd_detect_type[track_type]
            self.cmd_stitch = self.cmd_stitch_type[track_type]
            self.cmd_edit = self.cmd_edit_type[track_type]

            if os.stat(candidate_file).st_size > 0:
                if os.stat(tracked_file).st_size > 0:
                    self.logger.debug(f"Plotting {os.path.basename(tracked_file)}")

                    self._read_write_and_plot_tracks(
                        tracked_file,
                        self.processed_files[timestamp]["slp"],
                        self.time_range.split('-')[0],
                        self.time_range.split('-')[1],
                        self.variable_units,
                        title_prefix=f"{self.um_runid} {self.resolution_code} "
                        f"{self.time_range}",
                        title_suffix=f"{track_type} tracks",
                        write_to_netcdf = True
                        )
                else:
                    self.logger.error(
                        f"candidate file has data but there are "
                        f"no tracks in {tracked_file}"
                    )
            else:
                self.logger.error(f"candidate file is empty " f"{candidate_file}")

        # Produce an annual summary if this is the first period in a year
        if self._is_new_year(timestamp, timestamp_end):
            annual_tracks = self._run_annual_stitch(timestamp)
            for index, track_type in enumerate(self.track_types):
                self.cmd_detect = self.cmd_detect_type[track_type]
                self.cmd_stitch = self.cmd_stitch_type[track_type]
                self.cmd_edit = self.cmd_edit_type[track_type]

                self._read_write_and_plot_tracks(
                    annual_tracks[index],
                    self.processed_files[timestamp]["slp"],
                    timestamp[:4]+'0101',
                    timestamp[:4]+'1231',
                    self.variable_units,
                    title_prefix=f"{self.um_runid} {self.resolution_code} "
                    f"{timestamp[:4]}",
                    title_suffix=f"{track_type} annual tracks",
                    write_to_netcdf = True
                )

            # plot annual tracks (also fills in gaps in tracks)
            # create netcdf output file for tracks
        else:
            self.logger.debug("Not running annual stitch")

        # Produce an wholerun summary if this is the end of the run
        self.logger.debug(f"Running wholerun {self.lastcycle} {timestamp} {self.is_last_cycle}")
        if self.is_last_cycle == 'true':
            wholerun_tracks = self._run_wholerun_stitch()
            for index, track_type in enumerate(self.track_types):
                self.cmd_detect = self.cmd_detect_type[track_type]
                self.cmd_stitch = self.cmd_stitch_type[track_type]
                self.cmd_edit = self.cmd_edit_type[track_type]

                self._read_write_and_plot_tracks(
                    wholerun_tracks[index],
                    self.processed_files[timestamp]["slp"],
                    self.startdate[:8],
                    self.enddate[:8],
                    self.variable_units,
                    title_prefix=f"{self.um_runid} {self.resolution_code} "
                    f"{self.startdate[:8]} - {self.enddate[:8]}",
                    title_suffix=f"{track_type} wholerun tracks",
                    write_to_netcdf = True
                )
        else:
            self.logger.debug("Not running wholerun stitch")

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

            # need the first input file not to be orography, since that file has to have a time coordinate
            fnames = []
            for key in self.processed_files[timestamp]:
                fnames.append(self.processed_files[timestamp][key])
            fnames.remove(self.processed_files[timestamp]['orog'])
                
            in_file_list = os.path.join(self.outdir, 'in_file_list_detect.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(fnames) + ';'+self.processed_files[timestamp]['orog']
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
                    f"EXCEPTION found in TempestExtreme detect output\n" f"{sts.stdout}"
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
                self.outdir, f"{self.um_runid}_track_{self.time_range}_{track_type}.txt"
            )
            self._stitch_file(candidatefile, trackedfile, track_type)
            tracked_files.append(trackedfile)

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
                self.outdir, f"{self.um_runid}_track_year_{previous_year}_{track_type}.txt"
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
                f"{self.um_runid}_candidate_fullrun_{start_period}_{end_period}_{track_type}.txt",
            )
            with open(wholerun_candidate, "w") as out_file:
                for candidate_file in candidate_files:
                    with open(candidate_file) as in_file:
                        for line in in_file:
                            out_file.write(line)

            wholerun_track = os.path.join(
                self.outdir, f"{self.um_runid}_track_fullrun_{start_period}_{end_period}_{track_type}.txt"
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

    def _run_node_file_editor(self, tracked_files, timestamp):
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

        for index, track_type in enumerate(self.track_types):
            tracking_phase_commands = self._construct_command(track_type)
            tracked_file = tracked_files[index]
            edited_file = tracked_file[:-4] + "_edited_profile.txt"

            cmd_io = '{} --in_nodefile {} --out_nodefile {} '.format(
                    self.tc_editor_script,
                    tracked_file,
                    edited_file)
            in_file_list = os.path.join(self.outdir, 'in_file_list_edit.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(self.processed_files[timestamp][r] for r in self.processed_files[timestamp])
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)
            cmd_io += '--in_data_list '+in_file_list+' '

            cmd_edit = cmd_io + tracking_phase_commands["profile"]
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
                    f"EXCEPTION found in TempestExtreme editor output\n" f"{sts.stdout}"
                )
                raise RuntimeError(msg)
            edited_files.append(edited_file)

        return edited_files

    def _is_new_year(self, timestamp, timestamp_end):
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

    def _is_date_after(self, timeref, timetest):
        """
        Check if timetest is after timestamp.

        :returns: True if timetest is after timeref.
        :rtype: bool
        """

        if int(timetest[:8]) > int(timeref[:8]):
            return True
        else:
            return False

    def _write_dot_track_file(self, timestamp, timestamp_end, dot_file = 'do_tracking'):
        """
        Write a file indicating that this timestep needs to be tracked

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param str dot_file: The first part of the string of a filename to indicate which time periods still need tracking

        """
        do_tracking_file = os.path.join(self.output_directory, dot_file+'.'+timestamp+'-'+timestamp_end)
        if not os.path.exists(do_tracking_file):
            os.system('touch '+do_tracking_file)

    def _remove_dot_track_file(self, timestamp, timestamp_end, dot_file = 'do_tracking'):
        """
        Remove a file indicating that this timestep needs to be tracked

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param str dot_file: The first part of the string of a filename to indicate which time periods still need tracking

        """
        do_tracking_file = os.path.join(self.output_directory, dot_file+'.'+timestamp+'-'+timestamp_end)
        if os.path.exists(do_tracking_file):
            os.system('rm '+do_tracking_file)

    def _tidy_data_files(self, timestamp, timestamp_end, dot_file = 'do_tracking', f_remove = 'processed'):
        """
        Remove input files and tracking dot file for this timestamp (tidy up)

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param str dot_file: The first part of the string of a filename to indicate which time periods still need tracking
        :param str f_remove: An indicator of which files need to be deleted. 
                   processed = the (regridded) files read by the tracking code
                   source    = the files produced by the model
        """
        self.logger.info(f"Tidy up input files")
        files_remove = []
        if f_remove == 'processed':
            files_remove = self.processed_files[timestamp]
            files_remove.pop('orog', None)
        elif f_remove == 'source':
            if self.delete_source:
                files_remove = self.source_files[timestamp]
                files_remove.pop('orog', None)

        for f in files_remove:
            cmd_rm = 'rm '+files_remove[f]
            self.logger.info(f"cmd {cmd_rm}")
            sts = subprocess.run(
                cmd_rm,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
            self.logger.debug(sts.stdout)
            if "EXCEPTION" in sts.stdout:
                msg = (
                    f"EXCEPTION found in tidying up input data \n" f"{sts.stdout}"
                )
                raise RuntimeError(msg)
            
        if f_remove == 'processed':
            self._remove_dot_track_file(timestamp, timestamp_end, dot_file = dot_file)
        self.logger.debug(f"removed dot file {timestamp}")

    def _file_pattern(self, timestart, timeend, varname, um_stream = 'pt', frequency = '6h'):
        """
        Derive the input nc filenames from the file pattern, assuming a 
        um model filenaming pattern as here, or could be other patterns 
        for other models/platforms (which would need to be added)

        :param str timestart: The timestep of the start of the data period to process
        :param str timeend: The timestep of the end of the data period to process
        :param str um_stream: The name of the um output stream (output file identification)
        :param str frequency: The frequency of the input data (in hours, needs to include "h"), used to determine file naming
        :returns: a filename given the inputs to the pattern
        :rtype: str
        """
        if self.frequency == None:
            file_freq = frequency
        else:
            file_freq = str(self.frequency)+'h'

        if self.um_file_pattern != '':
            fname = self.um_file_pattern.format(
                self.um_runid,
                file_freq,
                timestart, 
                timeend,
                um_stream,
                varname
        )
        return fname
            

    def _construct_command(self, track_type):
        """
        Read the TempestExtreme command line parameters from the configuration.

        :param str track_type: The name of the type of tracking to run, possible values:
               detect, stitch, profile
        :returns: A dictionary with keys of `detect` and `stitch` and the
            values for each of these is a string containing the command
            line parameters for each of these TempestExtreme steps. The
            parameters are sorted into alphabetical order in each line.
        :rtype: dict
        """

        commands = {}
        for step in ["detect", "stitch", "profile"]:
            step_config = self.app_config.section_to_dict(f"{track_type}_{step}")

            step_arguments = [
                f"--{parameter} {step_config[parameter]}"
                for parameter in sorted(list(step_config.keys()))
                if step_config[parameter]
            ]
            commands[step] = " ".join(step_arguments)
        return commands

    def _generate_data_files(self, timestamp, regrid_resol = None):
        """
        Identify and then fix the grids and var_names in the input files.
        The time_range and frequency attributes are set when this method runs.

        :param str timestamp: The timestep of the start of the data period to process
        :param str regrid_resol: The resolution string is regridding is required
        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """
        timestamp_day = timestamp
        self.logger.debug(f"time_stamp_day {timestamp_day}")
        fname = self._file_pattern(timestamp_day+'*', '*', 'slp', um_stream = 'pt')
        file_search = os.path.join(
            self.input_directory, fname
        )
        self.logger.debug(f"file_search {file_search}")
        files = sorted(glob.glob(file_search))
        self.logger.debug(f"files {files}")
        if not files:
            msg = f"No input files found for glob pattern {file_search}"
            self.logger.error(msg)
            raise RuntimeError(msg)

        time_ranges = []
        unique_ranges = []
        frequencies = []
        unique_frequencies = []
        for filename in files:
            time_range = os.path.basename(filename).split("_")[3]
            frequency = os.path.basename(filename).split("_")[2]
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

        source_files, processed_files, variable_units = self._process_input_files(regrid_resol = regrid_resol)

        return source_files, processed_files, variable_units

    def _process_input_files(self, regrid_resol = None):
        """
        Identify and then fix the grids and var_names in the input files.

        :param str regrid_resol: The resolution string to be used if regridding is required
        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """
        filetypes_required = self.variables_input
        source_filenames = {}
        processed_filenames = {}
        variable_units = {}

        variables_required = {}
        # these variables need to have a new var_name, either because the default from the UM is confusing or unknown, and these names are needed for the variable name inputs for the TempestExtremes scripts
        variables_rename = ['zg', 'rv', 'rvT63', 'uas', 'vas', 'viwve', 'viwvn', 'ta', 'rvT42']
        for var in self.variables_input:
            variables_required[var] = {'fname': var}
            if var in variables_rename:
                variables_required[var].update({'varname_new': var})

        reference_name = self._file_pattern(self.time_range.split('-')[0], self.time_range.split('-')[1], variables_required["slp"]["fname"], um_stream = 'pt')
        reference_path = os.path.join(self.input_directory, reference_name)
        reference = iris.load_cube(reference_path)
        variable_units['slp'] = reference.units

        # Identify the grid and orography file
        if regrid_resol == None:
            longitude_size = reference.shape[-1]
            resolution = longitude_size // 2
            self.resolution_code = f"N{resolution}"

            processed_filenames["orog"] = os.path.join(
                self.orography_dir, f"orog_HadGEM3-GC31-{self.resolution_code}e.nc"
            )
        else:
            resolution_code = f"N{regrid_resol}"
            processed_filenames["orog"] = os.path.join(
                self.orography_dir, f"orog_HadGEM3-GC31-{regrid_resol}e.nc"
            )
        cube_orog = iris.load_cube(processed_filenames["orog"])
        variable_units['orog'] = cube_orog.units
        reference = cube_orog

        for filetype in filetypes_required:
            filename = self._file_pattern(self.time_range.split('-')[0], self.time_range.split('-')[1], variables_required[filetype]["fname"], um_stream = 'pt')

            input_path = os.path.join(self.input_directory, filename)
            if not os.path.exists(input_path):
                msg = f"Unable to find expected input file {input_path}"
                self.logger.error(msg)
                raise RuntimeError(msg)

            output_path = os.path.join(self.outdir, filename)
            if not os.path.exists(os.path.dirname(output_path)):
                os.makedirs(os.path.dirname(output_path))

            # apart from slp, regrid the data to the slp grid (all input data for TempestExtremes needs to be on the same grid)
            if 'slp' in filetype:
                if not os.path.exists(output_path):
                    if regrid_resol == None:
                        shutil.copyfile(input_path, output_path)
                    else:
                        cube = iris.load_cube(input_path)
                        regridded = cube.regrid(cube_orog, iris.analysis.Linear())
                        iris.save(regridded, output_path)
                    
            else:
                # regrid u and v to t grid and rename variable if necessary
                cube = iris.load_cube(input_path)
                variable_units[filetype] = cube.units
                if filetype == 'uas':
                    variable_units['sfcWind'] = cube.units
                regridded = cube.regrid(reference, iris.analysis.Linear())
                if "varname_new" in variables_required[filetype]:
                    regridded.var_name = variables_required[filetype]["varname_new"]
                iris.save(regridded, output_path)

            source_filenames[filetype] = input_path
            processed_filenames[filetype] = output_path

        self.logger.debug(f"Orography file {processed_filenames['orog']}")

        return source_filenames, processed_filenames, variable_units

    def _read_write_and_plot_tracks(
        self,
        tracked_file,
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
        :param dict variable_units: The units of each output variable, as derived from the input data files
        :param str title_prefix: The title for the plot.
        :param str title_suffix: The title for the plot.
        :param bool include_num_tracks: Include the number of tracks within the plot title
        :param bool write_to_netcdf: Whether to write tracks to netcdf file or not
        """

        # storms returned as dictionary with keys: length, lon, lat, year, month, day, hour, step
        # if we define extra output variables, then need to define their coordinate position so that we can read them into the storm variable
        coords_new = {}
        if self.output_vars_extra != None:
            for ic, coord in enumerate(self.output_vars_default):
                coords_new[coord] = ic+4
            for ic1, coord1 in enumerate(self.output_vars_extra):
                coords_new[coord1] = ic1+len(self.output_vars_default)+4

        storms = get_trajectories(tracked_file, nc_file_in, self.frequency, coords_new = coords_new)
        self.logger.debug(f"got storms {len(storms)}")
        self.logger.debug(f"got storms {storms[0]}")

        # write the storms to a netcdf file
        if write_to_netcdf:
        # write out netcdf files 
            with Dataset(nc_file_in, 'r') as x:
                time = x.variables['time']
                calendar = time.calendar
                calendar_units = self.calendar_units
            nc_file_out = self.um_runid+'_'+tracked_file[:-4]+'.nc'
            nc_file_out = os.path.join(os.path.dirname(tracked_file), os.path.basename(tracked_file)[:-4]+'.nc')
            self.logger.debug(f"open netcdf file {nc_file_out}")
            save_trajectories_netcdf(
                self.outdir, 
                nc_file_out, 
                storms, 
                calendar, 
                calendar_units, 
                variable_units, 
                self.frequency, 
                self.um_suiteid, 
                self.resolution_code, 
                self.cmd_detect, 
                self.cmd_stitch, 
                self.output_vars_default, 
                output_vars_extra=self.output_vars_extra, 
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

    def _get_app_options(self):
        """Get commonly used configuration items from the config file"""

        self.input_directory = self.app_config.get_property("common", "input_directory")
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
        self.slp_std_name = self.app_config.get_property("common", "slp_std_name")
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
        self.variables_input = eval(self.app_config.get_property("common", "variables_input"))
        self.output_vars_default = eval(self.app_config.get_property("common", "output_vars_default"))
        self.output_vars_extra = eval(self.app_config.get_property("common", "output_vars_extra"))
        self.um_file_pattern = self.app_config.get_property("common", "um_file_pattern")
        self.regrid_resolutions = eval(self.app_config.get_property("common", "regrid_resolutions"))

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
            self.um_suiteid= os.environ["SUITEID_OVERRIDE"]
        except:
            self.um_suiteid= os.environ["CYLC_SUITE_NAME"]
        self.cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
        self.time_cycle = os.environ["TIME_CYCLE"]
        self.previous_cycle = os.environ["PREVIOUS_CYCLE"]
        self.next_cycle = os.environ["NEXT_CYCLE"]
        self.startdate = os.environ["STARTDATE"]
        self.enddate = os.environ["ENDDATE"]
        self.lastcycle = os.environ["LASTCYCLE"]
        self.is_last_cycle = os.environ["IS_LAST_CYCLE"]

