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
        self.input_files = {}
        self.variable_units = {}
        self.calendar_units = 'days since 1950-01-01 00:00:00'
        self.regrid_resolutions = ['N96']
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
        self.logger.debug(f"time_stamp_day {timestamp_day}")
        self.logger.debug(f"enddate {self.enddate}")
        self.logger.debug(f"lastcycle {self.lastcycle}")
        self.logger.debug(f"is_last_cycle {self.is_last_cycle}")

        dot_file = 'do_tracking'
        self._write_dot_track_file(timestamp_day, dot_file=dot_file)

        dot_tracking_files = glob.glob(os.path.join(self.output_directory, dot_file+'*'))
        self.logger.debug(f"dot_tracking_files {dot_tracking_files}")
        if len(dot_tracking_files) > 0:
            for do_track_file in dot_tracking_files:
                ftimestamp_day = do_track_file.split('.')[1]
                # do not want to do calculations on data after the current cycle date
                if self._is_date_after(timestamp_day, ftimestamp_day):
                    continue
                fname = self._file_pattern(ftimestamp_day+'*', '*', 'slp', um_stream = 'pt')
                self.logger.debug(f"fname {fname}")
                file_search = os.path.join(
                    self.input_directory, fname
                )
                self.logger.debug(f"file_search {file_search}")
                if len(glob.glob(file_search)) > 0:
                    self.outdir = self.output_directory
                    input_files, variable_units = self._generate_data_files(ftimestamp_day)
                    self.input_files[ftimestamp_day] = input_files
                    self.variable_units = variable_units
                    self._run_steps(ftimestamp_day)
                    if self.delete_input:
                        self._tidy_data_files(ftimestamp_day)

                    if self.regrid_resolutions is not None:
                        for regrid_resol in self.regrid_resolutions:
                            self.outdir = self.output_directory+'_'+regrid_resol
                            input_files, variable_units = self._generate_data_files(ftimestamp_day, regrid_resol = regrid_resol)
                            self.input_files[ftimestamp_day] = input_files
                            self.variable_units = variable_units
                            self._run_steps(ftimestamp_day)
                            if self.delete_input:
                                self._tidy_data_files(ftimestamp_day)

                    # if this timestep has worked OK, then need to remove the dot_file and the data
                else:
                    self.logger.error(f"no files to process for timestamp " f"{ftimestamp_day}")
        else:
            self.logger.error(f"no dot files to process ")
 
    def _run_steps(self, timestamp):
        # Run the detection and stitching for this time period
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
                    if self.plot_tracks:
                        self.logger.debug(f"Plotting {os.path.basename(tracked_file)}")

                        self._read_and_plot_tracks(
                            tracked_file,
                            self.input_files[timestamp]["slp"],
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
        if self._is_new_year(timestamp):
            annual_tracks = self._run_annual_stitch(timestamp)
            for index, track_type in enumerate(self.track_types):
                self.cmd_detect = self.cmd_detect_type[track_type]
                self.cmd_stitch = self.cmd_stitch_type[track_type]
                self.cmd_edit = self.cmd_edit_type[track_type]

                self._read_and_plot_tracks(
                    annual_tracks[index],
                    self.input_files[timestamp]["slp"],
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

                self._read_and_plot_tracks(
                    wholerun_tracks[index],
                    self.input_files[timestamp]["slp"],
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

        :returns: The path to the candidate files (as a list ordered by track
            type) and details of the processed input files (as a dict).
        :rtype: tuple
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
            for key in self.input_files[timestamp]:
                fnames.append(self.input_files[timestamp][key])
            fnames.remove(self.input_files[timestamp]['orog'])
                
            in_file_list = os.path.join(self.outdir, 'in_file_list_detect.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(fnames) + ';'+self.input_files[timestamp]['orog']
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

        :returns: The string paths to the annual stitched file in a list.
        :rtype: list
        """

        tracked_files = []
        for track_type in self.track_types:
            previous_year = int(timestamp[:4]) - 1
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
        :param dict input_files: Details of the processed input files.
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
            #cmd_io += '--in_data "'+';'.join(self.input_files[timestamp][r] for r in self.input_files[timestamp])+'" '
            in_file_list = os.path.join(self.outdir, 'in_file_list_edit.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(self.input_files[timestamp][r] for r in self.input_files[timestamp])
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

    def _is_new_year(self, timestamp):
        """
        Check if this is the first cycle in a new year.

        :returns: True if this is the first submission period in a year.
        :rtype: bool
        """

        if timestamp[:4] != self.previous_cycle[:4]:
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

    def _write_dot_track_file(self, timestamp, dot_file = 'do_tracking'):
        """
        Write a file indicating that this timestep needs to be tracked

        """
        do_tracking_file = os.path.join(self.output_directory, dot_file+'.'+timestamp)
        if not os.path.exists(do_tracking_file):
            os.system('touch '+do_tracking_file)

    def _remove_dot_track_file(self, timestamp, dot_file = 'do_tracking'):
        """
        Remove a file indicating that this timestep needs to be tracked

        """
        do_tracking_file = os.path.join(self.output_directory, dot_file+'.'+timestamp)
        if os.path.exists(do_tracking_file):
            os.system('rm '+do_tracking_file)

    def _tidy_data_files(self, timestamp, dot_file = 'do_tracking'):
        """
        Remove input files and tracking dot file for this timestamp (tidy up)

        """
        self.logger.info(f"Tidy up input files")
        input_files_remove = self.input_files[timestamp]
        input_files_remove.pop('orog', None)
        for f in input_files_remove:
            cmd_rm = 'rm '+input_files_remove[f]
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
            
        self._remove_dot_track_file(timestamp, dot_file = dot_file)
        self.logger.debug(f"removed dot file {timestamp}")

    def _file_pattern(self, timestart, timeend, varname, um_stream = 'pt', frequency = '6h'):
        """
        Derive the input nc filenames from the file pattern 
        can be um as here, or could be other patterns for other models/platforms

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

        :param str track_type: The name of the type of tracking to run.
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

        input_files, variable_units = self._process_input_files(regrid_resol = regrid_resol)

        return input_files, variable_units

    def _process_input_files(self, regrid_resol = None):
        """
        Identify and then fix the grids and var_names in the input files.

        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """
        filetypes_required = self.variables_input
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

            processed_filenames[filetype] = output_path

        self.logger.debug(f"Orography file {processed_filenames['orog']}")

        return processed_filenames, variable_units

    def _read_and_plot_tracks(
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
        Read and then plot the tracks. The title is the specified prefix, the
        number of tracks (if requested) and then the suffix.

        :param str tracked_file: The path to the track file.
        :param str nc_file: The path to an input data netCDF file, which is
            used to gather additional information about the dates and calendars
            used in the data.
        :param str title_prefix: The title for the plot.
        :param str title_suffix: The title for the plot.
        :param bool include_num_tracks: Include the number of tracks loaded?
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
            self._create_netcdf(self.outdir, nc_file_out, storms, calendar, calendar_units, variable_units,
                                    startperiod=timestart, endperiod=timeend)

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
        self.delete_input = self.app_config.get_bool_property(
            "common", "delete_input"
        )
        self.plot_tracks = self.app_config.get_bool_property("common", "plot_tracks")
        # track_types is a Python list so eval converts str to list
        self.track_types = eval(self.app_config.get_property("common", "track_types"))
        self.variables_input = eval(self.app_config.get_property("common", "variables_input"))
        self.output_vars_default = eval(self.app_config.get_property("common", "output_vars_default"))
        self.output_vars_extra = eval(self.app_config.get_property("common", "output_vars_extra"))
        self.um_file_pattern = self.app_config.get_property("common", "um_file_pattern")

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
        self.startdate = os.environ["STARTDATE"]
        self.enddate = os.environ["ENDDATE"]
        self.lastcycle = os.environ["LASTCYCLE"]
        self.is_last_cycle = os.environ["IS_LAST_CYCLE"]

    def _define_netcdf_metadata(self, var, variable_units):
        long_name = 'unknown'
        description = 'unknown'
        units = '1'

        if 'slp' in var:
            standard_name = 'air_pressure_at_mean_sea_level'
            long_name = 'Sea Level Pressure'
            description = 'Sea level pressure for tracked variable'
            units = variable_units['slp']
        elif 'sfcWind' in var:
            standard_name = 'wind_speed'
            long_name = 'Near-surface Wind Speed'
            description = 'near-surface (usually 10 metres) wind speed'
            units = variable_units['sfcWind']
        elif 'orog' in var:
            standard_name = 'surface_altitude'
            long_name = 'Surface Altitude'
            description = 'Surface altitude (height above sea level)'
            units = variable_units['orog']
        elif 'wind' in var:
            standard_name = 'wind_speed'
            units = variable_units['sfcWind']
        elif 'rv' in var:
            standard_name = 'relative_vorticity'
            units = 's-1'
        elif 'zg' in var:
            standard_name = 'geopotential_height'
            long_name = 'Geopotential Height'
            description = 'Geopotential height difference'
            units = 'm'

        return standard_name, long_name, description, units

    def _create_netcdf(self, directory, savefname, storms, calendar, time_units, variable_units,
                      startperiod=None, endperiod=None):
        """
        Create netcdf file for the tracks. 
        May need metadata from a model nc file, so may need to create at a time when these are available
        """
        print ('making netCDF of outputs')
        
        self.savefname = savefname
        nc = Dataset(self.savefname, 'w', format='NETCDF4')
        nc.title = 'Tempest TC tracks'
        nc.directory = directory
        nc.tracked_data_frequency = self.frequency
        #nc.TRACK_DURATION_MIN = np.float64(self.TRACK_DURATION_MIN)

        nc.mo_runid = self.um_suiteid
        nc.grid = self.resolution_code
        nc.start_date = startperiod
        nc.end_date = endperiod
        nc.institution_id = 'MOHC'
        nc.algorithm = 'TempestExtremes_v2'
        nc.algorithm_ref = 'Ullrich and Zarzycki 2017; Zarzycki and Ullrich 2017; Ullrich et al. 2020'
        nc.detect_cmd = self.cmd_detect
        nc.stitch_cmd = self.cmd_stitch

        record_length = 0; tracks = 0
        for storm in storms:
            tracks += 1
            storm_length = storm['length']
            record_length += storm_length

        nc.createDimension('tracks', size = tracks) # unlimited
        nc.createDimension('record', size = record_length)

        nc.createVariable('FIRST_PT', np.int32, ('tracks'))   
        nc.createVariable('NUM_PTS', np.int32, ('tracks'))  
        nc.createVariable('TRACK_ID', np.int32, ('tracks'))  
        nc.createVariable('index', np.int32, ('record'))  
        nc.createVariable('time', 'f8', ('record'))
        nc.createVariable('lon', 'f4', ('record'))
        nc.createVariable('lat', 'f4', ('record'))

        if self.output_vars_extra != None:
            output_vars_all = self.output_vars_default.copy()
            output_vars_all.extend(self.output_vars_extra)

        for var in output_vars_all:
            nc.createVariable(var, 'f8', ('record'))

        nc.variables['FIRST_PT'].units = 'ordinal'
        nc.variables['FIRST_PT'].long_name = 'first_pt'
        nc.variables['FIRST_PT'].description = 'Index to first point of this track number'

        nc.variables['NUM_PTS'].units = 'ordinal'
        nc.variables['NUM_PTS'].long_name = 'num_pts'
        nc.variables['NUM_PTS'].description = 'Number of points for this track'

        nc.variables['TRACK_ID'].units = 'ordinal'
        nc.variables['TRACK_ID'].long_name = 'track_id'
        nc.variables['TRACK_ID'].description = 'Tropical cyclone track number'

        nc.variables['index'].units = 'ordinal'
        nc.variables['index'].long_name = 'track_id'
        nc.variables['index'].description = 'Track sequence number (0 - length of track - 1)'

        nc.variables['lat'].units = 'degrees_north'
        nc.variables['lat'].standard_name = 'latitude'
        nc.variables['lat'].long_name = 'latitude'
        nc.variables['lat'].description = 'Latitude (degrees north) associated with tracked variable'

        nc.variables['lon'].units = 'degrees_east'
        nc.variables['lon'].standard_name = 'longitude'
        nc.variables['lon'].long_name = 'longitude'
        nc.variables['lon'].description = 'Longitude (degrees east) associated with tracked variable'

        nc.variables['time'].units = time_units
        nc.variables['time'].calendar = calendar
        nc.variables['time'].standard_name = 'time'
        nc.variables['time'].long_name = 'time'

        for var in output_vars_all:
            standard_name, long_name, description, v_units = self._define_netcdf_metadata(var, variable_units)
            self.logger.debug(f"var, units {var} {v_units} ")
            nc.variables[var].standard_name = standard_name
            nc.variables[var].long_name = long_name
            nc.variables[var].description = description
            nc.variables[var].units = str(v_units)

        # read the storms and write the values to the file
        # track: first_pt, num_pts, track_id
        # record: lat, lon, time, slp, index(0:tracklen-1)
        first_pt = []; num_pts = []; track_id = []
        lon = []; lat = []; time = []; index = []

        variables_to_write = {}
        for var in output_vars_all:
            variables_to_write[var] = []

        first_pt_index = 0
        for ist, storm in enumerate(storms):
            first_pt.append(first_pt_index)
            num_pts.append(storm['length'])
            track_id.append(ist)
            first_pt_index += storm['length']

            for ipt in range(storm['length']):
                tunit = utime(time_units, calendar = calendar)
                t1 = tunit.date2num(datetime(storm['year'][ipt], storm['month'][ipt], storm['day'][ipt], storm['hour'][ipt]))
                time.append(t1)
                index.append(ipt)
                lon.append(storm['lon'][ipt])
                lat.append(storm['lat'][ipt])
                for var in output_vars_all:
                    variables_to_write[var].append(storm[var][ipt])
        
        self.logger.debug(f"first_pt {first_pt} ")
        self.logger.debug(f"tracks, record_length {tracks} {record_length} ")
        self.logger.debug(f"len(first_pt) {len(first_pt)} ")
        self.logger.debug(f"len(lon) {len(lon)} ")
        self.logger.debug(f"len(variables_to_write(slp)) {len(variables_to_write['slp'])} ")
        #self.logger.debug(f"variables_to_write[slp] {variables_to_write['slp']} ")
        # now write variables to netcdf 
        nc.variables['FIRST_PT'][:] = first_pt
        nc.variables['NUM_PTS'][:] = num_pts
        nc.variables['TRACK_ID'][:] = track_id
        nc.variables['index'][:] = index
        nc.variables['lon'][:] = lon
        nc.variables['lat'][:] = lat
        nc.variables['time'][:] = time
        for var in output_vars_all:
            self.logger.debug(f"var {var} ")
            nc.variables[var][:] = variables_to_write[var]

        nc.close()

    def _write2netcdf(self, rtime, stopper=0):
        """
        Write tracks to netcdf file - for the original ocean eddy code, the inactive (finished) eddies were written, not quite sure the equivalent here.
        'ncind' is important because prevents writing of
        already written tracks.
        Each inactive track is 'emptied' after saving
        
        rtime - current timestamp
        stopper - dummy value (either 0 or 1)
        """
        print ('Writing to netCDF')
        
        rtime += stopper
        tracks2save = np.array([self.get_inactive_tracks(rtime)])
        HBR = self.HOURS_BTWN_RECORDS

        if tracks2save.any(): 

            with Dataset(self.savedir, 'a') as nc:

                for i in np.nditer(tracks2save):

                    # saved2nc is a flag indicating if track[i] has been saved
                    if (not self.tracklist[i].saved2nc) and \
                       (np.all(self.tracklist[i].ocean_time)):

                        tsize = len(self.tracklist[i].lon)

                        if (tsize >= self.TRACK_DURATION_MIN / HBR) and tsize >= 1.:
                            lon = np.array([self.tracklist[i].lon])
                            lat = np.array([self.tracklist[i].lat])
                            slp = np.array([self.tracklist[i].slp])
                            wind10m = np.array([self.tracklist[i].wind10m])
                            track = np.full(tsize, self.ch_index, dtype=np.int32)
                            num_pts = np.arange(tsize, dtype=np.int32)
                            first_pt = np.arange(tsize, dtype=np.int32)
                            track_id = np.arange(tsize, dtype=np.int32)
                            time = np.array([self.tracklist[i].time])
                            #time = num2date(time, units = \
                            #    nc.variables['time'].units, \
                            #    calendar = nc.variables['time'].calendar)


                            tend = self.ncind + tsize
                            nc.variables['lon'][self.ncind:tend] = lon
                            nc.variables['lat'][self.ncind:tend] = lat
                            nc.variables['slp'][self.ncind:tend] = lp
                            nc.variables['time'][self.ncind:tend] = time

                            self.tracklist[i].saved2nc = True
                            self.ncind += tsize
                            self.ch_index += 1
                            nc.sync()

 
