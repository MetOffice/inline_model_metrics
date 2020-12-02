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
import iris

from afterburner.apps import AbstractApp

from tempest_helper import (
    count_trajectories,
    get_trajectories,
    plot_trajectories_cartopy,
)


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
        self.FILLVAL = 9999

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

        dot_file = 'do_tracking'
        self._write_dot_track_file(timestamp_day, dot_file = dot_file)

        dot_tracking_files = glob.glob(os.path.join(self.output_directory, dot_file+'*'))
        self.logger.debug(f"dot_tracking_files {dot_tracking_files}")
        if len(dot_tracking_files) > 0:
            for do_track_file in dot_tracking_files:
                ftimestamp_day = do_track_file.split('.')[1]
                file_search = os.path.join(
                    self.input_directory, f"atmos_{self.um_runid}*{ftimestamp_day}??-*-slp.nc"
                )
                self.logger.debug(f"file_search {file_search}")
                if len(glob.glob(file_search)) > 0:
                    #filenames = self._generate_data_files(ftimestamp_day)
                    self._run_steps(ftimestamp_day)
                    # if this timestep has worked OK, then need to remove the dot_file
                    # self.remove_dot_track_file(ftimestamp_day, dot_file = dot_file)
                else:
                    self.logger.error(f"no files to process for timestamp " f"{ftimestamp_day}")
        else:
            self.logger.error(f"no dot files to process ")
 
    def _run_steps(self, timestamp):
        # Run the detection and stitching for this time period
        candidate_files, input_files = self._run_detection(timestamp)
        tracked_files = self._run_stitching(candidate_files)
        edited_files = self._run_node_file_editor(tracked_files, input_files)

        # write out netcdf files 
        for index, track_type in enumerate(self.track_types):
            tracked_file = tracked_files[index]
            if os.stat(tracked_file).st_size > 0:
                ncfile = tracked_file+'.nc'
                nrecords = 10
                if not os.path.exists(ncfile):
                    self.logger.error(f"open netcdf file {ncfile}")
                    nc_out = self._create_netcdf(self.output_directory, ncfile, nrecords,
                      STARTDATE=timestamp, ENDDATE=timestamp, MODEL=self.um_runid)

        # Run the plotting for this time period (if required and
        # data available)
        for index, track_type in enumerate(self.track_types):
            candidate_file = candidate_files[index]
            tracked_file = tracked_files[index]
            if os.stat(candidate_file).st_size > 0:
                if os.stat(tracked_file).st_size > 0:
                    if self.plot_tracks:
                        self.logger.debug(f"Plotting {os.path.basename(tracked_file)}")

                        self._read_and_plot_tracks(
                            tracked_file,
                            input_files["pslfile"],
                            title_prefix=f"{self.um_runid} {self.resolution_code} "
                            f"{self.time_range}",
                            title_suffix=f"{track_type} tracks",
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
            annual_tracks = self._run_annual_stitch()
            if self.plot_tracks:
                for index, track_type in enumerate(self.track_types):
                    self._read_and_plot_tracks(
                        annual_tracks[index],
                        input_files["pslfile"],
                        title_prefix=f"{self.um_runid} {self.resolution_code} "
                        f"{timestamp[:4]}",
                        title_suffix=f"{track_type} annual tracks",
                    )

            # plot annual tracks (also fills in gaps in tracks)
            # create netcdf output file for tracks
        else:
            self.logger.debug("Not running annual stitch")

        # Produce an wholerun summary if this is the end of the run
        if self._enddate[0:8] == timestamp[0:8]:
            wholerun_tracks = self._run_wholerun_stitch()
            if self.plot_tracks:
                for index, track_type in enumerate(self.track_types):
                    self._read_and_plot_tracks(
                        wholerun_tracks[index],
                        input_files["pslfile"],
                        title_prefix=f"{self.um_runid} {self.resolution_code} "
                        f"{timestamp[:4]}",
                        title_suffix=f"{track_type} annual tracks",
                    )

            # plot annual tracks (also fills in gaps in tracks)
            # create netcdf output file for tracks
        else:
            self.logger.debug("Not running annual stitch")

    def _run_detection(self, timestamp):
        """
        Run the Tempest detection.

        :returns: The path to the candidate files (as a list ordered by track
            type) and details of the processed input files (as a dict).
        :rtype: tuple
        """
        self.logger.debug(f"cwd {os.getcwd()}")

        filenames = self._generate_data_files(timestamp)

        candidate_files = []

        for track_type in self.track_types:
            self.logger.debug(f"Runing {track_type} detection")
            candidatefile = os.path.join(
                self.output_directory,
                f"candidate_file_{timestamp}_{track_type}.txt",
            )
            self.logger.debug(f"candidatefile {candidatefile}")

            if not self.extended_files:
                # identify candidates
                # TODO there is a file difference here with the run version
                cmd_io = '{} --in_data "{};{};{};{};{}" --out {} '.format(
                    self.tc_detect_script,
                    filenames["pslfile"],
                    filenames["zgfile"],
                    filenames["ufile"],
                    filenames["vfile"],
                    filenames["topofile"],
                    candidatefile,
                )
            else:
                cmd_io = '{} --in_data "{};{};{};{};{};{}" ' "--out {} ".format(
                    self.tc_detect_script,
                    filenames["pslfile"],
                    filenames["zgfile"],
                    filenames["u10mfile"],
                    filenames["v10mfile"],
                    filenames["rvT63file"],
                    filenames["topofile"],
                    candidatefile,
                )

            tracking_phase_commands = self._construct_command(track_type)
            cmd_detect = cmd_io + tracking_phase_commands["detect"]
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

        return candidate_files, filenames

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
            self.logger.debug(f"Runing {track_type} stitching")

            candidatefile = candidate_files[index]
            trackedfile = os.path.join(
                self.output_directory, f"track_file_{self.time_range}_{track_type}.txt"
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
                self.output_directory,
                f"candidate_file_{previous_year}*_{track_type}.txt",
            )
            candidate_files = sorted(glob.glob(candidate_pattern))
            annual_candidate = os.path.join(
                self.output_directory,
                f"candidate_year_{previous_year}_{track_type}.txt",
            )
            with open(annual_candidate, "w") as out_file:
                for candidate_file in candidate_files:
                    with open(candidate_file) as in_file:
                        for line in in_file:
                            out_file.write(line)

            annual_track = os.path.join(
                self.output_directory, f"track_year_{previous_year}_{track_type}.txt"
            )
            self._stitch_file(annual_candidate, annual_track, track_type)
            tracked_files.append(annual_track)

        return tracked_files

    def _run_wholerun_stitch(self, timestamp):
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
                self.output_directory,
                f"candidate_file_????????_{track_type}.txt",
            )
            candidate_files = sorted(glob.glob(candidate_pattern))
            wholerun_candidate = os.path.join(
                self.output_directory,
                f"candidate_year_{start_period}_{end_period}_{track_type}.txt",
            )
            with open(wholerun_candidate, "w") as out_file:
                for candidate_file in candidate_files:
                    with open(candidate_file) as in_file:
                        for line in in_file:
                            out_file.write(line)

            wholerun_track = os.path.join(
                self.output_directory, f"track_year_{start_period}_{end_period}_{track_type}.txt"
            )
            self._stitch_file(wholerun_candidate, annual_track, track_type)
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

    def _run_node_file_editor(self, tracked_files, input_files):
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

            if not self.extended_files:
                cmd_io = (
                    '{} --in_data "{};{};{};{};{}" --in_nodefile {} '
                    "--out_nodefile {} ".format(
                        self.tc_editor_script,
                        input_files["pslfile"],
                        input_files["zgfile"],
                        input_files["u10mfile"],
                        input_files["v10mfile"],
                        input_files["topofile"],
                        tracked_file,
                        edited_file,
                    )
                )
            else:
                cmd_io = (
                    '{} --in_data "{};{};{};{};{};{};{};{};{};{};{};'
                    '{};{}" in_nodefile {} --out_nodefile {} '.format(
                        self.tc_editor_script,
                        input_files["pslfile"],
                        input_files["zgfile"],
                        input_files["ufile"],
                        input_files["vfile"],
                        input_files["rvfile"],
                        input_files["u10mfile"],
                        input_files["v10mfile"],
                        input_files["ws10mfile"],
                        input_files["viwvefile"],
                        input_files["viwvnfile"],
                        input_files["tafile"],
                        input_files["rvT63file"],
                        input_files["topofile"],
                        tracked_file,
                        edited_file,
                    )
                )

            cmd_edit = cmd_io + tracking_phase_commands["profile"]
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

    def _write_dot_track_file(self, timestamp, dot_file = 'do_tracking'):
        """
        Write a file indicating that this timestep needs to be tracked

        """
        do_tracking_file = os.path.join(self.output_directory, dot_file+'.'+timestamp)
        if not os.path.exists(do_tracking_file):
            os.system('touch '+do_tracking_file)

    def _remove_dot_track_file(self, timestamp, dot_file = 'do_tracking'):
        """
        Write a file indicating that this timestep needs to be tracked

        """
        do_tracking_file = os.path.join(self.output_directory, dot_file+'.'+timestamp)
        if os.path.exists(do_tracking_file):
            os.system('rm '+do_tracking_file)

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

    def _generate_data_files(self, timestamp):
        """
        Identify and then fix the grids and var_names in the input files.

        The time_range and frequency attributes are set when this method runs.

        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """
        #timestamp_day = self.cylc_task_cycle_time[:8]
        timestamp_day = timestamp
        self.logger.debug(f"time_stamp_day {timestamp_day}")
        file_search = os.path.join(
            self.input_directory, f"atmos_{self.um_runid}*{timestamp_day}??-*-slp.nc"
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

        filenames = self._process_input_files()

        return filenames

    def _process_input_files(self):
        """
        Identify and then fix the grids and var_names in the input files.

        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """
        filetypes_required = [
            "pslfile",
            "zgfile",
            "ufile",
            "vfile",
            "ws10mfile",
            "rvfile",
            "rvT63file",
            "u10mfile",
            "u10mfile",
            "v10mfile",
            "viwvefile",
            "viwvnfile",
            "tafile",
        ]
        processed_filenames = {}

        variables_required = {}
        variables_required["pslfile"] = {"fname": "slp"}
        variables_required["zgfile"] = {"fname": "zg", "varname_new": "zg"}
        variables_required["ufile"] = {"fname": "ua"}
        variables_required["vfile"] = {"fname": "va"}
        variables_required["ws10mfile"] = {"fname": "wsas"}
        variables_required["rvfile"] = {"fname": "rv", "varname_new": "rv"}
        variables_required["rvT63file"] = {"fname": "rvT63", "varname_new": "rvT63"}
        variables_required["u10mfile"] = {"fname": "u10m", "varname_new": "uas"}
        variables_required["v10mfile"] = {"fname": "v10m", "varname_new": "vas"}
        variables_required["viwvefile"] = {"fname": "viwve", "varname_new": "viwve"}
        variables_required["viwvnfile"] = {"fname": "viwvn", "varname_new": "viwvn"}
        variables_required["tafile"] = {"fname": "ta", "varname_new": "ta"}

        filename_format = "atmos_{}a_{}h_{}_pt-{}.nc"

        reference_name = filename_format.format(
            self.um_runid,
            self.frequency,
            self.time_range,
            variables_required["pslfile"]["fname"],
        )
        reference_path = os.path.join(self.input_directory, reference_name)
        reference = iris.load_cube(reference_path)

        for filetype in filetypes_required:
            filename = filename_format.format(
                self.um_runid,
                self.frequency,
                self.time_range,
                variables_required[filetype]["fname"],
            )

            input_path = os.path.join(self.input_directory, filename)
            if not os.path.exists(input_path):
                msg = f"Unable to find expected input file {input_path}"
                self.logger.error(msg)
                raise RuntimeError(msg)

            output_path = os.path.join(self.output_directory, filename)

            if variables_required[filetype]["fname"] in [
                "ua",
                "va",
                "wsas",
                "rv",
                "rvT63",
                "u10m",
                "v10m",
            ]:
                # regrid u and v to t grid and rename variable if necessary
                cube = iris.load_cube(input_path)
                regridded = cube.regrid(reference, iris.analysis.Linear())
                if "varname_new" in variables_required[filetype]:
                    regridded.var_name = variables_required[filetype]["varname_new"]
                iris.save(regridded, output_path)
            elif variables_required[filetype]["fname"] in ["ta", "zg"]:
                # rename variables only - could do with ncrename instead
                cube = iris.load_cube(input_path)
                if "varname_new" in variables_required[filetype]:
                    cube.var_name = variables_required[filetype]["varname_new"]
                iris.save(cube, output_path)
            else:
                if not os.path.exists(output_path):
                    shutil.copyfile(input_path, output_path)

            processed_filenames[filetype] = output_path

        # Identify the grid and orography file
        longitude_size = reference.shape[-1]
        resolution = longitude_size // 2
        self.resolution_code = f"N{resolution}"

        processed_filenames["topofile"] = os.path.join(
            self.orography_dir, f"orog_HadGEM3-GC31-{self.resolution_code}e.nc"
        )
        self.logger.debug(f"Orography file {processed_filenames['topofile']}")

        return processed_filenames

    def _read_and_plot_tracks(
        self,
        tracked_file,
        nc_file,
        title_prefix="",
        title_suffix="TempestExtremes TCs",
        include_num_tracks=True,
    ):
        """
        Read and then plot the tracks. The title is the specified prefix, the
        number of tracks (if requested) and then the suffix.

        :param str tracked_file: The path to the track file.
        :param str nc_file: The path to an input data netCDF file, which is
            used to gather additional information about the dates and calendars
            used in the data.
        :param str title_suffix: The title for the plot.
        :param str title_suffix: The title for the plot.
        :param bool include_num_tracks: Include the number of tracks loaded?
        """

        storms = get_trajectories(tracked_file, nc_file, self.frequency)
        self.logger.debug(f"got storms {storms}")

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
        plot_trajectories_cartopy(storms, filename, title=title_full)

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
        self.psl_std_name = self.app_config.get_property("common", "psl_std_name")
        self.orography_dir = self.app_config.get_property("common", "orography_dir")
        self.extended_files = self.app_config.get_bool_property(
            "common", "extended_files"
        )
        self.plot_tracks = self.app_config.get_bool_property("common", "plot_tracks")
        # track_types is a Python list so eval converts str to list
        self.track_types = eval(self.app_config.get_property("common", "track_types"))

    def _get_environment_variables(self):
        """
        Get the required environment variables from the suite. A list and
        explanation of the required environment variables is included in the
        documentation.
        """
        self.um_runid = os.environ["RUNID"]
        self.cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
        self.time_cycle = os.environ["TIME_CYCLE"]
        self.previous_cycle = os.environ["PREVIOUS_CYCLE"]
        self.startdate = os.environ["STARTDATE"]
        self.enddate = os.environ["ENDDATE"]

    def _create_netcdf(self, directory, savefname, nrecords,
                      STARTDATE=None, ENDDATE=None,
                      MODEL=None):
        """
        Create netcdf file for the tracks. 
        May need metadata from a model nc file, so may need to create at a time when these are available
        """
        print ('making netCDF of outputs')
        
        self.savefname = savefname
        nc = Dataset(self.savefname, 'w', format='NETCDF4')
        nc.title = 'Tempest TC tracks'
        nc.directory = directory
        #nc.HOURS_BTWN_RECORDS = np.float64(self.HOURS_BTWN_RECORDS)
        #nc.TRACK_DURATION_MIN = np.float64(self.TRACK_DURATION_MIN)

        nc.MODEL = MODEL
        nc.YEAR_MIN = np.int32(STARTDATE[0:4])
        nc.YEAR_MAX = np.int32(ENDDATE[0:4])
        nc.MON_MIN = np.int32(STARTDATE[4:6])
        nc.MON_MAX = np.int32(ENDDATE[4:6])
        calendar = '360_day'

        nc.createDimension('tracks', size=None) # unlimited
        nc.createDimension('record', size = nrecords)

        nc.createVariable('first_pt', np.int32, ('tracks'), fill_value=self.FILLVAL)   
        nc.createVariable('num_pts', np.int32, ('tracks'), fill_value=self.FILLVAL)  
        nc.createVariable('track_id', np.int32, ('tracks'), fill_value=self.FILLVAL)  
        nc.createVariable('index', np.int32, ('record'), fill_value=self.FILLVAL)  
        nc.createVariable('time', 'f8', ('record'), fill_value=self.FILLVAL)
        nc.createVariable('lon', 'f4', ('record'), fill_value=self.FILLVAL)
        nc.createVariable('lat', 'f4', ('record'), fill_value=self.FILLVAL)
        nc.createVariable('slp', 'f4', ('record'), fill_value=self.FILLVAL)
        nc.createVariable('wind10m', 'f4', ('record'), fill_value=self.FILLVAL)
        nc.createVariable('surface_altitude', 'f4', ('record'), fill_value=self.FILLVAL)

        nc.variables['first_pt'].units = 'ordinal'
        #nc.variables['first_pt'].min_val = np.int32(0)
        #nc.variables['first_pt'].max_val = np.int32(0)
        nc.variables['first_pt'].long_name = 'first_pt'
        nc.variables['first_pt'].description = 'Index to first point of this track number'

        nc.variables['num_pts'].units = 'ordinal'
        nc.variables['num_pts'].long_name = 'num_pts'
        # Add option here to set length of intervals.....
        nc.variables['num_pts'].description = 'Number of points for this track'

        nc.variables['track_id'].units = 'ordinal'
        nc.variables['track_id'].long_name = 'track_id'
        nc.variables['track_id'].description = 'Tropical cyclone track number'

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

        time_units = 'days since 1950-01-01'
        nc.variables['time'].units = time_units
        nc.variables['time'].calendar = calendar
        nc.variables['time'].standard_name = 'time'
        nc.variables['time'].long_name = 'time'

        nc.variables['slp'].units = 'hPa'
        nc.variables['slp'].missing_value = 1.0e+25
        nc.variables['slp'].standard_name = 'air_pressure_at_mean_sea_level'
        nc.variables['slp'].long_name = 'Sea Level Pressure'
        nc.variables['slp'].description = 'Sea level pressure for tracked variable'

        nc.variables['wind10m'].units = 'm s-1'
        nc.variables['wind10m'].missing_value = 1.0e+25
        nc.variables['wind10m'].standard_name = 'awind_speed'
        nc.variables['wind10m'].long_name = 'Near-surface Wind Speed'
        nc.variables['wind10m'].description = 'near-surface (usually 10 metres) wind speed'

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

 
