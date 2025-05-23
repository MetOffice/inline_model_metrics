# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
from abc import ABCMeta, abstractmethod
import os
import glob
import subprocess

import iris

from afterburner.apps import AbstractApp


class TempestError(Exception):
    """
    Custom TempestExtremes tracking exception
    """
    pass


class TempestExtremesAbstract(AbstractApp, metaclass=ABCMeta):
    """
    An abstract base TempestExtremes class tracker to run inline with a climate
    model.
    """

    def __init__(self, arglist=None, **kwargs):
        super().__init__(**kwargs)

        # Instance attributes set later from Rose config
        self.time_range = None
        self.frequency = None
        self.cmd_detect_type = {}
        self.cmd_stitch_type = {}
        self.cmd_detectblobs_type = {}
        self.cmd_stitchblobs_type = {}
        self.cmd_blobstats_type = {}
        self.cmd_edit_type = {}
        self.cmd_varproc_type = {}
        self.cmd_detect = None
        self.cmd_detectblobs = None
        self.cmd_stitchblobs = None
        self.cmd_blobstats = None
        self.cmd_stitch = None
        self.cmd_edit = None
        self.source_files = {}
        self.processed_files = {}
        self.processed_files_psl = {}
        self.variable_units = {}
        self.calendar_units = "days since 1950-01-01 00:00:00"
        self.regrid_resolutions = None
        self.outdir = None
        self.column_names = {}
        self._is_new_year = False
        self._old_year_value = None
        self._archived_files_dir = "archived_files"
        self.variables_rename = []
        self.outputcmd_detect_default = ""
        self.in_fmt_stitch_default = ""
        self.outputcmd_detect_default1 = ""
        self.in_fmt_stitch_default1 = ""
        self.out_fmt_nodeedit1_default = ""
        self.out_fmt_nodeedit2_default = ""

    @abstractmethod
    def run(self, *args, **kwargs):
        """
        Run the app
        """
        pass

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

    def _set_tracking_date(self, timestamp,
                              dot_file='current_tracking_date'):
        """
        Write a file indicating that this timestep is where the tracking has reached

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param str dot_file: The first part of the string of a filename to indicate
        :                    which time periods still need tracking

        """
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        current_tracking_file = os.path.join(self.outdir, dot_file+'.' + timestamp)
        previous_tracking_file = glob.glob(os.path.join(self.outdir, dot_file+'.*'))

        if len(previous_tracking_file) == 1:
            date_previous = os.path.basename(previous_tracking_file[0]).split('.')[1]
            if _is_date_after_or_equal(timestamp, date_previous):
                os.remove(previous_tracking_file[0])
                os.system('touch ' + current_tracking_file)
                return "DoTracking"
            else:
                return "AlreadyComplete"
        elif len(previous_tracking_file) == 0:
            os.system('touch ' + current_tracking_file)
            return "DoTracking"
        else:
            raise Exception('Too many previous tracking file dates')

    def _get_tracking_date(self, timestamp,
                              dot_file='current_tracking_date'):
        """
        Write a file indicating that this timestep is where the tracking has reached

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param str dot_file: The first part of the string of a filename to indicate
        :                    which time periods still need tracking

        """
        previous_tracking_file = glob.glob(os.path.join(self.outdir, dot_file+'.*'))
        if len(previous_tracking_file) == 1:
            date_previous = os.path.basename(previous_tracking_file[0]).split('.')[1]
            if _is_date_after_or_equal(timestamp, date_previous):
                return "DoTracking"
            else:
                return "AlreadyComplete"
        elif len(previous_tracking_file) == 0:
            return "DoTracking"
        else:
            raise Exception('Too many previous tracking file dates')

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

    def _tidy_data_files(self, timestamp, timestamp_end, var_list,
                         f_remove=False):
        """
        Remove processed input files for this timestamp (tidy up)

        :param str timestamp: The timestep of the start of the data period to process
        :param str timestamp_end: The timestep of the end of the data period to process
        :param list var_list: List of variable names for file pattern deletion
        :param str f_remove: An indicator of which files need to be deleted.
                   processed = the (regridded) files read by the tracking code
                   source    = the files produced by the model
        """
        self.logger.info(f"Tidy up input files")
        files_remove = []

        if f_remove:
            for var in var_list:
                f = self._file_pattern_processed(timestamp, timestamp_end, var,
                                                 frequency=self.data_frequency)
                if os.path.exists(os.path.join(self.outdir, f)):
                    self.logger.info(f"Deleting {os.path.join(self.outdir, f)}")
                    os.remove(os.path.join(self.outdir, f))

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

    def _file_pattern(self, timestart, timeend, varname, stream="pt",
                      frequency="6h"):
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

        if self.input_file_pattern != '':
            # file format from postproc
            fname = self.input_file_pattern.format(
                runid=self.runid,
                frequency=file_freq,
                date_start=timestart,
                date_end=timeend,
                stream=stream,
                variable=varname
            )
            #else:
            #    # file format from direct STASH to netcdf conversion
            #    fname = self.input_file_pattern.format(
            #        runid=self.runid,
            #        stream=stream,
            #        date_start=timestart,
            #        variable=varname
            #    )
        #self.logger.info(f"fname from pattern {fname} {stream} {timestart} "
        #                f"{timeend} {varname}")
        return fname.strip('"')

    def _file_pattern_processed(self, timestart, timeend, varname,
                                frequency='6h'):
        """
        For processed files, we know what the filenames look like, so
        search specifically

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
            runid=self.runid,
            frequency=file_freq,
            date_start=timestart,
            date_end=timeend,
            variable=varname
        )

        #self.logger.info(f"fname from _file_pattern_processed {fname} {timestart} " + \
        #                 "{timeend} {varname}")
        return fname.strip('"')

    def _construct_command(self, track_type):
        """
        Read the TempestExtreme command line parameters from the configuration.
        :param str track_type: The name of the type of tracking to run, possible
            values: detect, stitch, nodeedit
        :returns: A dictionary with keys of `detect`, `stitch`, `nodeedit`,
            `detectblobs`, `nodefilefilter` `varproc1`, `varproc2`, and the values for each of these is
            a string containing the command line parameters for each of these
            TempestExtreme steps. The parameters are sorted into alphabetical
            order in each line.
        :rtype: dict
        """

        # These are the first and last column headers in the standard
        # Tempest output files
        column_initial = "grid_x,grid_y,"
        column_final = ",year,month,day,hour"

        commands = {}
        fmt_value = {}
        for step in ["detect", "stitch", "nodeedit", "detectblobs", "nodefilefilter", "varproc1", "varproc2", "stitchblobs", "blobstats"]:
            try:
                step_config = self.app_config.section_to_dict(f"{track_type}_{step}")
                step_arguments = []
                for parameter in sorted(list(step_config.keys())):
                    #print('step_config ',step, step_config, parameter)
                    if step_config[parameter]:
                        if "_default1" in step_config[parameter]:
                            if "outputcmd" in step_config[parameter]:
                                param_value = self.outputcmd_detect_default1
                            elif "in_fmt_stitch" in step_config[parameter]:
                                param_value = self.in_fmt_stitch_default1
                                fmt_value['stitch'] = self.in_fmt_stitch_default1
                            elif "out_fmt_nodeedit1" in step_config[parameter]:
                                param_value = self.out_fmt_nodeedit1_default
                                fmt_value['nodeedit'] = self.out_fmt_nodeedit1_default
                            elif "out_fmt_nodeedit2" in step_config[parameter]:
                                param_value = self.out_fmt_nodeedit2_default
                                fmt_value['nodeedit'] = self.out_fmt_nodeedit2_default
                        elif "_default2" in step_config[parameter]:
                            if "outputcmd" in step_config[parameter]:
                                param_value = self.outputcmd_detect_default2
                            elif "in_fmt_stitch" in step_config[parameter]:
                                param_value = self.in_fmt_stitch_default2
                                fmt_value['stitch'] = self.in_fmt_stitch_default2
                        elif "_default3" in step_config[parameter]:
                            if "out_fmt_nodeedit" in step_config[parameter]:
                                param_value = self.out_fmt_nodeedit_default3
                                fmt_value['nodeedit'] = self.out_fmt_nodeedit_default3
                            elif "in_fmt_nodeedit" in step_config[parameter]:
                                param_value = self.in_fmt_nodeedit_default3
                        else:
                            param_value = step_config[parameter]
                            if "in_fmt" in parameter:
                                fmt_value["stitch"] = step_config[parameter]
                            elif "out_fmt" in parameter:
                                fmt_value["nodeedit"] = step_config[parameter]

                        # arguments that take no values
                        if "regional" in parameter or "tagonly" in parameter or "findblobs" in parameter or "out_headers" in parameter or "out_fulltime" in parameter:
                            if step_config[parameter]:
                                step_arguments.append(f"--{parameter} ")
                        else:
                            step_arguments.append(f"--{parameter} {param_value}")

                commands[step] = " ".join(step_arguments)
                self.logger.debug(f"step, commands {step} {commands[step]}")
            except:
                commands[step] = None

            # set up the column names of the track output file, to be used for
            # naming the storm keys
            if step == "stitch" and commands[step] is not None:
                col_names = column_initial + fmt_value["stitch"].strip('\"') +\
                            column_final
                self.column_names[track_type+"_stitch"] = {}
                names = col_names.split(',')
                for im, name in enumerate(names):
                    self.column_names[track_type+"_"+step][name] = im
            if step == "nodeedit" and commands[step] is not None:
                col_names = column_initial + fmt_value["nodeedit"].strip('\"') +\
                            column_final
                self.column_names[track_type+"_nodeedit"] = {}
                names = col_names.split(',')
                for im, name in enumerate(names):
                    self.column_names[track_type+"_"+step][name] = im

        return commands

    def _command_parameter(self, track_type, command, argument):
        """
        Read the TempestExtreme command line parameters from the configuration,
        and return the requested value of the argument for the given command.
        :param str track_type: The name of the type of tracking to run, possible
            values: detect, stitch, nodeedit
        :param str command: The name of the TempestExtremes command (as 
            defined in the input namelist
        :param str track_type: The name of the argument
        :returns: The parameter value of the given argument.
        :rtype: str
        """

        step = command
        step_config = self.app_config.section_to_dict(f"{track_type}_{step}")
        if argument in sorted(list(step_config.keys())):
            value = step_config[argument]
        else:
            value = None
        return value

    def _check_time_coord(self, fnames):
        """
        Check that file has latitude and longitude coordinates called
        latitude, longitude, if not change the names

        :param list fnames: filenames in which to check coord names
        """
        for fname in fnames:
            cube = iris.load_cube(fname)
            time_name = cube.coord('time').var_name
            if time_name != 'time':
                cmd = os.path.join(self.ncodir, "ncks") \
                      + " -6 " + fname + " " + fname + ".nc3"
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)
                cmd = os.path.join(self.ncodir, "ncrename") + \
                        " -d " + time_name + ",time -v " + time_name +\
                        ",time " + fname+".nc3"
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)

                cmd = "mv " + fname + ".nc3" + " " + fname
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)

    def _identify_processed_files(self, time_start, time_end, grid_resol="native", variables=None):
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

        if variables is None:
            variables = self.variables_rename

        # Identify the grid and orography file
        processed_filenames["orog"] = os.path.join(
            self.outdir, "orography.nc"
        )
        cube = iris.load_cube(processed_filenames["orog"])
        variable_units["orog"] = cube.units

        if not os.path.exists(os.path.dirname(self.outdir)):
            raise Exception("Processed file directory does not exist, should come "\
                    "from pre-processing " + self.outdir)

        for var_name in variables:
            # identify the processed path filename similar to CMIP6 naming,
            # will be standard regardless of the input filename structure
            # varname_new, freq, time,
            output_path = self._file_pattern_processed(time_start,
                                                       time_end,
                                                       var_name,
                                                       self.data_frequency)
            output_path = os.path.join(self.outdir, output_path)

            #self.logger.debug(f"read file {output_path}")
            if os.path.exists(output_path):
                processed_filenames[var_name] = output_path
                cube = iris.load_cube(output_path)
                cube = iris.util.squeeze(cube)
                variable_units[var_name] = cube.units
                if var_name == "uas":
                    variable_units["sfcWind"] = cube.units
                    variable_units["wind"] = cube.units
            else:
                raise Exception("Processed file does not exist "+output_path)

            self._check_time_coord([output_path])

        self.logger.debug(f"Orography file {processed_filenames['orog']}")

        return processed_filenames, variable_units

    def _get_app_options(self):
        """Get commonly used configuration items from the config file"""

        self.input_directory = self.app_config.get_property("common",
                                                            "input_directory")
        self.output_directory = self.app_config.get_property(
            "common", "output_directory"
        )
        self.orography_dir = self.app_config.get_property("common", "orography_dir")
        self.delete_processed = self.app_config.get_bool_property(
            "common", "delete_processed"
        )
        self.delete_source = self.app_config.get_bool_property(
            "common", "delete_source"
        )
        self.plot_tracks = self.app_config.get_bool_property("common", "plot_tracks")
        #self.nodeedit_vars = eval(self.app_config.get_property("common",
        #                                                       "nodeedit_vars"))
        self.input_file_pattern = self.app_config.get_property("common",
                                                            "input_file_pattern")
        self.file_pattern_processed = self.app_config.get_property("common",
                                                            "file_pattern_processed")
        self.regrid_resolutions = \
            eval(self.app_config.get_property("common", "regrid_resolutions"))
        self.data_frequency = self.app_config.get_property("common",
                                                            "data_frequency")
        self.outputcmd_detect_default1 = self.app_config.get_property("common",
                                                            "outputcmd_detect_default1")
        self.outputcmd_detect_default2 = self.app_config.get_property("common",
                                                            "outputcmd_detect_default2")
        self.in_fmt_stitch_default1 = self.app_config.get_property("common",
                                                            "in_fmt_stitch_default1")
        self.in_fmt_stitch_default2 = self.app_config.get_property("common",
                                                            "in_fmt_stitch_default2")
        self.in_fmt_nodeedit_default3 = self.app_config.get_property("common",
                                                            "in_fmt_nodeedit_default3")
        self.out_fmt_nodeedit_default3 = self.app_config.get_property("common",
                                                            "out_fmt_nodeedit_default3")

    def _get_environment_variables(self):
        """
        Get the required environment variables from the suite. A list and
        explanation of the required environment variables is included in the
        documentation.
        """

        # TODO check that the variables here match those in the documentation
        try:
            self.suiteid = os.environ["SUITEID_OVERRIDE"]
        except:
            self.suiteid = os.environ["CYLC_SUITE_NAME"]
        try:
            self.runid = self.suiteid.split('-')[1]
        except:
            self.runid = self.suiteid

        self.cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
        self.enddate = os.environ["ENDDATE"]
        self.inline_tracking = os.environ["INLINE_TRACKING"]
        self.lastcycle = os.environ["LASTCYCLE"]
        self.ncodir = os.environ["NCODIR"]
        self.next_cycle = os.environ["NEXT_CYCLE"]
        self.previous_cycle = os.environ["PREVIOUS_CYCLE"]
        self.resolution_code = os.environ["RESOL_ATM"]
        self.startdate = os.environ["STARTDATE"]
        self.time_cycle = os.environ["TIME_CYCLE"]
        self.tm2_cycle = os.environ["TM2_CYCLE"]
        self.tp2_cycle = os.environ["TP2_CYCLE"]
        track_by_year = os.environ["TRACK_BY_YEAR"]
        track_at_end = os.environ["TRACK_AT_END"]
        is_last_cycle = os.environ["IS_LAST_CYCLE"]

        if is_last_cycle == "true":
            self.is_last_cycle = True
        else:
            self.is_last_cycle = False

        if "T" in track_by_year or "t" in track_by_year:
            self.track_by_year = True
        else:
            self.track_by_year = False

        if "T" in track_at_end or "t" in track_at_end:
            self.track_at_end = True
        else:
            self.track_at_end = False

def _is_date_after(timetest, timeref):
    """
    Check if timetest is after timeref.

    :returns: True if timetest is after timeref.
    :rtype: bool
    """

    if int(timetest[:8]) > int(timeref[:8]):
        return True
    else:
        return False

def _is_date_after_or_equal(timetest, timeref):
    """
    Check if timetest is equal or after timestamp.

    :returns: True if timetest is equal or after timeref.
    :rtype: bool
    """

    if int(timetest[:8]) >= int(timeref[:8]):
        return True
    else:
        return False
