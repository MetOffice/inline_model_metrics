# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
from abc import ABCMeta, abstractmethod
import glob
import os
import re
import shutil
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
        self.resolution_code = None
        self.cmd_detect_type = {}
        self.cmd_stitch_type = {}
        self.cmd_edit_type = {}
        self.cmd_detect = None
        self.cmd_detectblobs = None
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
        self.variables_rename = []
        self.outputcmd_detect_default = ""
        self.in_fmt_stitch_default = ""
        self.out_fmt_profile1_default = ""
        self.out_fmt_profile2_default = ""

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
        #TODO an identical _generate_data_files() is called in preprocess. This one should be removed
        source_files, processed_files, variable_units = self._generate_data_files(timestamp)

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
            values: detect, stitch, profile
        :returns: A dictionary with keys of `detect`, `stitch`, `profile`
            and the values for each of these is a string containing the
            command line parameters for each of these TempestExtreme steps.
            The parameters are sorted into alphabetical order in each line.
        :rtype: dict
        """

        # These are the first and last column headers in the standard
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
        for var in self.variables_input:
            variables_required[var] = {'fname': var}
            if var in self.variables_rename:
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
            # TODO self.resolution_code is set in two places in the code
            self.resolution_code = f"N{resolution}"

            # HadGEM specific - make generic
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

        #TODO self.resolution_code is set in two places in the code
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

        # TODO check that the variables here match those in the documentation
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
