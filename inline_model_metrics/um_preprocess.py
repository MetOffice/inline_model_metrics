# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import iris
# windspharm needs installing: see https://ajdawson.github.io/windspharm/latest/
#from windspharm.iris import VectorWind

from .tempest_common import (TempestExtremesAbstract)

#from afterburner.apps import AbstractApp

class TempestError(Exception):
    """
    Custom TempestExtremes tracking exception
    """
    pass


class UMTempestPreprocess(TempestExtremesAbstract):
    """
    Preprocess Unified Model (Met Office HadGEM model) data for the
    TempestExtremes trackers.
    """

    def __init__(self, arglist=None, **kwargs):
        super().__init__(version="1.0", **kwargs)
        self._parse_args(arglist, desc="Inline Model Metrics Pre-Processor")
        self._parse_app_config()
        self._set_message_level()

        # Instance attributes set later
        self.time_range = None
        self.frequency = None
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

        output_dir_native = self.output_directory+'_native'
        if not os.path.exists(output_dir_native):
            try:
                os.makedirs(output_dir_native)
            except PermissionError:
                msg = f"Unable to create output directory {output_dir_native}"
                self.logger.error(msg)
                sys.exit(1)
        output_dir_native_archive = os.path.join(self.output_directory + '_native', self._archived_files_dir)
        if not os.path.exists(output_dir_native_archive):
            try:
                os.makedirs(output_dir_native_archive)
            except PermissionError:
                msg = f"Unable to create output directory {output_dir_native_archive}"
                self.logger.error(msg)
                sys.exit(1)

        # initialise variables that might not get set if no detection step
        self.outdir = self.output_directory + "_" + "native"

        # find out what the psl variable is for the input data
        if self.variables_rename[0] == "psl":
            psl_input_var = self.variables_input[0]
        else:
            psl_input_var = "psl"

        self.logger.debug(
            f"CYLC_TASK_CYCLE_TIME {self.cylc_task_cycle_time}, "
            f"runid {self.runid}, psl_input_var {psl_input_var}"
        )

        timestamp_day = self.cylc_task_cycle_time[:8]
        timestamp_endday = self.next_cycle[:8]

        condition = self._set_tracking_date(timestamp_day)
        if condition == "AlreadyComplete":
            return

        # this section of code processes data from the current timestep
        current_time = timestamp_day
        # find the relevant input data using the given file pattern
        fname = self._file_pattern(current_time + "*", "*", psl_input_var,
                                   frequency="*", stream=self.um_stream)
        file_search = os.path.join(self.input_directory, fname)
        self.logger.debug(f"file_search {file_search}")

        # pre-process the input files and produce standard processed files for
        # later use
        if glob.glob(file_search):
            for regrid_resol in self.regrid_resolutions:
                self.outdir = self.output_directory + '_' + regrid_resol
                source_files, processed_files, variable_units = \
                    self._generate_data_files(timestamp_day, timestamp_endday,
                                              psl_input_var, grid_resol=regrid_resol)
                #self._produce_derived_diagnostics(source_files, processed_files)

    # def _file_pattern(self, timestart, timeend, varname,
    #                   frequency="6h", stream="pt"):
    #     """
    #     Derive the input nc filenames from the file pattern, assuming a
    #     um model filenaming pattern as here, or could be other patterns
    #     for other models/platforms (which would need to be added)
    #
    #     :param str timestart: The timestep of the start of the data period to process
    #     :param str timeend: The timestep of the end of the data period to process
    #     :param str um_stream: The name of the um output stream (output file
    #     :                     identification)
    #     :param str frequency: The frequency of the input data (in hours, needs to
    #     :                     include "h"), used to determine file naming
    #     :returns: a filename given the inputs to the pattern
    #     :rtype: str
    #     """
    #     if self.frequency is None:
    #         file_freq = frequency
    #     else:
    #         file_freq = str(self.frequency)+"h"
    #
    #     if self.input_file_pattern != '':
    #         # file format based on input pattern
    #         fname = self.input_file_pattern.format(
    #             runid=self.runid,
    #             frequency=file_freq,
    #             date_start=timestart,
    #             date_end=timeend,
    #             variable=varname,
    #             stream=stream
    #         )
    #     self.logger.info(f"fname from _file_pattern {fname} {self.um_stream} {timestart}" +\
    #                     f" {timeend} {varname}")
    #     return fname.strip('"')
    #
    # def _file_pattern_processed(self, timestart, timeend, varname,
    #                   frequency="6h"):
    #     """
    #     For processed files, we know what the filenames look like, so
    #     search specifically
    #
    #     :param str timestart: The timestep of the start of the data period to process
    #     :param str timeend: The timestep of the end of the data period to process
    #     :param str frequency: The frequency of the input data (in hours, needs to
    #     :                     include "h"), used to determine file naming
    #     :returns: a filename given the inputs to the pattern
    #     :rtype: str
    #     """
    #     if self.frequency is None:
    #         file_freq = frequency
    #     else:
    #         file_freq = str(self.frequency)+"h"
    #
    #     fname = self.file_pattern_processed.format(
    #         runid=self.runid,
    #         frequency=file_freq,
    #         date_start=timestart,
    #         date_end=timeend,
    #         variable=varname
    #     )
    #
    #     self.logger.info(f"fname from _file_pattern_processed {fname} {timestart} " + \
    #                      "{timeend} {varname}")
    #     return fname.strip('"')

    def _generate_data_files(self, timestamp, timestamp_end,
                             psl_var, grid_resol="native"):
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
        fname = self._file_pattern(timestamp_day+"*", "*", psl_var,
                                   frequency="*", stream=self.um_stream)
        file_search = os.path.join(
            self.input_directory, fname
        )
        files = sorted(glob.glob(file_search))
        self.logger.debug(f"file_search in _generate_data_files {file_search}, " + \
                          "files {files}")
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
                time_range = time_range.replace(".nc", "")
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
            self._process_input_files(timestamp, timestamp_end,
                                      psl_var, grid_resol=grid_resol)

        return source_files, processed_files, variable_units

    def _check_time_coord(self, fnames):
        """
        Check that file has time coordinates called
        time, if not change the name

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
                        " -d " + time_name + ",time -v " + time_name + ",time " +\
                      fname+".nc3"
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)

                cmd = "mv " + fname + ".nc3" + " " + fname
                self.logger.debug(f"cmd {cmd}")
                subprocess.call(cmd, shell=True)

    def _check_um_resolution(self, cube, fname):
        """
        For UM model specifically, check that input resolution code agrees with
        input file size
        For UM, the n???e number is half of the number of zonal grid points
        """
        longitude_size = cube.shape[-1]
        cube_resolution = longitude_size // 2
        input_resolution = int(self.resolution_code[1:-1])
        print("input, real resol ",self.resolution_code, cube_resolution)
        if cube_resolution != input_resolution:
            msg = f"Resolution of input file {fname} not consistent with resolution " + \
                "code {self.resolution_code}"
            self.logger.error(msg)
            raise RuntimeError(msg)

    def _process_input_files(self, time_start, time_end, psl_var,
                             grid_resol="native"):
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
        for iv, var in enumerate(self.variables_input):
            variables_required[var] = {"fname": var}
            #if var in self.variables_rename:
            var_rename = self.variables_rename[iv]
            variables_required[var].update({"varname_new": var_rename})

        reference_name = self._file_pattern(self.time_range.split("-")[0],
                                            self.time_range.split("-")[1],
                                            variables_required[psl_var]["fname"],
                                            stream=self.um_stream)
        reference_path = os.path.join(self.input_directory, reference_name)
        reference = iris.load_cube(reference_path)
        variable_units[psl_var] = reference.units

        # Identify the grid and orography file
        if grid_resol == "native":
            # check resolution code in Met Office UM mdoel
            if self.suiteid[0:2] == "u-":
                self._check_um_resolution(reference, reference_path)

            processed_filenames["orog"] = os.path.join(
                self.orography_dir, f"orography-{self.resolution_code}.nc"
            )
            shutil.copyfile(processed_filenames["orog"],
                            os.path.join(self.outdir, "orography.nc"))
        else:
            processed_filenames["orog"] = os.path.join(
                self.orography_dir, f"orography-{grid_resol}.nc"
            )
        cube_orog = iris.load_cube(processed_filenames["orog"])
        variable_units["orog"] = cube_orog.units
        reference = cube_orog
        if cube_orog.var_name != "surface_altitude":
            cube_orog.var_name = "surface_altitude"
        cube_orog = iris.util.squeeze(cube_orog)
        iris.save(cube_orog, os.path.join(self.outdir, "orography.nc"))

        for var in self.variables_input:
            filename = self._file_pattern(self.time_range.split("-")[0],
                                          self.time_range.split("-")[1],
                                          variables_required[var]["fname"],
                                          stream=self.um_stream)

            #input_path = os.path.join(self.input_directory, filename)
            #if not os.path.exists(input_path):
            fin_path = os.path.join(self.input_directory, filename)
            file_search = glob.glob(fin_path)
            if len(file_search) != 1:
                msg = f"Unable to find expected input file {input_path}"
                self.logger.error(msg)
                raise RuntimeError(msg)
            else:
                input_path = file_search[0]

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

            # apart from psl, regrid the data to the psl grid (all input data for
            # TempestExtremes needs to be on the same grid)
            # regrid u and v to t grid and rename variable if necessary
            cube = iris.load_cube(input_path)
            cube = iris.util.squeeze(cube)
            ndims = cube.shape
            variable_units[var] = cube.units
            if var == "uas":
                variable_units["sfcWind"] = cube.units
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
                    regridded_level.var_name = regridded.var_name +\
                                               "_" + str(int(level))
                    iris.save(regridded_level, output_path[:-3] + "_" +
                              str(int(level)) + ".nc")
                    output_paths.append(output_path[:-3]+"_" +
                                        str(int(level)) + ".nc")
            check_paths = [output_paths] if isinstance(output_paths, str)\
                                        else output_paths
            self._check_time_coord(check_paths)

            source_filenames[var] = input_path
            processed_filenames[var] = output_paths

        self.logger.debug(f"Orography file {processed_filenames['orog']}")

        return source_filenames, processed_filenames, variable_units

    def _produce_derived_diagnostics(self, dir_proc, processed_filenames):
        """
        Produce any derived diagnostics needed e.g. streamfunction

        :param str dir_proc: Directory containing processed variables
        :param dict processed_filenames: Dictionary of processed variables
            and filenames
        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        """
        
        # first calculate the 850 hPa streamfunction
        if "ua_850" in processed_filenames and "va_850" in processed_filenames:
            uwnd = iris.load_cube(processed_filenames["ua_850"])
            vwnd = iris.load_cube(processed_filenames["va_850"])
            if 'mask' in dir(uwnd.data): uwnd.data = np.ma.filled(uwnd.data, 0)
            if 'mask' in dir(vwnd.data): vwnd.data = np.ma.filled(vwnd.data, 0)

            # Create a VectorWind instance to handle the computation of streamfunction
            w = VectorWind(uwnd, vwnd)

            # Compute the streamfunction
            sf = w.streamfunction()
            sf.var_name = "sf_850"
            fname_ua_850 = processed_filenames["ua_850"]
            fname_sf_850 = fname_ua_850.replace("ua_850", "sf_850")
            iris.save(sf, fname_sf_850)


    def _get_app_options(self):
        """Get commonly used configuration items from the config file"""

        self.input_directory = self.app_config.get_property("common",
                                                            "input_directory")
        self.output_directory = self.app_config.get_property(
            "common", "output_directory"
        )
        self.orography_dir = self.app_config.get_property("common", "orography_dir")
        self.variables_input = eval(self.app_config.get_property("common",
                                                                 "variables_input"))
        self.variables_rename = eval(self.app_config.get_property("common",
                                                                 "variables_rename"))
        self.input_file_pattern = self.app_config.get_property("common",
                                                            "input_file_pattern")
        self.file_pattern_processed = self.app_config.get_property("common",
                                                            "file_pattern_processed")
        self.regrid_resolutions = eval(self.app_config.get_property(
            "common",
            "regrid_resolutions"
        ))
        self.data_frequency = self.app_config.get_property("common",
                                                            "data_frequency")
        try:
            self.um_stream = eval(self.app_config.get_property("common",
                                                            "um_stream"))
        except:
            self.um_stream = "pt"

    def _get_environment_variables(self):
        """
        Get the required environment variables from the suite. A list and
        explanation of the required environment variables is included in the
        documentation.
        """
        try:
            self.suiteid = os.environ["SUITEID_OVERRIDE"]
        except:
            self.suiteid = os.environ["SUITE"]
        try:
            self.runid = self.suiteid.split('-')[1]
        except:
            self.runid = self.suiteid
        self.resolution_code = os.environ["RESOL_ATM"]
        print('resol in ',self.resolution_code)
        self.cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
        self.next_cycle = os.environ["NEXT_CYCLE"]
