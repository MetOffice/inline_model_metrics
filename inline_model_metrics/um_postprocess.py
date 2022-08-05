# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import glob
import os
import re
import shutil
import subprocess
import sys

from afterburner.apps import AbstractApp

from .tempest_common import (TempestExtremesAbstract)

class TempestError(Exception):
    """
    Custom TempestExtremes tracking exception
    """
    pass


class UMTempestPostprocess(AbstractApp):
    """
    Postprocess Unified Model (Met Office HadGEM model) data for the
    TempestExtremes trackers.
    """

    def __init__(self, arglist=None, **kwargs):
        super().__init__(version="1.0", **kwargs)
        self._parse_args(arglist, desc="Inline Model Metrics Post-Processor")
        self._parse_app_config()
        self._set_message_level()

        # Instance attributes set later
        self.time_range = None
        self.frequency = None
        self.resolution_code = None
        self.regrid_resolutions = None
        self.outdir = None
        self._archived_files_dir = "archived_files"

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
        self.outdir = self.output_directory + "_" + "native"

        self.logger.debug(
            f"CYLC_TASK_CYCLE_TIME {self.cylc_task_cycle_time}, "
            f"um_runid {self.um_runid}"
        )

        timestamp_day = self.cylc_task_cycle_time[:8]
        timestamp_endday = self.next_cycle[:8]
        timestamp_tm2 = self.tm2_cycle[:8]

        # this section of code processes data from the current timestep
        current_time = timestamp_day

        for regrid_resol in self.regrid_resolutions:
            self.outdir = self.output_directory + '_' + regrid_resol
            # process the AR files before archiving
            self._process_ar_for_archive(os.path.join(self.outdir,
                                                      self._archived_files_dir))
            if self.um_archive_to_mass:
                # run the archiving on any .arch files that exist
                self._archive_tracking_to_mass(os.path.join(self.outdir,
                                                        self._archived_files_dir))

        # delete source data if required
        if self.delete_source:
            for var in self.variables_input:
                fname = self._file_pattern(timestamp_tm2 + "*", "*", var,
                                       um_stream="pt", frequency="*")
                file_name = os.path.join(self.input_directory, fname)
                self.logger.debug(f"deleting source file_name {file_name}")
                if os.path.exists(file_name):
                    os.remove(file_name)

    def _archive_tracking_to_mass(
        self,
        archive_directory
    ):
        """
        Archive any files with .arch to Met Office MASS system
        :param str archive_directory: The directory to look for any .arch files
        """

        moosedir = "moose:/crum/{}/{}/"
        archive_files = glob.glob(os.path.join(archive_directory, "*.arch"))
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
                    os.remove(fname)
                    os.remove(fname_arch)
                self.logger.debug(sts.stdout)

    def _process_ar_for_archive(
        self,
        archive_directory,
        nc_compression="1",
        netcdf_format="4"
    ):
        """
        Process any AR files before archiving
        :param str archive_directory: The directory to look for any .arch files
        :param str nc_compression: Compression level for netcdf files
        :param netcdf_format: Format of netcdf files (netcdf4 by default)
        """

        archive_files = glob.glob(os.path.join(archive_directory, "*ARmask*.arch"))
        if len(archive_files) > 0:
            for fname_arch in archive_files:
                fname = fname_arch[:-5]
                if fname[-2:] == "nc":
                    # convert to netcdf4 and compress
                    cmd = "nccopy -" + netcdf_format + " -d " + nc_compression +\
                          " " + fname + " " + fname + ".nc4"
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
                        os.remove(fname)
                        os.rename(fname+".nc4", fname)

                if fname[-2:] == "nc":
                    # convert time coordinate to unlimited
                    cmd = os.path.join(self.ncodir, "ncks") + " --mk_rec_dmn time " + \
                          fname + " " + fname[:-3] + "_unlimited.nc"
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
                        os.remove(fname)
                        os.rename(fname[:-3] + "_unlimited.nc", fname)

    def _file_pattern(self, timestart, timeend, varname, um_stream="pt",
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
        if self.frequency is None:
            file_freq = frequency
        else:
            file_freq = str(self.frequency)+"h"

        if self.input_file_pattern != '':
            if "atmos" in self.input_file_pattern:
                # file format from postproc
                fname = self.input_file_pattern.format(
                    runid=self.um_runid,
                    frequency=file_freq,
                    date_start=timestart,
                    date_end=timeend,
                    stream=um_stream,
                    variable=varname
                )
            else:
                # file format from direct STASH to netcdf conversion
                fname = self.input_file_pattern.format(
                    runid=self.um_runid,
                    stream=um_stream,
                    date_start=timestart,
                    variable=varname
                )
        self.logger.info(f"fname from pattern {fname} {um_stream} {timestart} "
                        f"{timeend} {varname}")
        return fname.strip('"')

    def _get_app_options(self):
        """Get commonly used configuration items from the config file"""

        self.input_directory = self.app_config.get_property("common",
                                                            "input_directory")
        self.output_directory = self.app_config.get_property(
            "common", "output_directory"
        )
        self.orography_dir = self.app_config.get_property("common", "orography_dir")
        #self.delete_processed = self.app_config.get_bool_property(
        #    "common", "delete_processed"
        #)
        #self.delete_source = self.app_config.get_bool_property(
        #    "common", "delete_source"
        #)
        self.variables_input = eval(self.app_config.get_property("common",
                                                                 "variables_input"))
        self.variables_rename = eval(self.app_config.get_property("common",
                                                                 "variables_rename"))
        self.delete_source = self.app_config.get_bool_property(
            "common", "delete_source"
        )
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
        self.um_archive_to_mass = self.app_config.get_bool_property("common",
                                                            "um_archive_to_mass")

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
        self.next_cycle = os.environ["NEXT_CYCLE"]
        self.tm2_cycle = os.environ["TM2_CYCLE"]
        self.ncodir = os.environ["NCODIR"]
