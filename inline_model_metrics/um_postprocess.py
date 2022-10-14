# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import glob
import os
import re
import shutil
import subprocess
import sys

#from afterburner.apps import AbstractApp

from .tempest_common import (TempestExtremesAbstract)

class TempestError(Exception):
    """
    Custom TempestExtremes tracking exception
    """
    pass

class UMTempestPostprocess(TempestExtremesAbstract):
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
            f"runid {self.runid}"
        )

        timestamp_day = self.cylc_task_cycle_time[:8]
        timestamp_previous = self.previous_cycle[:8]
        timestamp_tm2 = self.tm2_cycle[:8]

        # Check whether the cylc date is after the most recent post-processed file
        condition = self._get_tracking_date(timestamp_day)
        if condition == "AlreadyComplete":
            return

        # this section of code processes data from the current timestep
        current_time = timestamp_day

        for regrid_resol in self.regrid_resolutions:
            self.outdir = self.output_directory + '_' + regrid_resol
            if self.um_archive_to_mass:
                # run the archiving on any .arch files that exist
                self._archive_tracking_to_mass(os.path.join(self.outdir,
                                                        self._archived_files_dir))
            if self.delete_processed:
                if not self.is_last_cycle == "true":
                    self._tidy_data_files(timestamp_tm2,
                            timestamp_previous, self.variables_rename)

        # delete source data if required
        if self.delete_source:
            for var in self.variables_input:
                fname = self._file_pattern(timestamp_tm2 + "*", "*", var,
                                       stream=self.um_stream, frequency="*")
                file_name = os.path.join(self.input_directory, fname)
                files_exist = glob.glob(file_name)
                if len(files_exist) > 0:
                    for f in files_exist:
                        self.logger.info(f"deleting source file_name {f}")
                        os.remove(f)

    def _archive_tracking_to_mass(
        self,
        archive_directory
    ):
        """
        Archive any files with .arch to Met Office MASS system
        :param str archive_directory: The directory to look for any .arch files
        """

        #moosedir = "moose:/crum/{}/{}/"
        archive_files = glob.glob(os.path.join(archive_directory, "*.arch"))
        if len(archive_files) > 0:
            if self.ensemble == "":
                moo_path = self.moo_dir.format(self.suiteid)
            else:
                moo_path = self.moo_dir.format(self.suiteid, self.ensemble)

            for fname_arch in archive_files:
                archive_error = False
                fname = fname_arch[:-5]
                if fname[-2:] == "nc":
                    mass_stream = "any.nc.file"
                else:
                    mass_stream = "ady.file"
                if os.stat(fname).st_size == 0:
                    self.logger.debug(f"File is zero length, no archive {fname}")
                    return

                cmd = "moo put -F " + fname + " " + moo_path + "/" + mass_stream
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
                    archive_error = True
                else:
                    os.remove(fname)
                    os.remove(fname_arch)
                self.logger.debug(sts.stdout)

                if archive_error:
                    if os.path.isdir(self.backup_data_directory):
                        backup_data_path = os.path.join(self.backup_data_directory,
                                                    self.suiteid)
                        if not os.path.exists(backup_data_path):
                                os.makedirs(backup_data_path)
                        cmd = "cp " + fname + " " + backup_data_path
                        self.logger.debug(f"Archive cmd {cmd}")
                    else:
                        self.logger.debug(f"Archive directory does not exist " + \
                                          f"{self.backup_data_directory}")
                    raise RuntimeError(msg)

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
        self.um_archive_to_mass = self.app_config.get_bool_property("common",
                                                            "um_archive_to_mass")
        try:
            self.um_stream = eval(self.app_config.get_property("common",
                                                            "um_stream"))
        except:
            self.um_stream = "pt"

        try:
            self.backup_data_directory = self.app_config.get_property("common",
                                                            "backup_directory")
        except:
            self.backup_data_directory = ""

        try:
            self.ensemble = eval(self.app_config.get_property("common",
                                                            "ensemble"))
        except:
            self.ensemble = ""

        if self.ensemble == "":
            self.moo_dir = 'moose:/crum/{}'
        else:
            self.moo_dir = 'moose:/ens/{}/{}'

    # def _get_environment_variables(self):
    #     """
    #     Get the required environment variables from the suite. A list and
    #     explanation of the required environment variables is included in the
    #     documentation.
    #     """
    #     try:
    #         self.suiteid = os.environ["SUITEID_OVERRIDE"]
    #     except:
    #         self.suiteid = os.environ["CYLC_SUITE_NAME"]
    #     self.runid = self.suiteid.split('-')[1]
    #     self.cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
    #     self.next_cycle = os.environ["NEXT_CYCLE"]
    #     self.previous_cycle = os.environ["PREVIOUS_CYCLE"]
    #     self.tm2_cycle = os.environ["TM2_CYCLE"]
    #     self.ncodir = os.environ["NCODIR"]
    #     self.is_last_cycle = os.environ["IS_LAST_CYCLE"]
    #     self.inline_tracking = os.environ["INLINE_TRACKING"]
