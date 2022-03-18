# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import glob
import os
import subprocess
import sys

from .tempest_common import TempestExtremesAbstract


class TempestError(Exception):
    """
    Custom TempestExtremes tracking exception
    """
    pass


class TempestExtremesAR(TempestExtremesAbstract):
    """
    Run the TempestExtremes atmospheric river tracker inline with a climate
    model.
    """

    def __init__(self, arglist=None, **kwargs):
        super().__init__(version="1.0", **kwargs)
        self._parse_args(arglist, desc="Inline TempestExtreme tracking")
        self._parse_app_config()
        self._set_message_level()

        self.cmd_detectblobs = None

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
        for track_type in self.ar_types:
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
        dot_file = "do_ar_tracking"
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
                    fname = self._file_pattern_processed(ftimestamp_day + "*", "*", "viwve",
                                                         frequency=self.data_frequency)
                    file_search = os.path.join(self.outdir, fname)
                    if glob.glob(file_search):
                        processed_files, variable_units = \
                        self._identify_processed_files(ftimestamp_day, ftimestamp_endday,
                                            grid_resol=regrid_resol)
                        #self.source_files[ftimestamp_day] = source_files
                        self.processed_files[ftimestamp_day] = processed_files
                        self.processed_files_slp[regrid_resol] = \
                            self.processed_files[ftimestamp_day]["viwve"]
                        self.variable_units = variable_units

                        # run TempestExtremes detect blobs
                        candidate_files = self._run_detect_blobs(ftimestamp_day)
                        # if this timestep has worked OK, then need to remove the dot_file
                        # (the data is needed later)
                        self._remove_dot_track_file(timestamp_previous, timestamp_day)

                    else:
                        self.logger.debug(f"no files to process for timestamp "
                                          f"{ftimestamp_day}")
        else:
            self.logger.error(f"no dot files to process ")

        # at this point, I can delete the input data for the T-2 timestep
        if self.delete_processed:
            self._tidy_data_files(timestamp_tm2, timestamp_tm1)

        # if self.delete_source:
        #    self._tidy_data_files(timestamp_previous, timestamp_day,
        #    f_remove = 'source')

    def _run_detect_blobs(self, timestamp):
        """
        Run the Tempest blob detection.

        :param str timestamp: The timestep of data/tracking to process
        :returns: The path to the candidate files (as a list ordered by track
            type) and details of the processed input files (as a dict).
        :rtype: tuple
        :
        """
        self.logger.debug(f"cwd {os.getcwd()}")

        for ar_type in self.ar_types:
            self.logger.debug(f"Running {ar_type} detection")
            candidatefile = os.path.join(
                self.outdir,
                f"{self.um_runid}_AR_{timestamp}_{ar_type}.txt",
            )
            self.logger.debug(f"candidatefile {candidatefile}")

            #cmd_io = '{} --out {} '.format(
            #        self.ar_detectblobs_script,
            #        candidatefile)
            cmd_io = "{} ".format(self.ar_detectblobs_script)

            # need the first input file not to be orography, since that file has to
            # have a time coordinate
            fnames = []
            #for key in self.processed_files[timestamp]:
            for key in ["viwve", "viwvn"]:
                if isinstance(self.processed_files[timestamp][key], list):
                    fnames.extend(self.processed_files[timestamp][key])
                else:
                    fnames.append(self.processed_files[timestamp][key])

            in_file_list = os.path.join(self.outdir, 'in_file_list_detectblobs.txt')
            with open(in_file_list, "w") as fh:
                text_str = ';'.join(fnames)
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)

            fname = self.processed_files[timestamp]["viwve"]
            ARmask = fname.replace("viwve", "ARmask")

            cmd_io += '--in_data_list '+in_file_list+' --out '+ARmask+' '

            tracking_phase_commands = self._construct_command(ar_type)
            cmd_detectblobs = cmd_io + tracking_phase_commands["detectblobs"]
            self.cmd_detect_type[ar_type] = cmd_detectblobs
            self.logger.info(f"Detect command {cmd_detectblobs}")

            sts = subprocess.run(
                cmd_detectblobs,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
            self.logger.debug(sts.stdout)
            if "EXCEPTION" in sts.stdout:
                msg = (
                    f"EXCEPTION found in TempestExtreme detectblobs output\n"
                    f"{sts.stdout}"
                )
                raise RuntimeError(msg)

        return ARmask

    def _get_app_options(self):
        """Get commonly used configuration items from the config file"""

        super()._get_app_options()

        self.ar_detectblobs_script = self.app_config.get_property(
            "common", "ar_detectblobs_script"
        )
        # track_types is a Python list so eval converts str to list
        self.ar_types = eval(self.app_config.get_property("common", "ar_types"))
        self.variables_input = eval(self.app_config.get_property("common",
                                                                 "ar_variables_input"))
        self.variables_rename = eval(self.app_config.get_property("common",
                                                                 "ar_variables_rename"))
