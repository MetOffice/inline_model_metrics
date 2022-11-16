# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import glob
import os
import subprocess
import sys
import shutil

from .tempest_common import (TempestExtremesAbstract, _is_date_after,
                             _is_date_after_or_equal)


class TempestError(Exception):
    """
    Custom TempestExtremes tracking exception
    """
    pass


class FrontalIdentification(TempestExtremesAbstract):
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
            self.processed_files_psl[regrid_resol] = ""
        for track_type in self.frontid_types:
            self.cmd_detect_type[track_type] = ""
        self.outdir = self.output_directory + "_" + "native"

        self.logger.debug(
            f"CYLC_TASK_CYCLE_TIME {self.cylc_task_cycle_time}, "
            f"runid {self.runid}"
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

        # Check whether the cylc date is after the most recent post-processed file
        #condition = self._get_tracking_date(timestamp_day)
        #if condition == "AlreadyComplete":
        #    self.logger.info(f"This date already completed {timestamp_day}")
        #    return

        # First write a dot_file to document which timestamps are yet to be tracked
        dot_file = "do_frontal_id"
        candidate_files = []
        self._write_dot_track_file(timestamp_day, timestamp_endday, dot_file=dot_file)

        dot_tracking_files = sorted(glob.glob(os.path.join(self.outdir, dot_file+"*")))
        self.logger.debug(f"dot_tracking_files {dot_tracking_files}")

        # loop through all dot files
        # only remove dot file when the detection has run
        # only try to do detection on data with timestamp before the current time
        # (want to make sure that data has been written on the previous step, hence
        # no conflict with writing if postproc on the next step might be running)
        if dot_tracking_files:
            for do_track_file in dot_tracking_files:
                ftimestamp_day = do_track_file.split(".")[1].split("-")[0]
                ftimestamp_endday = do_track_file.split(".")[1].split("-")[1]

                if self.inline_tracking == "True":
                    self.logger.debug(f"running inline {self.inline_tracking}")
                    # do not want to do calculations on data after or equal to the
                    # current cycle date, unless it is also the last
                    if _is_date_after_or_equal(ftimestamp_day, timestamp_day) \
                            and not self.is_last_cycle == "true":
                        continue

                    # if timestamp_previous is before the start date then no work
                    if _is_date_after(self.startdate, timestamp_previous):
                        continue

                # find the relevant input data using the given file pattern
                fname = self._file_pattern_processed(ftimestamp_day+"*", "*", "psl",
                                           frequency=self.data_frequency)
                file_search = os.path.join(self.input_directory, fname)
                self.logger.debug(f"file_search {file_search}")

                for regrid_resol in self.regrid_resolutions:
                    self.outdir = self.output_directory+'_'+regrid_resol
                    fname = self._file_pattern_processed(ftimestamp_day + "*", "*",
                                                         "ta_850",
                                                         frequency=self.data_frequency)
                    file_search = os.path.join(self.outdir, fname)
                    if glob.glob(file_search):
                        processed_files, variable_units = \
                        self._identify_processed_files(ftimestamp_day,
                                                       ftimestamp_endday,
                                                       grid_resol=regrid_resol)
                        #self.source_files[ftimestamp_day] = source_files
                        self.processed_files[ftimestamp_day] = processed_files
                        self.processed_files_psl[regrid_resol] = \
                            self.processed_files[ftimestamp_day]["ta_850"]
                        self.variable_units = variable_units

                        # run TempestExtremes detect blobs
                        candidate_file = self._run_detect_fronts(ftimestamp_day)

                        # process the AR netcdf files
                        #self._process_fronts_for_archive(candidate_file)

                        # if this timestep has worked OK, then need to remove
                        # the dot_file
                        self._remove_dot_track_file(ftimestamp_day, ftimestamp_endday,
                                                    dot_file=dot_file)

                        # at this point, I can delete the processed input data
                        #if self.delete_processed:
                        #    self._tidy_data_files(timestamp_previous, timestamp_day,
                        #                          self.variables_rename)

                    else:
                        self.logger.debug(f"no files to process for timestamp "
                                          f"{ftimestamp_day}")
        else:
            self.logger.error(f"no dot files to process ")

        # Test if new year, if so then concatenate all the previous year data
        # into 1 file
        is_new_year = (timestamp_day[0:4] != timestamp_previous[0:4]) and \
                _is_date_after(timestamp_previous, self.startdate)

        if is_new_year:
            self._produce_annual_fronts_file(self.outdir, timestamp_previous[0:4])

        if self.is_last_cycle == "true" or is_new_year:
            # archive any remaining fronts data
            self._archive_fronts_data(self.outdir, is_new_year, timestamp_previous[0:4])

    def _read_frontid_options(self, front_type):
        """
        Read the front_id command line parameters from the configuration.
        :param str front_type: The name of the type of tracking to run, possible
            values: detect
        :returns: A dictionary with keys of `detect`, and the values for each
            of these is
            a string containing the command line parameters for each of these
            front-id steps. The parameters are sorted into alphabetical
            order in each line.
        :rtype: dict
        """

        commands = {}
        fmt_value = {}
        for step in ["detect"]:
            try:
                step_config = self.app_config.section_to_dict(f"{front_type}_{step}")
                step_arguments = []
                for parameter in sorted(list(step_config.keys())):
                    if step_config[parameter]:
                        param_value = step_config[parameter]
                        step_arguments.append(f"--{parameter} {param_value}")
                commands[step] = " ".join(step_arguments)
                self.logger.debug(f"step, commands {step} {commands[step]}")
            except:
                commands[step] = ""
        return commands

    def _run_detect_fronts(self, timestamp):
        """
        Run the front identification.

        :param str timestamp: The timestep of data/tracking to process
        :returns: The path to the candidate files (as a list ordered by track
            type) and details of the processed input files (as a dict).
        :rtype: tuple
        :
        """
        os.chdir(os.path.dirname(self.frontid_detect_script))
        self.logger.debug(f"cwd {os.getcwd()}")

        for frontid_type in self.frontid_types:
            self.logger.debug(f"Running {frontid_type} detection")
            front_all_raw_file = os.path.join(
                self.outdir,
                f"{self.runid}_fronts_raw_all_{timestamp}_{frontid_type}.nc",
            )
            self.logger.debug(f"candidatefile {front_all_raw_file}")

            # Identify all raw fronts from model data
            # put together the command for all frontal identification
            cmd_io = "{} {} ".format(self.frontid_Rexec, self.frontid_detect_script)

            if "rh_850" in self.processed_files[timestamp]:
                cmd_io += " -r "+self.processed_files[timestamp]["rh_850"] + \
                    " --rname rh_850"
            elif "q_850" in self.processed_files[timestamp]:
                cmd_io += " -q " + self.processed_files[timestamp]["q_850"] + \
                    " --qname q_850"
            elif "wb_850" in self.processed_files[timestamp]:
                cmd_io += " -w " + self.processed_files[timestamp]["wb_850"] + \
                    " --wname wb_850"
            else:
                raise Exception("No humidity file available "+timestamp)

            if "ta_850" in self.processed_files[timestamp]:
                cmd_io += " -t "+self.processed_files[timestamp]["ta_850"] + \
                          " --tname ta_850"
            elif "wb_850" in self.processed_files[timestamp]:
                pass
            else:
                raise Exception("No temperature file available "+timestamp)

            if "ua_850" in self.processed_files[timestamp] and \
                    "va_850" in self.processed_files[timestamp]:
                cmd_io += " -u "+self.processed_files[timestamp]["ua_850"] + \
                          " -v " + self.processed_files[timestamp]["va_850"] + \
                          " --uname ua_850 --vname va_850"

            frontid_options = self._read_frontid_options(frontid_type)
            cmd_io += " " + frontid_options["detect"]

            cmd_io += " " + front_all_raw_file

            cmd_frontid = cmd_io
            self.logger.info(f"Front detect command {cmd_frontid}")

            sts = subprocess.run(
                cmd_frontid,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
            self.logger.debug(sts.stdout)
            if "EXCEPTION" in sts.stdout:
                msg = (
                    f"EXCEPTION found in front detect output\n"
                    f"{sts.stdout}"
                )
                raise RuntimeError(msg)

        return front_all_raw_file

    def _process_fronts_for_archive(
        self,
        fname,
        nc_compression="1",
        netcdf_format="4"
    ):
        """
        Process any AR files before archiving
        :param str fname: The filename to operate on
        :param str nc_compression: Compression level for netcdf files
        :param netcdf_format: Format of netcdf files (netcdf4 by default)
        """
        if os.path.exists(fname):
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

    def _produce_annual_fronts_file(
        self,
        outdir,
        year,
    ):
        """
        Merge one year of fronts data
        :param str outdir: output directory which contains AR files
        :param str year: the year of data to concatenate
        """
        files_to_join = sorted(glob.glob(os.path.join(outdir,
                                                      "*fronts_raw_all_"+year+"????_*.nc")))
        if len(files_to_join) > 0:
            file_year = files_to_join[0].replace(year+"0101", "year_"+year)
            cmd_io = "{} {} ".format(self.frontid_Rexec,
                os.path.join(os.path.dirname(self.frontid_detect_script),
                             "merge_fronts.R "))
            cmd = cmd_io + " " + " ".join(files_to_join) + \
                    " " + file_year

            self.logger.info(f"merge_fronts cmd {cmd}")

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
                for f in files_to_join:
                    os.remove(f)

    def _archive_fronts_data(self, outdir, is_new_year, year):
        """
        Archive the required AR data
        :param str outdir: output directory from which to archive AR files
        """
        if self.is_last_cycle == "true" :
            files_to_archive = glob.glob(os.path.join(outdir, "*fronts_raw_all*.nc"))
        elif is_new_year:
            files_to_archive = glob.glob(os.path.join(outdir,
                                                      "*fronts_raw_all*year_" + year + "*.nc"))
        self.logger.info(f"files_to_archive {files_to_archive}")

        if not os.path.exists(os.path.join(outdir, self._archived_files_dir)):
            os.makedirs(os.path.join(outdir, self._archived_files_dir))

        if len(files_to_archive) > 0:
            for front_file in files_to_archive:
                # move front_file to archive directory
                front_archive_file = os.path.join(outdir, self._archived_files_dir,
                                                os.path.basename(front_file))
                if self.is_last_cycle:
                    shutil.copy(front_file, front_archive_file)
                else:
                    os.replace(front_file, front_archive_file)
                    os.remove(front_file)
                self.logger.info(f"archived {front_file} to {front_archive_file}")

                with open(front_archive_file + ".arch", "a"):
                    os.utime(front_archive_file + ".arch", None)

    def _get_app_options(self):
        """Get commonly used configuration items from the config file"""

        super()._get_app_options()

        self.frontid_detect_script = self.app_config.get_property(
            "common", "frontid_detect_script"
        )
        self.frontid_Rexec = self.app_config.get_property(
            "common", "frontid_Rexec"
        )
        # track_types is a Python list so eval converts str to list
        self.frontid_types = eval(self.app_config.get_property("common",
                                                               "frontid_types"))
        self.variables_input = eval(self.app_config.get_property("common",
                                                        "frontid_variables_input"))
        self.variables_rename = eval(self.app_config.get_property("common",
                                                        "frontid_variables_rename"))
