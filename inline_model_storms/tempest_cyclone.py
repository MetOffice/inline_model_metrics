# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import glob
import os
import subprocess
import sys
import shutil
import pandas as pd
import xarray as xr
import csv
from datetime import datetime

from netCDF4 import Dataset

from .tempest_common import (TempestExtremesAbstract, _is_date_after,
                             _is_date_after_or_equal)

from .load_trajectories import get_trajectories,get_trajectories_csv
from .save_trajectories import save_trajectories_netcdf
from .trajectory_manipulations import convert_date_to_step,fill_trajectory_gaps

from tempest_helper import (
    count_trajectories,
    #get_trajectories,
    plot_trajectories_cartopy,
    #save_trajectories_netcdf,
    storms_overlap_in_time,
    storm_overlap_in_space,
    write_track_line,
    remove_duplicates_from_track_files
)


class TempestError(Exception):
    """
    Custom TempestExtremes tracking exception
    """
    pass


class TempestExtremesCyclone(TempestExtremesAbstract):
    """
    Run the TempestExtremes cyclone tracker inline with a climate model.
    """

    def __init__(self, arglist=None, **kwargs):
        super().__init__(version="1.0", **kwargs)
        self._parse_args(arglist, desc="Inline TempestExtremes TC tracking")
        self._parse_app_config()
        self._set_message_level()

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

        self.file_format = os.path.join(
                "{outdir}", "{runid}_{app}_{period}{timestamp}-{timestamp_next}_{data_freq}_{track_type}.txt")
        self.file_format_nc = os.path.join(
                "{outdir}", "{runid}_{app}_{period}{timestamp}-{timestamp_next}_{data_freq}_{track_type}.nc")

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
        for track_type in self.track_types:
            self.cmd_detect_type[track_type] = ""
            self.cmd_stitch_type[track_type] = ""
        self.outdir = self.output_directory + "_" + "native"

        self.logger.debug(
            f"CYLC_TASK_CYCLE_TIME {self.cylc_task_cycle_time}, "
            f"runid {self.runid}"
        )

        timestamp_current = self.cylc_task_cycle_time[:8]
        timestamp_next = self.next_cycle[:8]
        startdate = self.startdate[:8]
        enddate = self.enddate[:8]
        timestamp_tm2 = self.tm2_cycle[:8]
        end_of_year = False
        if timestamp_current[:4] != timestamp_next[:4]:
            end_of_year = True
        self.logger.debug(f"timestamp_day at top {timestamp_current}")
        self.logger.debug(f"startdate {self.startdate}")
        self.logger.debug(f"enddate {self.enddate}")
        self.logger.debug(f"lastcycle {self.lastcycle}")
        self.logger.debug(f"previous_cycle {self.previous_cycle}")
        self.logger.debug(f"next_cycle {self.next_cycle}")
        self.logger.debug(f"tm2_cycle {timestamp_tm2}")
        self.logger.debug(f"is_last_cycle {self.is_last_cycle}")
        self.logger.debug(f"inline_tracking {self.inline_tracking}")
        self.logger.debug(f"track_by_year {self.track_by_year}")
        self.logger.debug(f"track_at_end {self.track_at_end}")

        if self.inline_tracking == "True":
            self.logger.debug(f"running inline {self.inline_tracking}")

        for regrid_resol in self.regrid_resolutions:
            self.outdir = self.output_directory+'_'+regrid_resol
            fname = self._file_pattern_processed(timestamp_current + "*", "*",
                                                "psl", frequency=self.data_frequency)
            file_search = os.path.join(self.outdir, fname)
            self.logger.debug(f"file_search {file_search}")
            if glob.glob(file_search):
                processed_files, variable_units = \
                    self._identify_processed_files(timestamp_current,
                            timestamp_next, grid_resol=regrid_resol)
                self.processed_files[timestamp_current] = processed_files
                self.processed_files_psl[regrid_resol] = \
                        self.processed_files[timestamp_current]["psl"]
                self.variable_units = variable_units

                # run TempestExtremes detection
                detect_file = self._run_detection(self.outdir, timestamp_current, timestamp_next)
                vor_proc1 = self._run_variable_processor(self.outdir, timestamp_current, timestamp_next, proc='varproc1')
                vor_proc2 = self._run_variable_processor(self.outdir, timestamp_current, timestamp_next, proc='varproc2')
                detectblobs_file = self._run_detectblobs(self.outdir, timestamp_current, timestamp_next, detect_file)
            else:
                self.logger.debug(f"No psl files found for detection {file_search}")

            self._current_year_value = timestamp_current[0:4]

            self.logger.debug(f"track_by_year {self.track_by_year}")

            # Run the tracking if this is the last file of the year or 
            # the end of the run
            if self.track_by_year and \
                      (end_of_year or self.is_last_cycle):
                self._run_tracking_and_editing(
                    self.outdir,
                    timestamp_current,
                    timestamp_next,
                    self.processed_files_psl[regrid_resol],
                    track_at_end=False,
                    is_lastcycle=self.is_last_cycle)

                self._run_blobtracking_and_stats(
                    self.outdir,
                    timestamp_current,
                    timestamp_next,
                    track_at_end=False,
                    is_lastcycle=self.is_last_cycle)

            # Run the tracking over whole period if this is the end of the run
            if self.track_at_end and self.is_last_cycle:
                self._run_tracking_and_editing(
                    self.outdir,
                    timestamp_current,
                    timestamp_next,
                    self.processed_files_psl[regrid_resol],
                    track_at_end=True,
                    is_lastcycle=self.is_last_cycle,
                    startdate=startdate,
                    enddate=enddate)

                self._run_blobtracking_and_stats(
                    self.outdir,
                    timestamp_current,
                    timestamp_next,
                    track_at_end=True,
                    is_lastcycle=self.is_last_cycle,
                    startdate=startdate,
                    enddate=enddate)

            self._archive_track_data(
                self.outdir,
                timestamp_current,
                timestamp_next)

            self._tidy_files(self.outdir,
                             timestamp_tm2[:8],
                             self.previous_cycle[:8])

    def _run_tracking_and_editing(
        self,
        outdir,
        timestamp_current,
        timestamp_next,
        nc_file,
        track_at_end=False,
        is_lastcycle=False,
        startdate=None,
        enddate=None
    ):
        """
        Prepare for tracking by setting up required filenames and calling the
        Tempest tracking (stitch)

        :param str timestamp_current: The current timestep of data/tracking to process
        :param bool track_at_end: If this is this the last cycle in the 
        :     model run, then do operations over whole period
        :param bool is_lastcycle: Is this the last cycle in the model run
        """

        # run the stitching and edit/nodeedit on native grid data
        # if we want to track by year, then we want to collect candidates from (if existing) Dec of previous year, Jan-Dec of this year, and then (after tracking) need to remove any tracks that finish during Dec of previous year, that should sort out the year overlap
        # timestamp_day tells us year and month - if month is 12, then this is last data of year, so can trigger end of year
        detect_files = self._find_detect_files(outdir, timestamp_current, track_at_end, is_lastcycle, startdate=startdate, enddate=enddate)

        tracked_files = self._run_stitch_alltimes(outdir, timestamp_current, timestamp_next, detect_files, nc_file)
        self.logger.debug(f"Call NodeFileEditor {tracked_files}")

        nodeeditor_files = self._run_nodeeditor_alltimes(outdir, timestamp_current, tracked_files, detect_files)

    def _run_blobtracking_and_stats(
        self,
        outdir,
        timestamp_current,
        timestamp_next,
        track_at_end=False,
        is_lastcycle=False,
        startdate=None,
        enddate=None
    ):
        """
        Prepare for tracking by setting up required filenames and calling the
        Tempest tracking (stitch)

        :param str timestamp_current: The current timestep of data/tracking to process
        :param bool track_at_end: If this is this the last cycle in the 
        :     model run, then do operations over whole period
        :param bool is_lastcycle: Is this the last cycle in the model run
        """

        # run the stitching and edit/nodeedit on native grid data
        # if we want to track by year, then we want to collect candidates from (if existing) Dec of previous year, Jan-Dec of this year, and then (after tracking) need to remove any tracks that finish during Dec of previous year, that should sort out the year overlap
        # timestamp_day tells us year and month - if month is 12, then this is last data of year, so can trigger end of year
        detectblobs_files = self._find_detectblobs_files(outdir, timestamp_current, track_at_end, is_lastcycle, startdate=startdate, enddate=enddate)

        blobstats_files = self._run_blobstats(outdir, timestamp_current, detectblobs_files)

        trackedblobs_files = self._run_stitchblobs_alltimes(outdir, timestamp_current, detectblobs_files)

    def _archive_track_data(self, outdir, timestamp, timestamp_next,
                            track_at_end=False):
        """
        Set up archiving of the required data by producing .arch file beside any
        file needing to be archived, and um_postprocess will do archiving
        Want to do year files, last cycle. May want to move old files out of way.
        Need to figure out when to add .arch, and how to note that a file has been archived already (so don't add another .arch). Could choose to only do this on last cycle (if not year split), but risks all data being lost. Maybe at end of year and last cycle.
        If file is archived, do we have to delete it? This determines when to note file needs archiving.
        Need to archive every candidate file and detectblobs file as soon as possible (to avoid loosing), so that tracking could be done later anyway. Hence the tarring of multiple files rather gets in way. But files small, so not ideal for MASS.
        Naming of track files - only if exactly the same time period will the archiving overwrite an existing file.
        :param str outdir: output directory
        """
        if not os.path.exists(os.path.join(outdir, self._archived_files_dir)):
            os.makedirs(os.path.join(outdir, self._archived_files_dir))

        tracked_types = ["track", "blobstats"]
        for track_type in self.track_types:
            if self._construct_command(track_type)["nodeedit"] is not None:
                trackname = "tracknodeedit"
            else:
                trackname = "track"

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detect", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            self.logger.debug(f"Archive, is new year? {self._is_new_year}")
            # archive the annual files or remaining files at end of run
            if self._is_new_year or track_at_end:
                # first move (if new year) or copy (if end of run) the files to the
                # archive directory, and then add .arch archive indicator
                if self._is_new_year:
                    year = self._old_year_value
                    file_command = "move"
                else:
                    year = self._current_year_value
                    file_command = "copy"
                subdir = outdir.split('/')[-1]

            detect_files_search = self.file_format.format(
                outdir=outdir,
                app='detect',
                runid=self.runid,
                period="",
                timestamp="????????", 
                timestamp_next="????????",
                data_freq=timefilter, 
                track_type=track_type)

            detectblobs_files_search = self.file_format_nc.format(
                outdir=outdir, 
                app='detectblobs',
                runid=self.runid, 
                period="", 
                timestamp="????????", 
                timestamp_next="????????", 
                data_freq=timefilter, 
                track_type=track_type)
            
            files = sorted(glob.glob(detect_files_search))
            self.logger.debug(f"Files to archive {detect_files_search} {files}")
            if len(files) > 0:
                for f in files:
                    archive_path = os.path.join(
                        os.path.dirname(f), self._archived_files_dir,
                        os.path.basename(f))
                    f_archived = archive_path+".archived"
                    f_archive = archive_path+".arch"
                    if not os.path.exists(f_archived):
                        with open(f_archive, "a"):
                            os.utime(f_archive, None)
                    else:
                        with open(f_archived, "a"):
                            os.utime(f_archived, None)
                        
            files = sorted(glob.glob(detectblobs_files_search))
            self.logger.debug(f"Files to archive {detectblobs_files_search} {files}")
            if len(files) > 0:
                for f in files:
                    archive_path = os.path.join(
                        os.path.dirname(f), self._archived_files_dir,
                        os.path.basename(f))
                    f_archived = archive_path+".archived"
                    f_archive = archive_path+".arch"
                    if not os.path.exists(f_archived):
                        with open(f_archive, "a"):
                            os.utime(f_archive, None)
                    else:
                        with open(f_archived, "a"):
                            os.utime(f_archived, None)

            for ttype in tracked_types:
                tracked_file = os.path.join(
                    outdir,
                    f"{self.runid}_{ttype}*_????????-????????_*_{track_type}"
                )
                self.logger.debug(f"track files to archive {tracked_file}")
                
                for f_ending in [".txt", ".csv", ".nc"]:
                    files = glob.glob(tracked_file+f_ending)
                    if len(files) > 0:
                        for f in files:
                            archive_path = os.path.join(
                                os.path.dirname(f), self._archived_files_dir,
                                os.path.basename(f))
                            f_archived = archive_path+".archived"
                            f_archive = archive_path+".arch"
                            if not os.path.exists(f_archived):
                                with open(f_archive, "a"):
                                    os.utime(f_archive, None)
                            else:
                                with open(f_archived, "a"):
                                    os.utime(f_archived, None)

    def _tidy_files(self, outdir, timestamp_tm2, timestamp_previous):
        """
        Find all files that will not be needed on next step, and remove them

        :param str outdir: directory for the files
        :param str timestamp: The current timestep of data/tracking to process
        """

        search = os.path.join(outdir, "*period*")
        period_files = glob.glob(search)
        if len(period_files) > 0:
            for f in period_files:
                os.remove(f)

        nodeedit_included = False
        for track_type in self.track_types:
            tracking_phase_commands = self._construct_command(track_type)["nodeedit"]
            if tracking_phase_commands is not None:
                nodeedit_included = True

        if nodeedit_included:
            variables_to_delete = self.variables_rename.copy()
            for var in self.nodeedit_vars:
                variables_to_delete.remove(var)
        else:
            variables_to_delete = self.variables_rename.copy()

        self._tidy_data_files(timestamp_tm2,
                              timestamp_previous, variables_to_delete,
                              f_remove=self.delete_processed)

        

    def _diff_month(self, d1, d2):
        return (d1.year - d2.year) * 12 + d1.month - d2.month

    def _find_detect_files(self, outdir, timestamp, track_at_end, 
                           is_lastcycle, startdate=None, enddate=None):
        """
        Find all candidate files required for stitching, based on year
        and whether to track by year or do whole period
        Ideally for the track_at_end, one would calculate the number
        of expected files and check it is correct

        :param str outdir: directory for the files
        :param str timestamp: The current timestep of data/tracking to process
        :param bool track_at_end: Do tracking at end of run over all files
        :param bool is_lastcycle: Is this the last cycle in the model run
        :returns: dictionary of candidates with keys of track_types
        :rtype: list
        """

        detects_by_tracktype = {}
        year = timestamp[:4]
        yearm1 = str(int(year)-1)
        self.logger.debug(f"find detects {year} {yearm1}")
        for track_type in self.track_types:

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detect", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            # On last cycle, find all candidates
            if track_at_end:
                detectfile_search = self.file_format.format(
                    outdir=outdir, 
                    runid=self.runid, 
                    app="detect",
                    period="", 
                    timestamp="????????", 
                    timestamp_next="????????", 
                    data_freq=timefilter, 
                    track_type=track_type)
                detects_all = sorted(glob.glob(detectfile_search))
                if enddate is not None:
                    detects = []
                    for d in detects_all:
                        sdate = os.path.basename(d).split("_")[2][:8]
                        edate = os.path.basename(d).split("_")[2].split("-")[1][:8]
                        if sdate >= startdate and edate <= enddate:
                            detects.append(d)
                    n_months_expected = self._diff_month(datetime(int(enddate[:4]), int(enddate[4:6]), int(enddate[6:])), datetime(int(startdate[:4]), int(startdate[4:6]), int(startdate[6:])))
                    self.logger.debug(f"n_months_expected {n_months_expected}, len(detects) {len(detects)}")
                else:
                    detects = detects_all
            # track_by_year
            else:
                # Find candidates for this year
                detectfile_search = self.file_format.format(
                    outdir=outdir, 
                    runid=self.runid, 
                    app="detect", 
                    period="", 
                    timestamp=year+"????", 
                    timestamp_next="????????", 
                    data_freq=timefilter, 
                    track_type=track_type)
                self.logger.debug(f"detects_thisyear_search {detectfile_search}")
                detects_thisyear = sorted(glob.glob(detectfile_search))
                first_month = int(os.path.basename(detects_thisyear[0]).split("_")[2].split("-")[0][4:6])
                last_month = int(os.path.basename(detects_thisyear[-1]).split("_")[2].split("-")[0][4:6])

                detectfile_search = self.file_format.format(
                    outdir=outdir, 
                    runid=self.runid, 
                    app="detect", 
                    period="",
                    timestamp=yearm1+"12??", 
                    timestamp_next="????????", 
                    data_freq=timefilter, 
                    track_type=track_type)
                detects_prev_dec = glob.glob(detectfile_search)
                self.logger.debug(f"detects_lastyear_search {detectfile_search}")
                self.logger.debug(f"detects_thisyear {len(detects_thisyear)} {detects_thisyear}")
                self.logger.debug(f"detects_prev_dec {detects_prev_dec}")
                detects = detects_thisyear
                if len(detects_prev_dec) == 1:
                    detects1 = detects_prev_dec + detects_thisyear
                    detects = detects1

                if year == self.startdate[:4]:
                    # this is the first year, so may not have 12 months
                    expected_months = int(timestamp[4:6]) - int(self.startdate[4:6]) + 1
                else:
                    expected_months = int(timestamp[4:6])

                # want to check that we have the right number of files
                n_files_expected = expected_months + len(detects_prev_dec)

                # if this is not the last cycle, then need at least 12 months
                if not is_lastcycle and len(detects) < n_files_expected:
                    raise Exception("Not correct number of detects for full year "+str(len(detects)))

                if n_files_expected != len(detects):
                    self.logger.debug(f"detects {detects}, len(detects) {len(detects)}, n_files_expected {n_files_expected}")
                    raise Exception("Not correct number of detects for this period "+len(detects))

            if len(detects) == 0:
                raise Exception("No detect files found "+track_type)

            detects_by_tracktype[track_type] = detects
            self.logger.debug(f"detects {detects}")
        return detects_by_tracktype

    def _find_detectblobs_files(self, outdir, timestamp, track_at_end, 
                                is_lastcycle, startdate=None, enddate=None):
        """
        Find all detect files required for blobstitching, based on year
        and whether to track by year or do whole period

        :param str outdir: directory for the files
        :param str timestamp: The current timestep of data/tracking to process
        :param bool track_at_end: If this is this the last cycle in the 
        :     model run, then do operations over whole period
        :param bool is_lastcycle: Is this the last cycle in the model run
        :returns: dictionary of detects with keys of track_types
        :rtype: list
        """

        detectblobs_by_tracktype = {}
        year = timestamp[:4]
        yearm1 = str(int(year)-1)
        self.logger.debug(f"find detectblobs {year} {yearm1}")
        for track_type in self.track_types:
            tracking_phase_commands = self._construct_command(track_type)["detectblobs"]
            if tracking_phase_commands is None:
                continue

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detectblobs", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            # On last cycle, find all candidates
            if track_at_end:
                blobs_search = self.file_format_nc.format(
                    outdir=outdir, 
                    app="detectblobs", 
                    runid=self.runid, 
                    period="", 
                    timestamp="????????", 
                    timestamp_next="????????", 
                    data_freq=timefilter, 
                    track_type=track_type)
                blobs_all = sorted(glob.glob(blobs_search))
                if enddate is not None:
                    blobs = []
                    for d in blobs_all:
                        sdate = os.path.basename(d).split("_")[2][:8]
                        edate = os.path.basename(d).split("_")[2].split("-")[1][:8]
                        if sdate >= startdate and edate <= enddate:
                            blobs.append(d)
                    n_months_expected = self._diff_month(datetime(int(enddate[:4]), int(enddate[4:6]), int(enddate[6:])), datetime(int(startdate[:4]), int(startdate[4:6]), int(startdate[6:])))
                    self.logger.debug(f"n_months_expected {n_months_expected}, len(blobs) {len(blobs)}")
                else:
                    blobs = blobs_all
            else:
                # Find candidates for this year
                blobs_search = self.file_format_nc.format(
                    outdir=outdir, 
                    app="detectblobs", 
                    runid=self.runid, 
                    period="", 
                    timestamp=year+"????", 
                    timestamp_next="????????", 
                    data_freq=timefilter, 
                    track_type=track_type)
                blobs_thisyear = sorted(glob.glob(blobs_search))
                first_month = int(os.path.basename(blobs_thisyear[0]).split("_")[2].split("-")[0][4:6])
                last_month = int(os.path.basename(blobs_thisyear[-1]).split("_")[2].split("-")[0][4:6])

                blobs_search = self.file_format_nc.format(
                    outdir=outdir, 
                    app="detectblobs", 
                    runid=self.runid, 
                    period="", 
                    timestamp=yearm1+"12??", 
                    timestamp_next="????????", 
                    data_freq=timefilter, 
                    track_type=track_type)
                blobs_prev_dec = glob.glob(blobs_search)

                self.logger.debug(f"blobs_thisyear {blobs_thisyear}")
                self.logger.debug(f"blobs_prev_dec {blobs_prev_dec}")
                blobs = blobs_thisyear
                if len(blobs_prev_dec) == 1:
                    blobs1 = blobs_prev_dec + blobs_thisyear
                    blobs = blobs1

                if year == self.startdate[:4]:
                    # this is the first year, so may not have 12 months
                    expected_months = int(timestamp[4:6]) - int(self.startdate[4:6]) + 1
                else:
                    expected_months = int(timestamp[4:6])

                # want to check that we have the right number of files
                n_files_expected = expected_months + len(blobs_prev_dec)

                if not is_lastcycle and len(blobs) < n_files_expected:
                    raise Exception("Not correct number of blobs for full year "+str(len(blobs)))

                if n_files_expected != len(blobs):
                    self.logger.debug(f"blobs {blobs}, len(blobs) {len(blobs)}, n_files_expected {n_files_expected}")
                    raise Exception("Not correct number of blobs for this period "+str(len(blobs)))
            if len(blobs) == 0:
                raise Exception("No detectblobs files found "+track_type)

            detectblobs_by_tracktype[track_type] = blobs
            self.logger.debug(f"detectblobs {blobs}")
        return detectblobs_by_tracktype

    def _remove_overlapping_tracks(self, trackfile, year):
        """
        Remove any tracks that end before the year for which these tracks 
        are valid

        :param str tracks: csv file containing tracks to be edited
        :param str year: Main year for these tracks (assumed Jan-Dec)
        :returns: edited tracks containing only those that exist in year
        """
        trackfile_amend = trackfile[:-4]+'_amend.csv'
        tracks = pd.read_csv(trackfile, skipinitialspace=True)
        tids = pd.unique(tracks["track_id"])
        self.logger.debug(f"track tids {tids}")
        year_target = int(year)
        yearm1 = year_target - 1

        # find tracks that end before the start of the current year
        tids_to_drop = []
        for tid in tids:
            track = tracks[tracks["track_id"] == tid]
            track_startyear = track.year.values[0]
            track_endyear = track.year.values[-1]
            if (track_startyear == yearm1 and track_endyear == yearm1):
                tids_to_drop.append(tid)

        # drop the tracks that end before the current year
        tracks_amend = tracks.copy()
        for tid in tids_to_drop:
            tracks_amend = tracks_amend.drop(tracks_amend[tracks_amend["track_id"] == tid].index)

        # reset the track_id values to between 0-ntracks
        tracks_zero = tracks_amend.copy()
        tids_new = pd.unique(tracks_amend["track_id"])
        self.logger.debug(f"track tids_new {tids_new}")
        for iz, tid in enumerate(tids_new):
            self.logger.debug(f"tid new {iz} {tid}")
            track = tracks_zero[tracks_zero["track_id"] == tid]
            tracks_zero.loc[track.index[0]:track.index[-1], 'track_id'] = iz

        tracks_zero.to_csv(trackfile_amend, index=False)

        # replace the amended file back into the original trackfile
        os.replace(trackfile_amend, trackfile)

    def _run_detection(self, outdir, timestamp, timestamp_next):
        """
        Run the Tempest detection.

        :param str timestamp: The timestep of data/tracking to process
        :param str timestamp_next: The next timestep of data/tracking to process
        :returns: The path to the detect files (as a list ordered by track
            type) and details of the processed input files (as a dict).
        :rtype: tuple
        :
        """
        self.logger.debug(f"cwd {os.getcwd()}")

        detect_files = []

        for track_type in self.track_types:
            self.logger.debug(f"Running {track_type} detection")

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detect", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            detectfile = self.file_format.format(
                outdir=outdir, 
                runid=self.runid, app="detect", 
                period="", timestamp=timestamp, 
                timestamp_next=timestamp_next, data_freq=timefilter,
                track_type=track_type)
            self.logger.debug(f"detectfile {detectfile}")

            cmd_io = '{} --out {} '.format(
                    self.tc_detect_script,
                    detectfile)

            # need the first input file not to be orography, since that file has to
            # have a time coordinate
            fnames = []
            for key in self.processed_files[timestamp]:
                if isinstance(self.processed_files[timestamp][key], list):
                    fnames.extend(self.processed_files[timestamp][key])
                else:
                    fnames.append(self.processed_files[timestamp][key])
            fnames.remove(self.processed_files[timestamp]['orog'])

            in_file_list = os.path.join(self.outdir, 'in_file_list_detect.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(fnames) + ';' + \
                           self.processed_files[timestamp]['orog']
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
                    f"EXCEPTION found in TempestExtreme detect output\n"
                    f"{sts.stdout}"
                )
                raise RuntimeError(msg)

            detect_files.append(detectfile)

        return detect_files

    def _stitch_file(self, detect_file, tracked_file, track_type):
        """
        Run the Tempest stitching on the specified file.

        :param str detect_file: The path of the detect input file.
        :param str tracked_file: The path of the stitched output file.
        :param str track_type: The name of the type of tracking required.
        """

        tracking_phase_commands = self._construct_command(track_type)

        self.logger.debug(f"tracked_file {tracked_file}")

        # stitch detects together
        cmd_stitch_io = "{} --in {} --out {} ".format(
            self.tc_stitch_script, detect_file, tracked_file
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

    def _run_stitch_alltimes(self, outdir, timestamp, timestamp_end, detect_files, nc_file):
        """
        Concatenate the detect files together and then
        stitch this file.

        :param str outdir: directory for the files
        :param str timestamp: The current timestep of data/tracking to process
        :param str timestamp: The previous timestep of data/tracking to process
        :returns: The string paths to the annual stitched file in a list.
        :rtype: list
        """

        tracked_files = {}
        year = timestamp[:4]
        yearm1 = str(int(year)-1)
        month = timestamp[4:6]
        for track_type in self.track_types:
            detects = detect_files[track_type]
            self.logger.debug(f"detects {detects}")

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detect", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            # find out format of output file and name track file accordingly
            fileformat = self._command_parameter(
                track_type, 
                "stitch", 
                "out_file_format")
            if fileformat is None or fileformat == "":
                fileformat = "gfdl"
            self.logger.debug(f"File format of stitch file {fileformat}")

            if len(detects) >= 1:
                start_date = os.path.basename(detects[0]).split("_")[2].split('-')[0]
                end_date = os.path.basename(detects[-1]).split("_")[2].split('-')[1]
                end_month = end_date[4:6]
                if end_month == "12":
                    self.logger.debug(f"end of year {end_date} {end_month}")

                detect_alltimes = self.file_format.format(
                    outdir=outdir, 
                    app="detect",
                    runid=self.runid, 
                    period="period_", 
                    timestamp=start_date, 
                    timestamp_next=end_date, 
                    data_freq=timefilter, 
                    track_type=track_type)

                with open(detect_alltimes, "w") as out_file:
                    for detect_file in detects:
                        with open(detect_file) as in_file:
                            for line in in_file:
                                out_file.write(line)

                track_file_alltimes = self.file_format.format(
                    outdir=outdir, 
                    app="track",
                    runid=self.runid, 
                    period="", 
                    timestamp=start_date, 
                    timestamp_next=end_date, 
                    data_freq=timefilter, 
                    track_type=track_type)[:-3]+fileformat.strip('\"')

                self._stitch_file(detect_alltimes, track_file_alltimes, track_type)
                #if self.track_by_year:
                #    self._remove_overlapping_tracks(track_alltimes, year)

                tracked_files[track_type] = track_file_alltimes

                # only run conversion to netcdf if file is csv
                if "csv" in fileformat:
                    self._read_write_and_plot_tracks(
                        track_file_alltimes,
                        track_type,
                        nc_file,
                        start_date,
                        end_date,
                        self.variable_units,
                        title_prefix=f"{self.runid} {self.resolution_code} "
                        f"{start_date}_{end_date}",
                        title_suffix=f"{track_type} tracks",
                        write_to_netcdf=True
                    )

        return tracked_files

    def _run_nodeeditor_alltimes(self, outdir, timestamp, tracked_files, detect_files, grid_resol='native'):
        """
        Run the Tempest node file editor on all of the tracked files.

        :param list tracked_files: The paths (as strings) of the tracked files
            produced by the stitching process. The files are in order of
            track type.
        :param str timestamp: The timestep of data/tracking to process
        :returns: The path to the detect file, the path to the tracked file,
            the path of the sea level pressure input netCDF file and , all as
            strings.
        :rtype: tuple
        """

        nodeedit_files = []
        self.logger.info(f"start NodeFileEditor")

        # for each timestamp in the detect files, we want to find the 
        # input data for the nodeedit and save it
        for track_type in self.track_types:
            tracking_phase_commands = self._construct_command(track_type)["nodeedit"]
            self.logger.info(f"NodeFileEditor {track_type} {tracking_phase_commands}")
            if tracking_phase_commands is None:
                continue

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detect", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            # find out format of output file and name track file accordingly
            fileformat = self._command_parameter(
                track_type, 
                "nodeedit", 
                "out_nodefile_format")
            if fileformat is None or fileformat == "":
                fileformat = "gfdl"
            self.logger.debug(f"File format of stitch file {fileformat}")

            detects = detect_files[track_type]

            processed_files_needed = []
            for det in detects:
                timestamp_start = os.path.basename(det).split('_')[2].split('-')[0]
                timestamp_end = os.path.basename(det).split('_')[2].split('-')[1]

                processed_files, variable_units = \
                        self._identify_processed_files(
                            timestamp_start,
                            timestamp_end,
                            grid_resol=grid_resol,
                            variables=self.nodeedit_vars)
                self.logger.info(f"processed_files {processed_files}")

                for var in self.nodeedit_vars:
                    processed_files_needed.append(processed_files[var])

            self.logger.info(f"processed_files_needed {processed_files_needed}")

            tracked_file = tracked_files[track_type]
            start_date = os.path.basename(tracked_file).split("_")[2].split("-")[0]
            end_date = os.path.basename(tracked_file).split("_")[2].split("-")[1]
            nodeedit_file_alltimes = self.file_format.format(
                outdir=outdir, 
                app="tracknodeedit",
                runid=self.runid, 
                period="", 
                timestamp=start_date, 
                timestamp_next=end_date, 
                data_freq=timefilter, 
                track_type=track_type)[:-3]+fileformat.strip('\"')

            cmd_io = '{} --in_nodefile {} --out_nodefile {} '.format(
                    self.tc_editor_script,
                    tracked_file,
                    nodeedit_file_alltimes)

            in_file_list = os.path.join(self.outdir, 'in_file_list_edit.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(processed_files_needed)
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)
            cmd_io += '--in_data_list '+in_file_list+' '

            cmd_edit = cmd_io + tracking_phase_commands
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
                    f"EXCEPTION found in TempestExtreme editor output\n"
                    f"{sts.stdout}"
                )
                raise RuntimeError(msg)
            nodeedit_files.append(nodeedit_file_alltimes)

            # only convert to netcdf if this is csv
            if "csv" in fileformat:
                self._read_write_and_plot_tracks(
                    nodeedit_file_alltimes,
                    track_type,
                    processed_files_needed[-1],
                    start_date,
                    end_date,
                    self.variable_units,
                    title_prefix=f"{self.runid} {self.resolution_code} "
                    f"{start_date}_{end_date}",
                    title_suffix=f"{track_type} tracks",
                    write_to_netcdf=True
                )

        return nodeedit_files

    def _run_variable_processor(self, outdir, timestamp, timestamp_next,
                                grid_resol='native', proc='varproc1'):
        """
        Run the Tempest variable processor on all of the tracked files.

        :param str timestamp: The timestep of data/tracking to process
        :param str timestamp_next: The next timestep of data/tracking to process
        :returns: The path to the detect file, the path to the tracked file,
            the path of the sea level pressure input netCDF file and , all as
            strings.
        :rtype: tuple
        """

        varproc_files = []

        processed_files_curr, variable_units_curr = \
            self._identify_processed_files(timestamp, timestamp_next,
                                           grid_resol=grid_resol)
        self.logger.info(f"processed_files_curr {processed_files_curr}")

        if proc == "varproc1":
            varproc_vars = self.varproc1_vars
        elif proc == "varproc2":
            varproc_vars = self.varproc2_vars

        for index, track_type in enumerate(self.track_types):
            varproc_commands = self._construct_command(track_type)[proc]
            if varproc_commands is None:
                continue

            processed_files_in = []
            # put the path names of the files with variables needed into list
            for var in varproc_vars[:-1]:
                processed_files_in.append(processed_files_curr[var])

            self.logger.info(f"processed_files_needed {processed_files_in}")

            in_file_list = os.path.join(outdir, 'in_file_list_'+proc+'.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(processed_files_in)
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)

            out_file_list = os.path.join(outdir, 'out_file_list_'+proc+'.txt')
            varout_path = processed_files_in[0].replace(varproc_vars[0], varproc_vars[-1])
            with open(out_file_list, 'w') as fh:
                text_str = varout_path
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)

            cmd_io = '{} --in_data_list {} --out_data_list {} '.format(
                    self.tc_varproc_script,
                    in_file_list,
                    out_file_list)

            cmd_varproc = cmd_io + varproc_commands
            self.cmd_varproc_type[track_type] = cmd_varproc
            self.logger.info(f"Varproc command {cmd_varproc}")

            sts = subprocess.run(
                cmd_varproc,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
            self.logger.debug(sts.stdout)
            if "EXCEPTION" in sts.stdout:
                msg = (
                    f"EXCEPTION found in TempestExtreme varproc output\n"
                    f"{sts.stdout}"
                )
                raise RuntimeError(msg)
            varproc_files.append(varout_path)
            self.variables_rename.append(varproc_vars[-1])

        return varproc_files

    def _run_detectblobs(self, outdir, timestamp, timestamp_next, 
                         detect_file,
                         grid_resol='native'):
        """
        Run the Tempest DetectBlobs on the vor_cyc, ua_925 and va_925.

        :param str timestamp: The timestep of data/tracking to process
        :param str timestamp_next: The next timestep of data/tracking to process
        :returns: The path to the detect file, the path to the tracked file,
            the path of the sea level pressure input netCDF file and , all as
            strings.
        :rtype: tuple
        """

        detectblobs_files = []

        for index, track_type in enumerate(self.track_types):
            detectblobs_commands = self._construct_command(track_type)["detectblobs"]
            if detectblobs_commands is None:
                continue

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detectblobs", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            processed_files_curr, variable_units_curr = \
                    self._identify_processed_files(
                        timestamp,
                        timestamp_next,
                        grid_resol=grid_resol,
                        variables=self.detectblobs_vars)
            self.logger.info(f"processed_files_curr {processed_files_curr}")
            processed_files_in = []
            for var in self.detectblobs_vars:
                processed_files_in.append(processed_files_curr[var])

            self.logger.info(f"processed_files_needed {processed_files_in}")

            in_file_list = os.path.join(self.outdir, 'in_file_list_detectblobs.txt')
            with open(in_file_list, 'w') as fh:
                text_str = ';'.join(processed_files_in)
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)

            for det in detect_file:
                if track_type in det:
                    detectblobs_file = self.file_format_nc.format(
                        outdir=outdir, 
                        runid=self.runid, 
                        app="detectblobs", 
                        period="", 
                        timestamp=timestamp, 
                        timestamp_next=timestamp_next, 
                        data_freq=timefilter, 
                        track_type=track_type)

            out_list = os.path.join(self.outdir, 'out_list_detectblobs.txt')
            with open(out_list, 'w') as fh:
                text_str = detectblobs_file
                self.logger.debug(f"file_list {text_str}")
                fh.write(text_str)

            cmd_io = '{} --in_data_list {} --out_list {} '.format(
                    self.tc_detectblobs_script,
                    in_file_list,
                    out_list)

            cmd_detectblobs = cmd_io + detectblobs_commands
            self.cmd_detectblobs_type[track_type] = cmd_detectblobs
            self.logger.info(f"Varproc command {cmd_detectblobs}")

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

            detectblobs_files.append(detectblobs_file)

        return detectblobs_files

    def _run_blobstats(self, outdir, timestamp, detectblobs_files,
                                grid_resol='native'):
        """
        Run the Tempest BlobStats on the detectblobs file.

        :param str timestamp: The timestep of data/tracking to process
        :returns: The path to the detect file, the path to the tracked file,
            the path of the sea level pressure input netCDF file and , all as
            strings.
        :rtype: tuple
        """

        blobstats_files = []
        for index, track_type in enumerate(self.track_types):
            blobstats_commands = self._construct_command(track_type)["blobstats"]
            if blobstats_commands is None:
                continue

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detectblobs", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            detectblobs = detectblobs_files[track_type]
            self.logger.debug(f"detectblobs_files {detectblobs_files}")

            if len(detectblobs) >= 1:
                start_date = os.path.basename(detectblobs[0]).split("_")[2].split('-')[0]
                end_date = os.path.basename(detectblobs[-1]).split("_")[2].split('-')[1]
                end_month = end_date[4:6]
                if end_month == "12":
                    self.logger.debug(f"end of year {end_date} {end_month}")

                detectblobs_alltimes = self.file_format_nc.format(
                    outdir=outdir, 
                    runid=self.runid, 
                    app="detectblobs", 
                    period="period_", 
                    timestamp=start_date, 
                    timestamp_next=end_date, 
                    data_freq=timefilter, 
                    track_type=track_type)

                cmd = "ncrcat -O "+' '.join(detectblobs)+" "+detectblobs_alltimes
                self.logger.debug(f"ncrcat detectblobs {cmd}")
                sts = self._run_cmd(cmd, check=False)

                blobstats_alltimes = self.file_format.format(
                    outdir=outdir, 
                    runid=self.runid, 
                    app="blobstats", 
                    period="", 
                    timestamp=start_date, 
                    timestamp_next=end_date, 
                    data_freq=timefilter, 
                    track_type=track_type)

            cmd_io = '{} --in_file {} --out_file {} '.format(
                    self.tc_blobstats_script,
                    detectblobs_alltimes,
                    blobstats_alltimes)

            cmd_blobstats = cmd_io + blobstats_commands
            self.cmd_blobstats_type[track_type] = cmd_blobstats
            self.logger.info(f"Blobstats command {cmd_blobstats}")

            sts = subprocess.run(
                cmd_blobstats,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
            self.logger.debug(sts.stdout)
            if "EXCEPTION" in sts.stdout:
                msg = (
                    f"EXCEPTION found in TempestExtreme blobstats output\n"
                    f"{sts.stdout}"
                )
                raise RuntimeError(msg)

            blobstats_files.append(blobstats_alltimes)

        return blobstats_files

    def _stitchblobs_file(self, detect_file, tracked_file, track_type):
        """
        Run the Tempest stitching on the specified file.

        :param str detect_file: The path of the detect input file.
        :param str tracked_file: The path of the stitched output file.
        :param str track_type: The name of the type of tracking required.
        """

        stitchblobs_commands = self._construct_command(track_type)["stitchblobs"]
        if stitchblobs_commands is None:
            return

        self.logger.debug(f"tracked_file {tracked_file}")

        # stitch detects together
        cmd_stitchblobs_io = "{} --in {} --out {} ".format(
            self.tc_stitchblobs_script, detect_file, tracked_file
        )
        cmd_stitchblobs = cmd_stitchblobs_io + stitchblobs_commands
        self.cmd_stitchblobs_type[track_type] = cmd_stitchblobs
        self.logger.info(f"Stitch command {cmd_stitchblobs}")

        sts = subprocess.run(
            cmd_stitchblobs,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=True,
        )
        self.logger.debug(f"sts err {sts.stdout}")

        if "EXCEPTION" in sts.stdout:
            msg = (
                f"EXCEPTION found in TempestExtreme stitchblobs output\n" f"{sts.stdout}"
            )
            raise RuntimeError(msg)

    def _run_stitchblobs_alltimes(self, outdir, timestamp, detectblobs_files):
        """
        Concatenate the detectblobs files together and then
        stitch this file.

        :param str outdir: directory for the files
        :param str timestamp: The current timestep of data/tracking to process
        :param str timestamp: The previous timestep of data/tracking to process
        :returns: The string paths to the annual stitched file in a list.
        :rtype: list
        """

        trackedblobs_files = {}
        year = timestamp[:4]
        yearm1 = str(int(year)-1)
        month = timestamp[4:6]
        for track_type in self.track_types:
            tracking_phase_commands = self._construct_command(track_type)["stitchblobs"]
            if tracking_phase_commands is None:
                continue

            # find out time frequency of the data used and name the detect 
            # file accordingly
            timefilter = self._command_parameter(
                track_type, 
                "detectblobs", 
                "timefilter")
            if timefilter is not None and len(timefilter) > 0:
                timefilter = timefilter.strip('\"')[:2]
            else:
                timefilter = self.data_frequency

            detectblobs = detectblobs_files[track_type]
            self.logger.debug(f"detectblobs_files {detectblobs_files}")

            if len(detectblobs) >= 1:
                start_date = os.path.basename(detectblobs[0]).split("_")[2].split('-')[0]
                end_date = os.path.basename(detectblobs[-1]).split("_")[2].split('-')[1]
                start_year = start_date[:4]
                end_year = end_date[:4]
                end_month = end_date[4:6]
                if end_month == "01":
                    self.logger.debug(f"end of year {end_date} {end_month}")

                detectblobs_alltimes = self.file_format_nc.format(
                    outdir=outdir, 
                    runid=self.runid, 
                    app="detectblobs",
                    period="period_",
                    timestamp=start_date, 
                    timestamp_next=end_date, 
                    data_freq=timefilter, 
                    track_type=track_type)

                cmd = "ncrcat -O "+' '.join(detectblobs)+" "+detectblobs_alltimes
                self.logger.debug(f"ncrcat detectblobs {cmd}")
                sts = self._run_cmd(cmd, check=False)

                trackblobs_alltimes = self.file_format_nc.format(
                    outdir=outdir, 
                    runid=self.runid, 
                    app="trackblobs",
                    period="",
                    timestamp=start_date, 
                    timestamp_next=end_date, 
                    data_freq=timefilter, 
                    track_type=track_type)

                self._stitchblobs_file(detectblobs_alltimes, trackblobs_alltimes, track_type)
                #if self.track_by_year:
                #    self._remove_overlapping_tracks(trackblobs_alltimes, year)

                trackedblobs_files[track_type] = trackblobs_alltimes

        return trackedblobs_files

    def _read_write_and_plot_tracks(
        self,
        tracked_file,
        track_type,
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
        :param dict variable_units: The units of each output variable, as derived from
        :                           the input data files
        :param str title_prefix: The title for the plot.
        :param str title_suffix: The title for the plot.
        :param bool include_num_tracks: Include the number of tracks within the plot
        :                               title
        :param bool write_to_netcdf: Whether to write tracks to netcdf file or not
        """

        # storms returned as dictionary with keys: length, lon, lat, year, month, day,
        # hour, step
        track_step = os.path.basename(tracked_file).split('_')[1]
        if track_step == 'track':
            step = 'stitch'
        elif track_step == 'tracknodeedit':
            step = 'nodeedit'
        #hr_frequency = int(self.data_frequency.strip('h'))
        timefreq = os.path.basename(tracked_file).split('_')[3]
        hr_frequency = int(timefreq.strip("h"))

        filetype = os.path.basename(tracked_file).split('.')[-1]
        if "csv" in filetype:
            storms, columns = get_trajectories_csv(tracked_file, nc_file_in, hr_frequency,
                                  self.column_names[track_type+'_'+step])
        else:
            storms = get_trajectories(tracked_file, nc_file_in, hr_frequency,
                                  self.column_names[track_type+'_'+step])

        self.logger.debug(f"got storms {len(storms)} from {tracked_file}")

        if len(storms) == 0:
            return
        # self.logger.debug(f"got storms {storms[0]}")

        self.logger.debug(f"netcdf file in {nc_file_in}")
        nc_file_out = os.path.join(
            os.path.dirname(tracked_file),
            os.path.basename(tracked_file).split(".")[0]+"_nogaps.nc")

        if "csv" in filetype:
            filled_csv = nc_file_out[:-3]+".csv"
            with open(filled_csv, "w", newline="") as f:
                for istorm, storm in enumerate(storms):
                    w = csv.writer(f, delimiter=',')
                    key_list = list(columns)
                    limit = len(storm[columns[0]])
                    if istorm == 0:
                        w.writerow(columns)
                    for  i in range(limit):
                        w.writerow([storm[x][i] for x in key_list])

        # write the storms to a netcdf file
        if write_to_netcdf:
            with Dataset(nc_file_in, "r") as x:
                time = x.variables["time"]
                calendar = time.calendar
                calendar_units = self.calendar_units

            self.logger.debug(f"open netcdf file {nc_file_out}")
            save_trajectories_netcdf(
                self.outdir,
                nc_file_out,
                storms,
                calendar,
                calendar_units,
                variable_units,
                timefreq,
                self.suiteid,
                self.resolution_code,
                self.cmd_detect_type[track_type],
                self.cmd_stitch_type[track_type],
                self.column_names[track_type+'_'+step],
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

        super()._get_app_options()

        self.tc_detect_script = self.app_config.get_property(
            "common", "tc_detect_script"
        )
        self.tc_stitch_script = self.app_config.get_property(
            "common", "tc_stitch_script"
        )
        self.tc_editor_script = self.app_config.get_property(
            "common", "tc_editor_script"
        )
        self.tc_varproc_script = self.app_config.get_property(
            "common", "tc_varproc_script"
        )
        self.tc_detectblobs_script = self.app_config.get_property(
            "common", "tc_detectblobs_script"
        )
        self.tc_blobstats_script = self.app_config.get_property(
            "common", "tc_blobstats_script"
        )
        self.tc_stitchblobs_script = self.app_config.get_property(
            "common", "tc_stitchblobs_script"
        )
        self.nodeedit_vars = eval(self.app_config.get_property("common",
                                                               "nodeedit_vars"))
        self.varproc1_vars = eval(self.app_config.get_property("common",
                                                               "varproc1_vars"))
        self.varproc2_vars = eval(self.app_config.get_property("common",
                                                               "varproc2_vars"))
        self.detectblobs_vars = eval(self.app_config.get_property("common",
                                                               "detectblobs_vars"))
        # track_types is a Python list so eval converts str to list
        self.track_types = eval(
            self.app_config.get_property("common", "track_types"))
        self.variables_rename = eval(self.app_config.get_property("common",
                                                            "tc_variables"))
