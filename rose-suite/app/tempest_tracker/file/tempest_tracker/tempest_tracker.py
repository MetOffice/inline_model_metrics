# (C) British Crown Copyright 2020, Met Office.
# Please see LICENSE for license details.
import glob
import os
import re
import shutil
import subprocess
import sys

import iris

from afterburner.apps import AbstractApp

from tempest_helper import (count_trajectories, get_trajectories,
                            plot_trajectories_cartopy)


class TempestTracker(AbstractApp):
    """
    Run the TempestExtreme tracker inline with a climate model
    """
    def __init__(self, arglist=None, **kwargs):
        super(TempestTracker, self).__init__(version='1.0', **kwargs)
        self._parse_args(arglist, desc="Inline TempestExtreme tracking")
        self._parse_app_config()
        self._set_message_level()

        # Instance attributes set later
        self.time_range = None
        self.frequency = None

    @property
    def cli_spec(self):
        """
        Defines the command-line interface specification for the app.
        """
        return [
            {'names': ['-c', '--config-file'],
             'help': 'Pathname of app configuration file'},
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
                msg = (f'Unable to create output directory '
                       f'{self.output_directory}')
                self.logger.error(msg)
                sys.exit(1)

        self.logger.debug(f'CYLC_TASK_CYCLE_TIME {self.cylc_task_cycle_time}, '
                          f'um_runid {self.um_runid}')

        # Run the tracking
        candidate_files, tracked_files, nc_file = self._run_tracking()

        # Run the plotting (if required and data available)
        for index, track_type in enumerate(self.track_types):
            candidate_file = candidate_files[index]
            tracked_file = tracked_files[index]
            if os.stat(candidate_file).st_size > 0:
                if os.stat(tracked_file).st_size > 0:
                    title = track_type + ' tracks'
                    if self.plot_tracks:
                        self._read_and_plot_tracks(
                            tracked_file, nc_file, title=title
                        )
                else:
                    self.logger.error(f'candidate file has data but there are '
                                      f'no tracks in {tracked_file}')
            else:
                self.logger.error(f'candidate file is empty '
                                  f'{candidate_file}')

    def _run_tracking(self):
        """
        Run the tracking.

        :returns: The path to the candidate file, the path to the tracked file,
            the path of the sea level pressure input netCDF file and , all as
            strings.
        :rtype: tuple
        """
        self.logger.error(f'cwd {os.getcwd()}')

        filenames = self._generate_data_files()

        candidate_files = []
        tracked_files = []

        for track_type in self.track_types:
            self.logger.debug(f'Runing {track_type} tracking')
            candidatefile = os.path.join(
                self.output_directory,
                f'candidate_file_{self.cylc_task_cycle_time}_{track_type}.txt'
            )
            self.logger.debug(f'candidatefile {candidatefile}')

            if not self.extended_files:
                # identify candidates
                cmd_io = '{} --in_data "{};{};{};{};{}" --out {} '.format(
                    self.tc_detect_script,
                    filenames['pslfile'],
                    filenames['zgfile'],
                    filenames['ufile'],
                    filenames['vfile'],
                    filenames['topofile'],
                    candidatefile
                )
            else:
                cmd_io = ('{} --in_data "{};{};{};{};{};{};{};{};{};{};{};{};{}" '
                          '--out {} '.format(
                    self.tc_detect_script,
                    filenames['pslfile'],
                    filenames['zgfile'],
                    filenames['ufile'],
                    filenames['vfile'],
                    filenames['rvfile'],
                    filenames['u10mfile'],
                    filenames['v10mfile'],
                    filenames['ws10mfile'],
                    filenames['viwvefile'],
                    filenames['viwvnfile'],
                    filenames['tafile'],
                    filenames['rvT63file'],
                    filenames['topofile'],
                    candidatefile
                )
                )

            tracking_phase_commands = self._construct_command(track_type)
            cmd_detect = cmd_io + tracking_phase_commands['detect']
            self.logger.info(f'Detect command {cmd_detect}')

            sts = subprocess.run(cmd_detect, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, universal_newlines=True,
                                 check=True)
            self.logger.debug(sts.stdout)
            if 'EXCEPTION' in sts.stdout:
                msg = (f'EXCEPTION found in TempestExtreme detect output\n'
                       f'{sts.stdout}')
                raise RuntimeError(msg)

            trackedfile = os.path.join(
                self.output_directory,
                f'track_file_{self.time_range}_{track_type}.txt'
            )

            # stitch candidates together
            cmd_stitch_io = '{} --in {} --out {} '.format(
                self.tc_stitch_script,
                candidatefile,
                trackedfile
            )
            cmd_stitch = cmd_stitch_io + tracking_phase_commands['stitch']
            self.logger.info(f'Stitch command {cmd_stitch}')
            sts = subprocess.run(cmd_stitch, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, universal_newlines=True,
                                 check=True)
            self.logger.debug(f'sts err {sts.stdout}')
            candidate_files.append((candidatefile))
            tracked_files.append(trackedfile)

        return candidate_files, tracked_files, filenames['pslfile']

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
        for step in ['detect', 'stitch']:
            step_config = self.app_config.section_to_dict(
                f'{track_type}_{step}'
            )

            step_arguments = [f'--{parameter} {step_config[parameter]}'
                              for parameter in sorted(list(step_config.keys()))
                              if step_config[parameter]]
            commands[step] = ' '.join(step_arguments)
        return commands

    def _generate_data_files(self):
        """
        Identify and then fix the grids and var_names in the input files.

        The time_range and frequency attributes are set when this method runs.

        :returns: A dictionary of the files found for this period and a string
            containing the period between samples in the input data.
        :rtype: dict
        """
        timestamp_day = self.cylc_task_cycle_time[:8]
        self.logger.debug(f'time_stamp_day {timestamp_day}')
        file_search = os.path.join(
            self.input_directory,
            f'atmos_{self.um_runid}*{timestamp_day}??-*-slp.nc')
        self.logger.debug(f'file_search {file_search}')
        files = sorted(glob.glob(file_search))
        self.logger.debug(f'files {files}')
        if not files:
            msg = f'No input files found for glob pattern {file_search}'
            self.logger.error(msg)
            raise RuntimeError(msg)

        time_ranges = []
        unique_ranges = []
        frequencies = []
        unique_frequencies = []
        for filename in files:
            time_range = os.path.basename(filename).split('_')[3]
            frequency = os.path.basename(filename).split('_')[2]
            time_ranges.append(time_range)
            frequencies.append(frequency)
            unique_ranges = list(set(time_ranges))
            unique_frequencies = list(set(frequencies))
        if len(unique_ranges) == 0:
            msg = 'No tracked_file periods found'
            self.logger.error(msg)
            raise RuntimeError(msg)
        elif len(unique_ranges) != 1:
            msg = 'No tracked_file periods found'
            self.logger.error(msg)
            raise RuntimeError(msg)
        else:
            self.time_range = unique_ranges[0]
        if len(unique_frequencies) == 0:
            msg = 'No tracked_file frequencies found'
            self.logger.error(msg)
            raise RuntimeError(msg)
        elif len(unique_frequencies) != 1:
            msg = 'No tracked_file frequencies found'
            self.logger.error(msg)
            raise RuntimeError(msg)
        else:
            frequency = unique_frequencies[0]
            components = re.match(r'(\d+)', frequency)
            if not components:
                msg = r'No digit found in frequency {frequency}'
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
        filetypes_required = ['pslfile', 'zgfile', 'ufile', 'vfile',
                              'ws10mfile', 'rvfile', 'rvT63file', 'u10mfile',
                              'u10mfile', 'v10mfile', 'viwvefile', 'viwvnfile',
                              'tafile']
        processed_filenames = {}

        processed_filenames['topofile'] = self.orography_file

        variables_required = {}
        variables_required['pslfile'] = {'fname': 'slp'}
        variables_required['zgfile'] = {'fname': 'zg', 'varname_new': 'zg'}
        variables_required['ufile'] = {'fname': 'ua'}
        variables_required['vfile'] = {'fname': 'va'}
        variables_required['ws10mfile'] = {'fname': 'wsas'}
        variables_required['rvfile'] = {'fname': 'rv', 'varname_new': 'rv'}
        variables_required['rvT63file'] = {'fname': 'rvT63',
                                           'varname_new': 'rvT63'}
        variables_required['u10mfile'] = {'fname': 'u10m', 'varname_new': 'uas'}
        variables_required['v10mfile'] = {'fname': 'v10m', 'varname_new': 'vas'}
        variables_required['viwvefile'] = {'fname': 'viwve',
                                           'varname_new': 'viwve'}
        variables_required['viwvnfile'] = {'fname': 'viwvn',
                                           'varname_new': 'viwvn'}
        variables_required['tafile'] = {'fname': 'ta',
                                        'varname_new': 'ta'}

        filename_format = 'atmos_{}a_{}h_{}_pt-{}.nc'

        reference_name = filename_format.format(
            self.um_runid,
            self.frequency,
            self.time_range,
            variables_required['pslfile']['fname']
        )
        reference_path = os.path.join(self.input_directory, reference_name)

        for filetype in filetypes_required:
            filename = filename_format.format(
                self.um_runid,
                self.frequency,
                self.time_range,
                variables_required[filetype]['fname']
            )

            input_path = os.path.join(self.input_directory, filename)
            if not os.path.exists(input_path):
                msg = f'Unable to find expected input file {input_path}'
                self.logger.error(msg)
                raise RuntimeError(msg)

            output_path = os.path.join(self.output_directory, filename)

            if (variables_required[filetype]['fname'] in
                    ['ua', 'va', 'wsas', 'rv', 'rvT63', 'u10m', 'v10m']):
                # regrid u and v to t grid and rename variable if necessary
                cube = iris.load_cube(input_path)
                reference = iris.load_cube(reference_path)
                regridded = cube.regrid(reference, iris.analysis.Linear())
                if 'varname_new' in variables_required[filetype]:
                    regridded.var_name = (variables_required[filetype]
                                                            ['varname_new'])
                iris.save(regridded, output_path)
            elif (variables_required[filetype]['fname'] in
                  ['ta', 'zg']):
                # rename variables only - could do with ncrename instead
                cube = iris.load_cube(input_path)
                if 'varname_new' in variables_required[filetype]:
                    cube.var_name = variables_required[filetype]['varname_new']
                iris.save(cube, output_path)
            else:
                if not os.path.exists(output_path):
                    shutil.copyfile(input_path, output_path)

            processed_filenames[filetype] = output_path

        return processed_filenames

    def _read_and_plot_tracks(self, tracked_file, nc_file, title=''):
        """
        Read and then plot the tracks

        :param str tracked_file: The path to the track file.
        :param str nc_file: The path to an input data netCDF file, which is
            used to gather additional information about the dates and calendars
            used in the data.
        :param str title: The title for the plot.
        """
        storms = get_trajectories(tracked_file, nc_file, self.frequency)
        num_trajectories = count_trajectories(storms)

        title_suffix = title if title else "TempestExtremes TCs"
        title_full = (f"{self.um_runid} {self.resolution_code} "
                      f"{self.time_range} {num_trajectories} {title_suffix}")
        filename = tracked_file[:-4]+'.png'
        plot_trajectories_cartopy(storms, filename, title=title_full)

    def _get_app_options(self):
        """Get commonly used configuration items from the config file"""
        self.resolution_code = self.app_config.get_property('common',
                                                            'resolution')

        self.input_directory = self.app_config.get_property('common',
                                                            'input_directory')
        self.output_directory = self.app_config.get_property('common',
                                                             'output_directory')
        self.tc_detect_script = self.app_config.get_property('common',
                                                             'tc_detect_script')
        self.tc_stitch_script = self.app_config.get_property('common',
                                                             'tc_stitch_script')
        self.psl_std_name = self.app_config.get_property('common',
                                                         'psl_std_name')
        self.orography_file = self.app_config.get_property('common',
                                                           'orography_file')
        self.extended_files = self.app_config.get_bool_property(
            'common',
            'extended_files'
        )
        self.plot_tracks = self.app_config.get_bool_property('common',
                                                             'plot_tracks')
        # track_types is a Python list so eval converts str to list
        self.track_types = eval(self.app_config.get_property('common',
                                                             'track_types'))

    def _get_environment_variables(self):
        """
        Get the required environment variables from the suite. A list and
        explanation of the required environment variables is included in the
        documentation.
        """
        self.um_runid = os.environ['RUNID']
        self.cylc_task_cycle_time = os.environ['CYLC_TASK_CYCLE_TIME']
