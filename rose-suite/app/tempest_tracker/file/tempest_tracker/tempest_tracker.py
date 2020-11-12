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
        metadata = self._collect_metadata(True)

        # Run the tracking
        (candidatefile, trackedfile, nc_file,
         tracking_period, frequency) = self._run_tracking(metadata)

        # Run the plotting (if required)
        if os.stat(candidatefile).st_size > 0:
            if os.stat(trackedfile).st_size > 0:
                title = self.track_type + ' tracks'
                if self.plot_tracks:
                    self._read_and_plot_tracks(
                        trackedfile, nc_file,
                        tracking_period,
                        frequency,
                        title=title
                    )
            else:
                self.logger.error(f'candidatefile has data but no tracks '
                                  f'{candidatefile}')
        else:
            self.logger.error(f'candidatefile is empty {candidatefile}')

    def _run_tracking(self, metadata):
        """
        Run the tracking.

        :param dict metadata: The collection of metadata.
        :returns: The path to the candidate file, the path to the tracked file,
            the path of the sea level pressure input netCDF file and the period
            between time points in the input file, all as strings and the
            frequency of time points in the file in hours as an integer.
        :rtype: tuple
        """
        self.logger.error(f'cwd {os.getcwd()}')

        filenames, tracking_period, frequency = self._generate_data_files()

        candidatefile = os.path.join(
            self.output_directory,
            f'candidate_file_{self.cylc_task_cycle_time}_{self.track_type}.txt'
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
                candidatefile)
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
                candidatefile)
            )

        tracking_phase_commands = self._construct_command()
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
            f'track_file_{tracking_period}_{self.track_type}.txt'
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

        return (candidatefile, trackedfile, filenames['pslfile'],
                tracking_period, frequency)

    def _collect_metadata(self, zg_available=False):
        """
        Collect metadata into a single dictionary.

        :param bool zg_available: True if the zg (geopotential height) variable
            is available to use.
        :returns: The collected metadata in a single dictionary
        :rtype: dict
        """
        metadata = {}
        metadata['model'] = self.um_runid
        metadata['resol'] = self.resolution_code
        metadata['model_name'] = metadata['model']
        metadata['u_var'] = 'x_wind'
        metadata['v_var'] = 'y_wind'
        if not zg_available:
            metadata['t_var'] = 'ta'
        else:
            metadata['t_var'] = 'unknown'
        metadata['lat_var'] = 'lat'
        metadata['lon_var'] = 'lon'
        metadata['topo_var'] = 'surface_altitude'
        return metadata

    def _construct_command(self):
        """
        Read the TempestExtreme command line parameters from the configuration.

        :returns: A dictionary with keys of `detect` and `stitch` and the
            values for each of these is a string containing the command
            line parameters for each of these TempestExtreme steps. The
            parameters are sorted into alphabetical order in each line.
        :rtype: dict
        """

        commands = {}
        for step in ['detect', 'stitch']:
            step_config = self.app_config.section_to_dict(
                f'{self.track_type}_{step}'
            )

            step_arguments = [f'--{parameter} {step_config[parameter]}'
                              for parameter in sorted(list(step_config.keys()))
                              if step_config[parameter]]
            commands[step] = ' '.join(step_arguments)
        return commands

    def _generate_data_files(self):
        """
        Identify and then fix the grids and var_names in the input files.

        :returns: A tuple of a dictionary of the files found for this period,
            a string containing the  start and end timestrings that the
            tracking was run over and the integer frequency in hours between
            time points.
        :rtype: tuple
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

        periods = []
        frequencies = []
        for filename in files:
            period = os.path.basename(filename).split('_')[3]
            frequency = os.path.basename(filename).split('_')[2]
            periods.append(period)
            frequencies.append(frequency)
            unique_periods = list(set(periods))
            unique_frequencies = list(set(frequencies))
        if len(unique_periods) == 0:
            msg = 'No tracked_file periods found'
            self.logger.error(msg)
            raise RuntimeError(msg)
        elif len(unique_periods) != 1:
            msg = 'No tracked_file periods found'
            self.logger.error(msg)
            raise RuntimeError(msg)
        else:
            period = unique_periods[0]
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
                frequency_value = int(components[1])

        filenames = self._process_input_files(period)

        return filenames, period, frequency_value

    def _process_input_files(self, period):
        """
        Identify and then fix the grids and var_names in the input files.

        :param str period:
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

        # TODO add names of files generated (and required) to documentation
        filename_format = 'atmos_{}a_6h_{}_pt-{}.nc'

        for filetype in filetypes_required:
            filename = filename_format.format(
                self.um_runid,
                period,
                variables_required[filetype]['fname']
            )

            input_path = os.path.join(self.input_directory, filename)
            if not os.path.exists(input_path):
                msg = f'Unable to find expected input file {input_path}'
                self.logger.error(msg)
                raise RuntimeError(msg)

            output_path = os.path.join(self.output_directory, filename)

            reference_name = filename_format.format(
                self.um_runid,
                period,
                variables_required['pslfile']['fname']
            )
            reference_path = os.path.join(self.input_directory, reference_name)

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

            # newer tempest can specify names of latitude, longitude
            # fix_lat_lon(fname_local)
            processed_filenames[filetype] = output_path

        return processed_filenames

    def _read_and_plot_tracks(self, trackedfile, nc_file, tracking_period,
                              frequency, title=''):
        """
        Read and then plot the tracks

        :param str trackedfile: The path to the track file.
        :param str nc_file: The path to an input data netCDF file, which is
            used to gatehr additional information about the dates and calendars
            used in the data.
        :param str tracking_period: The model output period.
        :param int frequency: The time in hours between data points.
        :param str title: The title for the plot.
        """
        storms = get_trajectories(trackedfile, nc_file, frequency)
        num_trajectories = count_trajectories(storms)

        title_suffix = title if title else "TempestExtremes TCs"
        title_full = (f"{self.um_runid} {self.resolution_code} "
                      f"{tracking_period} {num_trajectories} {title_suffix}")
        filename = trackedfile[:-4]+'.png'
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
        self.extended_files = self.app_config.get_bool_property('common',
                                                                'extended_files')
        self.plot_tracks = self.app_config.get_bool_property('common',
                                                             'plot_tracks')
        self.track_type = self.app_config.get_property('common', 'track_type')

    def _get_environment_variables(self):
        """
        Get the required environment variables from the suite. A list and
        explanation of the required environment variables is included in the
        documentation.
        """
        self.um_runid = os.environ['RUNID']
        self.cylc_task_cycle_time = os.environ['CYLC_TASK_CYCLE_TIME']

