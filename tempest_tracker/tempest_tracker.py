# (C) British Crown Copyright 2020, Met Office.
# Please see LICENSE for license details.
import os
import subprocess
import sys

from afterburner.apps import AbstractApp


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
        metadata = self.collect_metadata(True)

        # TODO please define
        nVars_stitch = 6 # lon, lat, psl, topo

        # TODO Move to config
        track_type = 'tc_slp'

        candidatefile, trackedfile, nc_file, tracking_period = self.run_tracking(
            metadata,
            tracking_type=track_type
        )

        if os.stat(candidatefile).st_size > 0:
            if os.stat(trackedfile).st_size > 0:
                title = track_type + ' tracks'
                self.read_and_plot_tracks(
                    trackedfile, nc_file, tracking_period, self.um_runid,
                    self.resolution_code, title=title, feature_type=track_type,
                    nVars_stitch=nVars_stitch
                )
            else:
                self.logger.error(f'candidatefile has data but no tracks '
                                  f'{candidatefile}')
        else:
            self.logger.error(f'candidatefile is empty {candidatefile}')

    def run_tracking(self, metadata, tracking_type='tc_slp'):
        """
        Run the tracking.

        :param dict metadata: The collection of metadata
        :param str tracking_type: The type of tracking to run
        """
        print('cwd ', os.getcwd())
        # nl_file = './namelist_'+tracking_type.lower()+'.nl'
        nl_file = os.path.join(os.getcwd(), tracking_type + '.nl')
        cmds = self.construct_command(nl_file, tracking_type)

        filenames, period, output_uv, psl_var = self.data_files(self.cylc_task_cylc_time, self.um_runid,
                                                           self.resolution_code, self.input_directory,
                                                           self.output_directory, metadata)

        candidatefile = os.path.join(self.output_directory,
                                     'candidate_file_' + self.cylc_task_cylc_time + '_' + tracking_type + '.txt')

        files_extend = False
        if not files_extend:
            # identify candidates
            cmd_io = '{} --in_data "{};{};{};{};{}" --out {} '.format(
                self.tc_detect_script, filenames['pslfile'], filenames['zgfile'],
                filenames['ufile'], filenames['vfile'], filenames['topofile'],
                candidatefile)
        else:

            cmd_io = '{} --in_data "{};{};{};{};{};{};{};{};{};{};{};{};{}" --out {} '.format(
                self.tc_detect_script, filenames['pslfile'], filenames['zgfile'],
                filenames['ufile'], filenames['vfile'], filenames['rvfile'],
                filenames['u10mfile'], filenames['v10mfile'],
                filenames['ws10mfile'],
                filenames['viwvefile'], filenames['viwvnfile'],
                filenames['tafile'],
                filenames['rvT63file'], filenames['topofile'], candidatefile)

        cmd_detect = cmd_io + cmds['detect']
        print('tc command ', cmd_detect)
        # python version 2 or 3
        if int(sys.version[0]) == 2:
            sts = subprocess.check_output(cmd_detect, shell=True,
                                          stderr=subprocess.STDOUT,
                                          universal_newlines=True)
            print('sts ', sts)
        else:
            sts = subprocess.run(cmd_detect, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
            print(sts.stdout)
            if 'EXCEPTION' in sts.stdout:
                raise Exception('Stop')

        print('candidatefile ', candidatefile)

        tracking_period = period

        trackedfile = os.path.join(self.output_directory,
                                   'track_file_' + tracking_period + '_' + tracking_type + '.txt')

        # stitch candidates together
        cmd_stitch_io = '{} --in {} --out {} '.format(self.tc_stitch_script,
                                                      candidatefile,
                                                      trackedfile)
        cmd_stitch = cmd_stitch_io + cmds['stitch']
        print('tc command ', cmd_stitch)
        if int(sys.version[0]) == 2:
            sts = subprocess.check_output(cmd_stitch, shell=True,
                                          stderr=subprocess.STDOUT,
                                          universal_newlines=True)
            print('sts ', sts)
        else:
            sts = subprocess.run(cmd_stitch, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
            print('sts err ', sts.stderr)

        parameter_outfile = os.path.join(self.output_directory, metadata[
            'model'] + '_' + self.resolution_code + '_' + tracking_period + '_parameter_output_' + tracking_type + '.txt')
        self.write_parameter_file(parameter_outfile, cmd_detect, cmd_stitch)

        return candidatefile, trackedfile, filenames['pslfile'], tracking_period

    def collect_metadata(self, zg_available=False):
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
        metadata['psl_var'] = 'air_pressure_at_sea_level'
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

    def construct_command(self, nl_file, tracking_type):
        """
        Read the TempestExtreme command line parameters from the configuration.
        """
        commands = {}
        # for wibble in self.app_config.
        return commands

        # nl_dict = {}
        # with open(nl_file, 'r') as infile:
        #     for line in infile.readlines():
        #         if '&' + tracking_type in line:
        #             # section = line.strip().split(':')[-1][:-1]
        #             section = line.strip().split('_')[-1]
        #             nl_dict[section] = {}
        #         elif '=' in line:
        #             line_val = line.strip().split('=', 1)
        #             parameter = line_val[0]
        #             # strip the comma from the end of the namelist line
        #             value = line_val[1][:-1]
        #             nl_dict[section][parameter] = value
        #
        # cmd_step = {}
        # cmds = {}
        # for step in nl_dict:
        #     print('step ', step)
        #     cmd_step[step] = []
        #     for elt in nl_dict[step]:
        #         if nl_dict[step][elt] != '':
        #             cmd_step[step].append('--' + elt + ' ' + nl_dict[step][elt])
        #     cmds[step] = ' '.join(cmd_step[step])
        #
        # print('cmds ', cmds)
        # return cmds

    def data_files(self, *args):
        raise NotImplementedError

    def read_and_plot_tracks(self, *args, **kwargs):
        raise NotImplementedError

    def write_parameter_file(self, *args):
        raise NotImplementedError

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

    def _get_environment_variables(self):
        """
        Get the required environment variables from the suite. A list and
        explanation of the required environment variables is included in the
        documentation.
        """
        self.um_runid = os.environ['RUNID']
        self.cylc_task_cycle_time = os.environ['CYLC_TASK_CYCLE_TIME']