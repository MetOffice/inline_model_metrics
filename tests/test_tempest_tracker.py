# (C) British Crown Copyright 2020, Met Office.
# Please see LICENSE for license details.
"""
Tests the tempest_tracker.TempestTracker Afterburner app
"""
import os
import shutil
import tempfile
import unittest

from tempest_tracker import TempestTracker


class TestTempestTracker(unittest.TestCase):
    """Tests the tempest_tracker.TempestTracker Afterburner app"""
    def setUp(self):
        self.runtime_dir = tempfile.mkdtemp()
        _fd, self.cfg_file = tempfile.mkstemp(suffix='.conf',
                                              dir=self.runtime_dir)
        os.environ['RUNID'] = 'ab123'
        os.environ['CYLC_TASK_CYCLE_TIME'] = '20000101T0000Z'
        self.basic_app_config = """
        [common]
        resolution=N96
        input_directory=${DATAM}
        output_directory=${DATAM}/tempest_tracking
        tc_detect_script=/home/d05/hadom/tempestextremes-master/bin/DetectNodes
        tc_stitch_script=/home/d05/hadom/tempestextremes-master/bin/StitchNodes
        track_type=tc_slp

        [tc_slp_detect]
        searchbymin="air_pressure_at_sea_level"
        mergedist=6.0
        closedcontourcmd="air_pressure_at_sea_level,200.0,5.5,0;_DIFF(zg(0),zg(2)),-6.0,6.5,1.0"
        outputcmd="air_pressure_at_sea_level,min,0;_VECMAG(x_wind,y_wind),max,2;_DIFF(zg(0),zg(2)),min,6.5;surface_altitude,max,0"
        latname="latitude"
        lonname="longitude"
        verbosity=
        
        [tc_slp_stitch]
        min_endpoint_dist=8.0
        """

    def tearDown(self):
        if os.path.isdir(self.runtime_dir):
            shutil.rmtree(self.runtime_dir, ignore_errors=True)

    def test_collect_metadata_with_zg(self):
        """Test that the metadata dict's populated from the config"""
        _create_app_config_file(self.cfg_file, self.basic_app_config)
        args = ['-c', self.cfg_file, '-q']
        app = TempestTracker(args)
        app._get_app_options()
        app._get_environment_variables()
        metadata = app._collect_metadata(True)
        expected = {
            'model': 'ab123',
            'resol': 'N96',
            'model_name': 'ab123',
            'u_var': 'x_wind',
            'v_var': 'y_wind',
            't_var': 'unknown',
            'lat_var': 'lat',
            'lon_var': 'lon',
            'topo_var': 'surface_altitude'
        }
        self.assertEqual(expected, metadata)

    def test_collect_metadata_without_zg(self):
        """Test that the metadata dict's populated from the config"""
        _create_app_config_file(self.cfg_file, self.basic_app_config)
        args = ['-c', self.cfg_file, '-q']
        app = TempestTracker(args)
        app._get_app_options()
        app._get_environment_variables()
        metadata = app._collect_metadata(False)
        expected = {
            'model': 'ab123',
            'resol': 'N96',
            'model_name': 'ab123',
            'u_var': 'x_wind',
            'v_var': 'y_wind',
            't_var': 'ta',
            'lat_var': 'lat',
            'lon_var': 'lon',
            'topo_var': 'surface_altitude'
        }
        self.assertEqual(expected, metadata)

    def test_construct_command(self):
        """Test that the command dict is created from the config file"""

        _create_app_config_file(self.cfg_file, self.basic_app_config)
        args = ['-c', self.cfg_file, '-q']
        app = TempestTracker(args)
        app._get_app_options()
        app._get_environment_variables()
        commands = app._construct_command()
        expected = {
            'detect': '--closedcontourcmd "air_pressure_at_sea_level,200.0,5.5,'
                      '0;_DIFF(zg(0),zg(2)),-6.0,6.5,1.0" '
                      '--latname "latitude" '                      
                      '--lonname "longitude" '
                      '--mergedist 6.0 '
                      '--outputcmd "air_pressure_at_sea_level,min,0;_VECMAG('
                      'x_wind,y_wind),max,2;_DIFF(zg(0),zg(2)),min,6.5;surface'
                      '_altitude,max,0" '
                      '--searchbymin "air_pressure_at_sea_level"',
            'stitch': '--min_endpoint_dist 8.0'
        }
        self.maxDiff = None  # Show the full diff if test fails
        self.assertEqual(expected, commands)


def _create_app_config_file(config_file, config_text):
    with open(config_file, 'w') as fh:
        fh.writelines([line.strip()+'\n' for line in config_text.split('\n')])


if __name__ == '__main__':
    unittest.main()