# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
"""
Tests the tempest_tracker.TempestTracker Afterburner app
"""
import os
import shutil
import tempfile
import unittest

from inline_model_metrics import TempestExtremesCyclone
from inline_model_metrics.tempest_common import _is_date_after, _is_date_after_or_equal


class TestTempestTracker(unittest.TestCase):
    """Tests the inline_model_metrics.TempestExtremesCyclone Afterburner app"""

    def setUp(self):
        self.runtime_dir = tempfile.mkdtemp()
        _fd, self.cfg_file = tempfile.mkstemp(suffix=".conf", dir=self.runtime_dir)
        os.environ["RUNID"] = "ab123"
        os.environ["CYLC_SUITE_NAME"] = "u-ab123"
        os.environ["CYLC_TASK_CYCLE_TIME"] = "20000101T0000Z"
        os.environ["NEXT_CYCLE"] = "20000201T0000Z"
        os.environ["TIME_CYCLE"] = "20000101T0000Z"
        os.environ["TM2_CYCLE"] = "19991101T0000Z"
        os.environ["PREVIOUS_CYCLE"] = "19991201T0000Z"
        os.environ["STARTDATE"] = "19700101T0000Z"
        os.environ["LASTCYCLE"] = "20191201T0000Z"
        os.environ["ENDDATE"] = "20200101T0000Z"
        os.environ["IS_LAST_CYCLE"] = "true"
        os.environ["NCODIR"] = "/some/dir"
        self.basic_app_config = """
        [common]
        resolution=N96
        input_directory=${DATAM}
        output_directory=${DATAM}/tempest_tracking
        nodeedit_vars=["psl", "uas", "vas"]
        regrid_resolutions=["native"]
        tc_detect_script=/home/d05/hadom/tempestextremes-master/bin/DetectNodes
        tc_stitch_script=/home/d05/hadom/tempestextremes-master/bin/StitchNodes
        track_types=['tc_psl']
        tc_variables_input=["psl", "uas"]
        tc_variables_rename=["psl", "uas"]
        um_file_pattern="{runid}a.{stream}{date_start}_{variable}.nc"
        file_pattern_processed = "{variable}_{frequency}_{runid}_{date_start}-{date_end}.nc"
        in_fmt_stitch_default="lon,lat,psl_min,sfcWind_max"
        out_fmt_profile1_default="lon,lat,psl_min,sfcWind_max"
        out_fmt_profile2_default="lon,lat,psl_min"

        [tc_psl_detect]
        searchbymin="air_pressure_at_sea_level"
        mergedist=6.0
        closedcontourcmd="air_pressure_at_sea_level,200.0,5.5,0;_DIFF(zg(0),zg(2)),-6.0,6.5,1.0"
        outputcmd="air_pressure_at_sea_level,min,0;_VECMAG(x_wind,y_wind),max,2;_DIFF(zg(0),zg(2)),min,6.5;surface_altitude,max,0"
        latname="latitude"
        lonname="longitude"
        verbosity=

        [tc_psl_stitch]
        in_fmt="in_fmt_stitch_default"
        min_endpoint_dist=8.0
        
        [tc_psl_profile]
        in_fmt="lon,lat,psl,wind10m,zgdiff,surface_altitude"
        out_fmt="out_fmt_profile1_default"
        """  # noqa

    def tearDown(self):
        if os.path.isdir(self.runtime_dir):
            shutil.rmtree(self.runtime_dir, ignore_errors=True)

    def test_construct_command(self):
        """Test that the command dict is created from the config file"""

        _create_app_config_file(self.cfg_file, self.basic_app_config)
        args = ["-c", self.cfg_file, "-q"]
        app = TempestExtremesCyclone(args)
        app._get_app_options()
        app._get_environment_variables()
        commands = app._construct_command("tc_psl")
        expected = {
            "detect": '--closedcontourcmd "air_pressure_at_sea_level,200.0,5.5,'
            '0;_DIFF(zg(0),zg(2)),-6.0,6.5,1.0" '
            '--latname "latitude" '
            '--lonname "longitude" '
            "--mergedist 6.0 "
            '--outputcmd "air_pressure_at_sea_level,min,0;_VECMAG('
            "x_wind,y_wind),max,2;_DIFF(zg(0),zg(2)),min,6.5;surface"
            '_altitude,max,0" '
            '--searchbymin "air_pressure_at_sea_level"',
            "stitch": '--in_fmt "lon,lat,psl_min,sfcWind_max" --min_endpoint_dist 8.0',
            "profile": '--in_fmt "lon,lat,psl,wind10m,zgdiff,surface_altitude" '
            '--out_fmt "lon,lat,psl_min,sfcWind_max"',
            "detectblobs": None,
            "nodefilefilter": None,
        }
        self.maxDiff = None  # Show the full diff if test fails
        self.assertEqual(expected, commands)

    def test_is_date_after_one_millenium(self):
        """Test _is_date_after()"""
        self.assertTrue(_is_date_after("20000101T0000Z", "19700101T0000Z"))

    def test_is_date_after_one_year(self):
        """Test _is_date_after()"""
        self.assertTrue(_is_date_after("19700101T0000Z", "19690101T0000Z"))

    def test_is_date_after_one_day(self):
        """Test _is_date_after()"""
        self.assertTrue(_is_date_after("19700102T0000Z", "19700101T0000Z"))

    def test_is_date_after_not_one_day(self):
        """Test _is_date_after()"""
        self.assertFalse(_is_date_after("19700101T0000Z", "19700102T0000Z"))

    def test_is_date_after_not_one_year(self):
        """Test _is_date_after()"""
        self.assertFalse(_is_date_after("19700101T0000Z", "19710101T0000Z"))

    def test_is_date_after_or_equal_equal(self):
        """Test _is_date_after()"""
        self.assertTrue(_is_date_after_or_equal("19700101T0000Z", "19700101T0000Z"))

    def test_is_date_after_or_equal_equal_day_diff_hour(self):
        """Test _is_date_after()"""
        self.assertTrue(_is_date_after_or_equal("19700101T0000Z", "19700101T0100Z"))

    def test_is_date_after_or_equal_one_year(self):
        """Test _is_date_after()"""
        self.assertTrue(_is_date_after_or_equal("19700101T0000Z", "19690101T0000Z"))

    def test_is_date_after_or_equal_not_one_year(self):
        """Test _is_date_after()"""
        self.assertFalse(_is_date_after_or_equal("19690101T0000Z", "19700101T0000Z"))


def _create_app_config_file(config_file, config_text):
    with open(config_file, "w") as fh:
        fh.writelines([line.strip() + "\n" for line in config_text.split("\n")])


if __name__ == "__main__":
    unittest.main()
