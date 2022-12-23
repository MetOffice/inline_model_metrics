# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
"""
Tests inline_model_metrics.common
"""
import subprocess

import unittest

from inline_model_metrics.common import run_cmd


class TestRunCmd(unittest.TestCase):
    """Tests the inline_model_metrics.common.run_cmd"""
    def test_success(self):
        sts = run_cmd("true")
        self.assertEqual(sts.returncode, 0)

    def test_fails_check_false(self):
        sts = run_cmd("false", check=False)
        self.assertEqual(sts.returncode, 1)

    def test_fails_check_true(self):
        self.assertRaises(subprocess.CalledProcessError, run_cmd, "false")

    def test_stderr_check_false_ignored(self):
        sts = run_cmd(">&2 echo 'bananas'", check=False)
        self.assertEqual(sts.returncode, 0)
        self.assertEqual(sts.stderr, "bananas\n")

    def test_stderr_check_true_exception(self):
        self.assertRaisesRegex(
            RuntimeError,
            "Error found in cmd output bananas",
            run_cmd,
            ">&2 echo 'bananas'",
            check=True
        )
