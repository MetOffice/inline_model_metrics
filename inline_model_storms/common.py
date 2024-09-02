# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import logging
import subprocess

logger = logging.getLogger(__name__)


def run_cmd(cmd, check=True):
    """
    Run the command 'cmd' in a shell.

    :param str cmd: The command to run.
    :param bool check: Raise an exception if there's a non-zero return code.
    :rtype: subprocess.CompletedProcess
    :return: The subprocess.run() return object
    """
    sts = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=check,
    )
    if check and sts.stderr:
        if "Warning" in sts.stderr:
            msg = f"Warning found in cmd output {sts.stderr}"
            logger.warning(msg)
        else:
            msg = f"Error found in cmd output {sts.stderr}"
            raise RuntimeError(msg)
    return sts
