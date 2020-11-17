#!/bin/bash
# (C) British Crown Copyright 2020, Met Office.
# Please see LICENSE for license details.
#
# SYNOPSIS
#
#   copy_code_to_suite.sh <suite_path>
#
# DESCRIPTION
#
#   Copy the Python code to the tempest_tracker app in the specified rose suite.
#   The existing Python code directory (app/tempest_tracker/file/tempest_tracker)
#   in the suite is first deleted.
#
#   Users must remember to run `rose suite-run --reload` to also copy the new
#   code to the host that the tracking is being run on.
#
# ARGUMENTS
#
#   suite_path
#      The path on the local machine that the Rose suite that is running the
#      tracking was checked out to.

SUITE_DIR=$1
if [ ! -d "$SUITE_DIR" ]; then
  echo "$SUITE_DIR does not exit"
  exit 1
fi

# Get current directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Get the directory that the code is thereofore from
REPO_DIR=$(readlink -f "$DIR/..")

SRC_DIR=$REPO_DIR/rose-suite/app/tempest_tracker/file/tempest_tracker
DEST_PARENT=$SUITE_DIR/app/tempest_tracker/file
DEST_DIR=$DEST_PARENT/tempest_tracker

if [ -d "$DEST_DIR" ]; then
  rm -rf "$DEST_DIR"
else
  echo "No existing $DEST_DIR to remove"
fi

if [ ! -d "$DEST_PARENT" ]; then
  mkdir -p "$DEST_PARENT"
fi
cp -r "$SRC_DIR" "$DEST_PARENT"

echo ""
echo "Don't forget to run \"rose suite-run --reload\" to copy the changes to \
the host that the tracking is being run on."