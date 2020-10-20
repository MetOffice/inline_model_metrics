#!/usr/bin/env bash
# run_tempest_tracking.sh
# set up enviromnent for TempestExtremes tracking (needs netcdf module load)
# then run it

set -eu

module load um_tools
module load scitools/production-os41-1
module load gcc/7.3.0
#module load cray-mpich/7.0.4
#module load hdf5/1.10.4/gnu/6.1.0
module load cray-netcdf/4.3.2

# Run script
python tempest_tracking.py "$@"
