# inline_model_storms
A Met Office tenten project Rose app to run model storms inline (or postprocess) with a climate model.

Documentation is available from: https://inline-model-storms.readthedocs.io/

The code consists of:
  within inline_model_storms:
    tempest_common.py - common code used by other parts of package
    tempest_cyclone.py - code to execute parts of the TempestExtremes detection and tracking code
    tempest_atmos_river.py - code to execute the TempestExtremes atmospheric river detection
    um_preprocess.py - code to preprocess the data obtained from the data archive, to format it for TempestExtremes using the namelists given in the suite apps
    um_postprocess - code for postprocessing, data archiving and tidying up
    load_trajectories.py - code to load the files produced by TempestExtremes
    save_trajectories.py - code to save the tracked files as netcdf
    trajectory_manipulations.py - routines to manipulate the tracks produced

together with a rose/cylc suite containing apps:
  tempest_get_data - suite-specific code to get data from Met Office MASS archive, including STASH codes etc within the file/ directory.
  tempest_preprocess - namelists to preprocess the data after getting it, to make it compatible with TempestExtremes with regards to variable names etc
  tempest_tracker - namelists to run various components of TempestExtremes with given input parameters
  tempest_atmosriver - namelists to run atmospheric river identification
  tempest_postprocess - namelists to run the postprocessing (data archiving, tidying up)

Currently the code saves and archives tracked output to the Met Office MASS storage. This code would need altering for other platforms.
