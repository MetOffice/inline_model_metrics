# inline_model_metrics
A Met Office tenten project Rose app to run model metrics inline with a climate model.

Documentation is temporarily available from: https://metoffice.github.io/inline_model_metrics/

The code consists of a preprocessor (currently specific to the UM output) which takes the input data from a source (model, reanalysis etc) and produces variables and filename suitable for use by the inline codes. 

There are initially two components of the inline tracking, tempest_cyclone and tempest_atmos_river. Tempest_cyclone uses the cyclone-tracking capability of TempestExtremes, which tempest_atmos_river uses TempestExtremes to track atmospheric rivers.

Other inline metrics code can be included here, will required changes to the pre-processing to provide the required input files/variables.

Currently the code saves and archives tracked output to the Met Office MASS storage. This needs to be split out into a postprocessor code to be more generic.
