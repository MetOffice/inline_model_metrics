Developers' Instructions
========================

Developing the code is complicated by the fact that this is a single Rose app
that will be used in many Rose suites. Care needs to be taken so that any changes
that are made in one suite are committed to the repository so that they can be
re-used in other suites.

Editing source code
###################

The safe method:

#. Make all changes to the code in your checked out copy of https://github.com/MetOffice/tenten_tempestextremes/
#. Copy the changes to the suite by running::

      bin/copy_code_to_suite.sh ~/roses/<suite-id>

#. Reload the suite to pick up the changes::

      rose suite-run --reload

#. Re-run the tracking task in the Cylc GUI, look at the output from the suite and
   go back to step 1 to make further changes.
#. When everything's working, commit and push the changes to GitHub.

The dangerous way:

#. Make your changes in the suite's work directory on the host where the tracking's
   running.
#. When everything's working, copy the changes back from the host where the tracking's
   running back to your local copy of https://github.com/MetOffice/tenten_tempestextremes/
   and commit and push the changes.

Running the Tests
#################

#. Load a Python environment in the typical way for your location.
#. Point to the Afterburner installation on your local machine and add this to your
   `PYTHONPATH`::

      export AFTERBURNER_HOME_DIR=/path/to/afterburner/software/turbofan/current
      export PYTHONPATH=$AFTERBURNER_HOME_DIR/lib/python:$PYTHONPATH

#. Download a local copy of tempest_helper (https://github.com/MetOffice/tempest_helper/)
   and add it to your `PYTHONPATH`::

      export PYTHONPATH=/path/to/tempest_helper:$PYTHONPATH

#. Add the tenten_tempestextremes code to your `PYTHONPATH`::

      export PYTHONPATH=/path/to/tenten_tempestextremes/rose-suite/app/tempest_tracker/file:$PYTHONPATH

#. Run the tests::

      pytest -vv


To add new track types
######################

#. Add an option for the new track type to the values in `[common=track_type]` in
   `rose-meta.conf`.
#. Select the new track type in the `track_type` setting in the `[common]` section
   of `rose-app.conf`.
#. Add `[<track_type>_detect]`  and `[<track_type>_stitch]` sections to `rose-app.conf`.