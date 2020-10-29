Installing
==========

The following code must already be installed and included in your path:

* afterburner (https://code.metoffice.gov.uk/trac/afterburner/)

Then perform the following steps:

#. Checkout the code
#. Copy the Rose app directory into your Rose suite:
`cp -r rose-suite/app/tempest_tracker ~/roses/<suite-id>/app/`
#. Modify `suite.rc` to add the changes shown in `rose-suite/suite_rc.diff`.
#. Add the path to the code to the `PYTHONPATH` in `rose-app.conf`
