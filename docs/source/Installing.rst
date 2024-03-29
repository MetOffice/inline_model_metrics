Installing
==========

The following code must already be installed on the host that the tracking will
be run on:

* afterburner (https://code.metoffice.gov.uk/trac/afterburner/)
* tempest_helper (https://github.com/MetOffice/tempest_helper/)
* inline_model_metrics (this software) (https://github.com/MetOffice/inline_model_metrics)
* TempestExtremes

Then perform the following steps:

#. If you want to plot the tracks, ensure that the Cartopy Natural Earth files
   are installed into `~/.local/share/cartopy/shapefiles/natural_earth/physical/`
   in your home directory on the host that the tracking will be run on, if that
   host doesn't have access to the Internet to download these files itself.
#. Checkout the code.
#. Copy the Rose app directory into your Rose suite::

      cp -r rose-suite/app/tempest_tracker ~/roses/<suite-id>/app/

#. Ensure that the  `AFTERBURNER_HOME_DIR` setting in the `app/tempest_tracker/rose-app.conf`
   file in your suite points to your Afterburner installation on the host that the tracking will
   be run on.
#. Ensure that the  `PYTHONPATH` setting in the `app/tempest_tracker/rose-app.conf`
   file in your suite points to your tempest_helper and inline_model_metrics installations
   on the host that the tracking will be run on.
#. Ensure that the `tc_detect_script`, `tc_stitch_script` and `tc_editor_script` settings
   in the `app/tempest_tracker/rose-app.conf` file in your suite points to your
   TempestExtremes installation on the host that the tracking will be run on.
#. Update the `orography_file` setting in the `app/tempest_tracker/rose-app.conf`
   file to refer to a suitable orography file for this simulation.
   TempestExtremes installation on the host that the tracking will be run on.
#. Modify `suite.rc` in your suite to add the changes shown in `rose-suite/suite_rc.diff`.

