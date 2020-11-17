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

Git workflow
############

#. Setup your SSH keys (https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/connecting-to-github-with-ssh)
#. Obtain a local copy of the suite::

      git clone git@github.com:MetOffice/tenten_tempestextremes.git

#. Change directory into the repository::

      cd tenten_tempestextremes

#. If you have previously cloned the repository then update your local copy to
   the latest version::

      git checkout main
      git pull origin main

#. Create a branch for your changes and change to this branch::

      git fetch origin
      git branch <branch-name> origin/main
      git checkout <branch-name>

#. When you're happy with your changes, check the changes that you've made::

      git status
      git diff

#. You can then commit these changes and push them back up to GitHub::

      git commit -am '<commit message>'
      git push origin <branch-name>

#. You can make as many commits and pushes as you want.
#. Create a pull request at https://github.com/MetOffice/tenten_tempestextremes/compare
   by in the compare box selecting your branch and clicking on "Create pull request".

#. Once the code has been reviewed then the pull request can be merged using the
   GitHub web page for that pull request.

#. After merging, change your local copy of the code back to the main branch, delete
   your local copy of the development branch and pull in the changes from GitHub's
   main branch to your local copy::

      git checkout main
      git branch -D  <branch-name>
      git pull origin main

To examine complex changes, you might also want to consider using graphical diff
rather than just the command line, in which case you should add the following lines
to your `~/.gitconfig` file::

   [diff]
        tool = tkdiff
   [difftool]
        prompt = False

and then, rather than just `git diff`, you can use::

   git difftool


Running the tests
#################

#. Load a Python environment in the usual way for your location.
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

Building the Documentation
##########################

In the checked out repository, make sure that your Python environment includes
Sphinx (standard scientific ones do)::

   export PYTHONPATH=/path/to/tenten_tempestextremes/rose-suite/app/tempest_tracker/file
   cd docs
   make clean && make html

The build documentation can then be viewed in your browser (replace Firefox
with the name of your browser if required)::

   firefox build/html/index.html

To add new track types
######################

#. Add an option for the new track type to the values in `[common=track_type]` in
   `rose-meta.conf`.
#. Select the new track type in the `track_type` setting in the `[common]` section
   of `rose-app.conf`.
#. Add `[<track_type>_detect]`  and `[<track_type>_stitch]` sections to `rose-app.conf`.