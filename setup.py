# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.
import os
import setuptools


def get_long_description():
    """Use the contents of README.md as the long description"""
    with open("README.md", "r") as fh:
        return fh.read()


def extract_version():
    """
    Retrieve version information from the  __init__.py module.
    """
    version = ""
    directory = os.path.dirname(__file__)
    filename = os.path.join(directory, "inline_model_metrics", "__init__.py")

    with open(filename) as fd:
        for line in fd:
            line = line.strip()
            if line.startswith("__version__"):
                try:
                    version = line.split("=")[1].strip(" \"'")
                except Exception:
                    pass
                break

    if not version:
        print(f"WARNING: Unable to parse version information from file: {filename}")
        version = "0.0.0"

    return version


setuptools.setup(
    name="inline_model_metrics",
    packages=["inline_model_metrics"],
    version=extract_version(),
    license="BSD 3-Clause License",
    description=(
        "inline_model_metrics is a package to to run model metrics inline with a climate model."
    ),
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    author="Met Office",
    author_email="jon.seddon@metoffice.gov.uk",
    url="https://github.com/MetOffice/inline_model_metrics",
    download_url="https://github.com/MetOffice/inline_model_metrics/releases",
    keywords=["climate", "tracking", "inline"],
    install_requires=["metoffice-afterburner"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
