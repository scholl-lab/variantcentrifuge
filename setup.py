# File: setup.py
# Location: variantcentrifuge/setup.py

"""
Setup script for variantcentrifuge.

This file configures how the package is built, installed, and what
dependencies are required.
"""

import os
from setuptools import setup, find_packages

# Read the README for the long description
this_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_dir, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="variantcentrifuge",
    version="0.1.0",
    description="A tool to filter, extract, and analyze variants from VCF files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/scholl-lab/variantcentrifuge",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "variantcentrifuge=variantcentrifuge.cli:main"
        ]
    },
    include_package_data=True,  # Ensure package data specified in MANIFEST.in is included
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
)