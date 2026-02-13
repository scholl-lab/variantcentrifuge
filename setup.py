# File: variantcentrifuge/setup.py
# Location: variantcentrifuge/variantcentrifuge/setup.py
"""
Setup script for variantcentrifuge.

This file configures how the package is built, installed, and what
dependencies are required.
"""

import os

from setuptools import find_packages, setup

# Load version from version.py without importing the module
version = {}
with open(os.path.join("variantcentrifuge", "version.py")) as f:
    exec(f.read(), version)

# Read the README for the long description
this_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_dir, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="variantcentrifuge",
    version=version["__version__"],
    description="A tool to filter, extract, and analyze variants from VCF files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Bernt Popp",
    author_email="bernt.popp.md@gmail.com",
    url="https://github.com/scholl-lab/variantcentrifuge",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "pandas",
        "jinja2",
        "scipy",
        "statsmodels",
        "numpy",
        "intervaltree",
        "psutil",
        "smart-open",
        "openpyxl",
    ],
    extras_require={
        "dev": [
            "pytest",
            "pytest-cov",
            "pytest-mock",
            "black==24.4.0",
            "flake8==7.0.0",
            "flake8-docstrings",
            "isort==5.13.2",
            "pre-commit",
        ],
    },
    entry_points={"console_scripts": ["variantcentrifuge=variantcentrifuge.cli:main"]},
    include_package_data=True,
    package_data={"variantcentrifuge": ["config.json"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
)
