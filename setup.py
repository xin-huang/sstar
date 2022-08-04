import os.path
from setuptools import setup

# The directory containing this file
HERE = os.path.abspath(os.path.dirname(__file__))

# The text of the README file
with open(os.path.join(HERE, "README.md")) as fid:
    README = fid.read()

# This call to setup() does all the work
setup(
    name="sstar",
    version="1.1.2",
    description="A Python package for detecting archaic introgression from population genetic data with S*",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/xin-huang/sstar",
    author="Xin Huang",
    author_email="xinhuang.res@gmail.com",
    license="Apache-2.0",
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
    packages=["sstar"],
    include_package_data=True,
    install_requires=[
        "demes",
        "numpy",
        "pandas",
        "rpy2",
        "scikit-allel",
        "scipy",
    ],
    entry_points={"console_scripts": ["sstar=sstar.__main__:main"]},
)
