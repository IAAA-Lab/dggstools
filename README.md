# DGGSTools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10659071.svg)](https://doi.org/10.5281/zenodo.10659071)

A python library and command line tool to manipulate raster and vector GIS data in the spatial framework provided by 
a DGGS (rHEALPix for now).

**Requirements**
- Python 3.10 or higher
- `pip` for installing Python packages

## Install the package and the command line tool
1. Download the latest package released on the GitHub Repository:

Get the .whl file from <https://github.com/IAAA-Lab/dggstools/releases/latest>.

2. Install the package and command line tool

```
pip install ./name-of-the-file-you-have-downloaded.whl
```

3. Now you can run dggstools in the command line:

```
dggstools --help
```

## Installation with Docker
We provide a Docker file for simplifying deployment and ensuring consistency across different development environments.

Build the image with:

```
docker build -t dggstools .
```

Run the tests with:

```
docker run --rm dggstools python -m unittest discover -s tests -p '*.py'
```

## Installation from sources 

To install `dggstools`, you can follow these steps right after you clone the Git repository (in its root directory):

1. Update the build tools:

```
pip install --upgrade pip          
pip install --upgrade build
pip install --upgrade wheel
```

2. Install the dependencies and build the `dggstools` package:

```
pip install .
```

3. Run the provided tests to see if everything is working (optional)
The package `dggstools` uses the `unittest` framework for testing. All the necessary data to run the tests are also included in the 
repository.

```
python -m unittest discover -s tests/data_tests -p '*.py'
```
