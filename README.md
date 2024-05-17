# DGGSTools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10659071.svg)](https://doi.org/10.5281/zenodo.10659071)

A python library (see the [latest version API docs](https://www.iaaa.es/dggstools/dggstools.html) and command line tool to manipulate raster and vector GIS data in the spatial framework provided by 
a DGGS (rHEALPix for now).

**Requirements**
- Python 3.10 or higher
- `pip` for installing Python packages

## Install the package and the command line tool
1. If you want the latest stable version (which is in PyPi):

```
pip install dggstools
```

2. If you prefer the latest non-stable version, you can download the latest package released on the GitHub Repository:

 - Get the .whl file from <https://github.com/IAAA-Lab/dggstools/releases/latest>.
 -  Install this wheel file:

 ```
 pip install ./name-of-the-file-you-have-downloaded.whl
 ```

3. In any case, once installed you can run dggstools in the command line:

```
dggstools --help
```

And you will see something like this:

```
Usage: dggstools [OPTIONS] COMMAND [ARGS]...                                                                                                                          
                                                                                                                                                                       
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --install-completion          Install completion for the current shell.                                                                                             │
│ --show-completion             Show completion for the current shell, to copy it or customize the installation.                                                      │
│ --help                        Show this message and exit.                                                                                                           │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ print-ras-rhpx-metadata  Takes a GeoTIFF file produced by dggstools and prints the metadata that dggstools stores in it. This metadata are necessary to store some  │
│                          rHEALPix system specific information.                                                                                                      │
│ print-vec-rhpx-metadata  Takes a GeoPackage file produced by dggstools and prints the metadata that dggstools stores in it. This metadata are necessary to store    │
│                          some rHEALPix system specific information, and some other information that can be useful if you want the original raster file back.        │
│ ras-rhpx-to-vec-rhpx     Transforms a rHEALPix GeoTIFF dataset produced by dggstools to a vector dataset in the GeoPackage format.                                  │
│ ras-to-rhpx-ras          Transforms a raster dataset in a common GIS format and reference system to a rHEALPix GeoTIFF. This includes: warping to the rHEALPix      │
│                          projection, resampling to one of the allowed rHEALPix resolutions (which depend on the rHEALPix system being used) and aligning to that    │
│                          rHEALPix grid.                                                                                                                             │
│ vec-ras-area-error       Takes a vector file and a rasterized rHEALPix version (as produced by the vec-to-rhpx-ras command) and:  - measures the area of each       │
│                          geometry in vector file;  - compares each of these areas with the areas of the cells which correspond to that geometry in the vector file. │
│                          This is an experimental, not thoroughly tested and barely documented command, and should be used just for testing purposes.                │
│ vec-rhpx-to-ras-rhpx     Transforms a vector dataset in rHEALPix produced by dggstools with the ras-rhpx-to-vec-rhpx command, to a raster dataset in GeoTIFF which  │
│                          is very similar to the one that was used as the original input to that operation.                                                          │
│ vec-to-rhpx-ras          Transforms a vector dataset with polygons in a common GIS format and reference system to a rHEALPix GeoTIFF. This GeoTIFF has a            │
│                          rasterization of the polygons following the constraints of the rHEALPix system (projection, valid resolution and grid alignment).          │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯ 
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
