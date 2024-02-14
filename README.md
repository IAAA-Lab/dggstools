# dggstools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10659071.svg)](https://doi.org/10.5281/zenodo.10659071)

Tools to manipulate raster and vector GIS data in the spatial framework provided by a DGGS.

**Requirements**
- Python 3.11 or higher
- `pip` for installing Python packages
- `GDAL` as dependency. 
- A minimum of 4 GB RAM (8 GB recommended for optimal performance)
- At least 10 GB of free disk space
- An active internet connection for downloading necessary Python packages
- Operating system: Windows, macOS, or Linux

Alternatively, we provide a Docker file for simplifying deployment and ensuring consistency across different development environments.

## Install and test instructions

`dggstools` uses the native `unittest` framework for testing.

All the necessary data for end-to-end testing is also included in the repository, ensuring you can validate the functionality of dggstools right after installation.


### Docker version

Build the image with:

```
docker build -t dggstools .
```

Run the tests with:

```
docker run --rm dggstools python -m unittest discover -s tests -p '*.py'
```

### Native python install 

To install `dggstools`, you can follow these steps:

1. Install the dependencies using pip and the provided `requirements.txt` file:

```
pip install -r requirements.txt
```

2. To install the `dggstools` package itself, run:

```
python setup.py install
```
