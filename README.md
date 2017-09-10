[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.888270.svg)](https://doi.org/10.5281/zenodo.888270)

A script to calculate reflectance and radiance from ASTER L1T images (hdf files). It also import the results into a 3-dimensional array using rasterio.
    
All bands are reprojected to a desired shape and CRS.

Source shape and transform are calculated directly from the input file.

## Dependencies

For now the dependencies are numpy, pandas, GDAL, rasterio and pyproj

## Note

Please, consider donation if you are profiting from this script.
author: Bruno Ruas de Pinho - brunorpinho10@gmail.com

## TODO

* Clean up the code
* Transform into a Python package

## Cite As

De Pinho, Bruno. (2017, September 10). brunorpinho/asterl1t-processing: A Python script to calculate reflectance and radiance from ASTER L1T images (hdf files). Zenodo. http://doi.org/10.5281/zenodo.888270
