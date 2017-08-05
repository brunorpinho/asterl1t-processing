
A script to calculate reflectance from "ASTER_L1T .hdf" files. It also import the results into a 3-dimensional array using rasterio.
    
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
