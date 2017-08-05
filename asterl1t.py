# Copyright (c) 2017, Bruno Ruas de Pinho
# brunorpinho10@gmail.com

import pandas as pd
import gdal
import rasterio
from rasterio.warp import reproject, Resampling
import datetime
import pyproj
import numpy as np

ucc = np.matrix([  
    [0.676,  1.688,    2.25,   np.nan],
    [0.708,  1.415,    1.89,   np.nan],
    [0.423,  0.862,    1.15,   np.nan],
    [0.423,  0.862,    1.15,   np.nan],
    [0.1087, 0.2174,   0.2900, 0.2900],
    [0.0348, 0.0696,   0.0925, 0.4090],
    [0.0313, 0.0625,   0.0830, 0.3900],
    [0.0299, 0.0597,   0.0795, 0.3320],
    [0.0209, 0.0417,   0.0556, 0.2450],
    [0.0159, 0.0318,   0.0424, 0.2650],
    [np.nan, 0.006822, np.nan, np.nan],
    [np.nan, 0.006780, np.nan, np.nan],
    [np.nan, 0.006590, np.nan, np.nan],
    [np.nan, 0.005693, np.nan, np.nan],
    [np.nan, 0.005225, np.nan, np.nan],
])

ucc_header = ['HGH', 'NOR', 'LG1', 'LG2']
band_index = ['1', '2', '3B', '3N', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13' ,'14']
ucc_dataframe = pd.DataFrame(ucc, index=band_index, columns=ucc_header)
ucc_dict = ucc_dataframe.transpose().to_dict()

res = np.concatenate([[15] * 4, [30] * 6, [90] * 5])
resolution_dict = dict(zip(band_index, res))

index_dict = dict(zip(band_index, np.arange(0, 15)))

irradiance_dict = {
    '1': 1848.99,
    '2': 1555.74,
    '3N': 1119.47,
    '4': 231.25,
    '5': 79.81,
    '6': 74.99,
    '7': 68.66,
    '8': 59.74,
    '9': 56.92
}

def process_aster_dataset(aster_file, 
                          destination,
                          dst_transform,
                          dst_crs,
                          resampling=Resampling.bilinear,
                          return_radiance=False):
    
    '''Calculates reflectance of "ASTER_L1T .hdf" files and import results into a 3-dimensional array using rasterio.
    
    All bands are reprojected to the desired shape.
    Source shape and transform is calculated from the input file.
    
    This is based on this official rasterio example:
        https://github.com/mapbox/rasterio/blob/master/docs/topics/reproject.rst 
        
    
    Parameters
    ----------
    
    aster_file: string or file
        A path to the ASTER L1T .hdf file.
        
    destination: 2d array
        The destination is a 2D array of shape MxN. This is used to
        specify the shape of the projected ASTER data.
    
    dst_transform: affine.Affine()
        Target affine transformation. This can be easily obtained using
        the rasterio.transform functions.
        
    dst_crs: CRS or dict
        Target coordinate reference system. Example: {'init': u'epsg:4326'}
        for WGS84 reference system.
        
    resampling: int
        Resampling method to use from rasterio.warp.Resampling.
            One of the following:
                Resampling.nearest,
                Resampling.bilinear,
                Resampling.cubic,
                Resampling.cubic_spline,
                Resampling.lanczos,
                Resampling.average,
                Resampling.mode    
        
    return_radiance: Boolean
        Return radiance (True) or reflectance (False).
        Default: False
    
    
    Returns
    ---------
    
    output: ndarray with shape MxNxB - reprojected reflectance
    
    Where:
    
    - M is the number of rows;
    - N is te number of columns;
    - B is the number of bands (9 or 14 depending of return_radiance):
    
    numpy B index |  0    1    2     3     4    5    6    7    8    9    10    11    12    13    14
    ------------------------------------------------------------------------------------------------
     ASTER  band  | '1', '2', '3B', '3N', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13' ,'14'
     
    
    Reference
    ---------
    Finn, Michael P., Matthew D. Reed, and Kristina H. Yamamoto. 
        "A straight forward guide for processing radiance and reflectance for EO-1 ALI, 
        Landsat 5 TM, Landsat 7 ETM+, and ASTER." Unpublished Report from USGS/Center 
        of Excellence for Geospatial Information Science 8 (2012).
    
        https://cegis.usgs.gov/soil_moisture/pdf/A%20Straight%20Forward%20guide%20for%20Processing%20Radiance%20and%20Reflectance_V_24Jul12.pdf
    
    
    Info about ASTER Level 1T products:
        https://lpdaac.usgs.gov/dataset_discovery/aster/aster_products_table/ast_l1t
        
    To download ASTER Level 1T products:
        https://earthexplorer.usgs.gov/ 
        
    '''

    aster = gdal.Open(aster_file)
    aster_sds = aster.GetSubDatasets()
    
    if return_radiance:
        aster_array = np.expand_dims(destination, 2).repeat(15, 2).copy()
    else:
        aster_array = np.expand_dims(destination, 2).repeat(10, 2).copy()

    wgs84_22s = pyproj.Proj({'init': u'epsg:32722'})

    solar_direction = np.array(aster.GetMetadata()['SOLARDIRECTION'].split(', '), dtype=np.float32)
    gain_dict = dict([g.replace('0', '').split('=')[1].split(', ') for g in aster.GetMetadata_List() if 'GAIN' in g])
    [gain_dict.__setitem__(b, 'NOR') for b in ['10', '11', '12', '13', '14']]

    ul = np.asarray(aster.GetMetadata()['UPPERLEFT'].split(', ')[::-1], dtype=np.float64)
    lr = np.asarray(aster.GetMetadata()['LOWERRIGHT'].split(', ')[::-1], dtype=np.float64)
    
    calendardate = aster.GetMetadata()['CALENDARDATE']
    yr,mo,day = [np.uint8(s) for s in [calendardate[0:4], calendardate[4:6], calendardate[6:]]]
    doy = datetime.datetime(yr, mo, day).timetuple().tm_yday

    d = (1 - 0.01672 * np.cos( np.radians( 0.9856 * ( doy - 4 ) ) ) )

    raster_list = [i[0] for i in aster_sds if 'ImageData' in i[0]]
    
    for raster in raster_list:

        band_str = raster.split('ImageData')[-1]
        index = index_dict[band_str]
        gain = gain_dict[band_str]
        
        print ('File: {} - Band: {}'.format(aster_file, band_str, aster_array.shape[-1]))        
        
        with rasterio.open(raster) as src:

            source = src.read(1)

            src_shape = src.shape
            src_crs = {'init': u'epsg:4326'} # WGS 84
            
            src_transform = rasterio.transform.from_bounds(ul[0], lr[1], lr[0], ul[1], src_shape[1], src_shape[0])

            dst = destination.copy()

            reproject(
                source,
                dst,
                src_transform=src_transform,
                src_crs=src_crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                resampling=resampling
            )

            dst = np.ma.masked_equal(dst, 0).filled(np.nan)
            dst = (dst - 1) * ucc_dict[band_str][gain]
            
            if return_radiance:

                aster_array[:,:, index] = dst #radiance
                
            else:
            
                if band_str in band_index[:10]:

                    reflectance_dst = (np.pi * dst * (d**2)) / (irradiance_dict[band_str] * np.sin(np.deg2rad(solar_direction[1])))

                    aster_array[:,:, index] = reflectance_dst  #reflectance
            
    return aster_array
