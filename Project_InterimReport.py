#!/usr/bin/env python
"""
ABE651 (Spring 2021)
Project | Interim report
Edited by Ana Morales

"""

# import required modules
import rasterio
from rasterio import plot
import matplotlib.pyplot as plt
import os
import geopandas as gpd
import earthpy as et
import earthpy.plot as ep
from shapely.geometry import mapping, shape
import fiona
from pyproj import Proj, transform
from osgeo import gdal
import numpy as np
from rasterio.plot import show_hist
from rasterio.plot import show
import rioxarray as rxr
import rasterstats
from matplotlib.ticker import FormatStrFormatter #For the 2 decimals in plots
from rasterstats import zonal_stats

# define function
def stack_bands(folderPath, rootname, stackBands):
    '''
    The function takes three arguments:
        1) folderPath: path in which the individual layers are available
        2) rootname: each individual layer has a rootname that is common to all the bands
        3) stackBands: path and name of the output stacked layer
    '''   
    #import bands as separate 1 band raster
    band2 = rasterio.open(folderPath+rootname+'_B02.jp2', driver='JP2OpenJPEG') #blue
    band3 = rasterio.open(folderPath+rootname+'_B03.jp2', driver='JP2OpenJPEG') #green
    band4 = rasterio.open(folderPath+rootname+'_B04.jp2', driver='JP2OpenJPEG') #red
    band8 = rasterio.open(folderPath+rootname+'_B08.jp2', driver='JP2OpenJPEG') #nir
    
    #stack bands 
    with rasterio.open(stackBands,'w',driver='Gtiff', width=band4.width, height=band4.height, 
                       count=4,crs=band4.crs,transform=band4.transform, dtype=band4.dtypes[0]) as bgrn:
        bgrn.write(band2.read(1),1) 
        bgrn.write(band3.read(1),2) 
        bgrn.write(band4.read(1),3)
        bgrn.write(band8.read(1),4) 
        bgrn.close()
   
    return stackBands

def info_raster(rasterPath):
    '''
    The function takes one argument, which is the path to the raster of interest.
    '''
    #Opening the raster of interest
    with rasterio.open(rasterPath,'r') as raster:
        #Coordinate reference system
        raster_crs= raster.crs
        print("The coordinate reference system for {} is {} \n".format(rasterPath, raster_crs))
    
        #number of raster bands
        raster_bands= raster.count
        print("The number of bands for {} is {} \n".format(rasterPath, raster_bands))
    
        #number of raster columns
        raster_columns= raster.width
        print("The width (number of columns) of this image {} is {} \n".format(rasterPath, raster_columns))

        #number of raster rows
        raster_rows= raster.height
        print("The height (number of rows) of this image {} is {} \n".format(rasterPath, raster_rows))
        
        raster.close()


def info_vector(polygonFile):
    '''
    The function takes one argument, which is the path to the shapefile (vector) of interest.
    '''
    #Opening the vector of interest using geopandas
    vector = gpd.read_file(polygonFile)
    
    #Coordinate reference system
    vector_crs= vector.crs
    print("The coordinate reference system for {} is {} \n".format(polygonFile, vector_crs))
    
    #Vector type
    vector_type= type(vector)
    print("The type of file for {} is {} \n".format(polygonFile, vector_type))
    
    #Spatial extent
    vector_extent= vector.total_bounds
    print("The spatial extent for {} is {} \n".format(polygonFile, vector_extent))
    
    
def raster_change_crs(raster_original, raster_reprojected):
    '''
    The function takes two arguments:
    1) raster_original: Path to the raster of interest.
    2) raster_reprojected: Path where the new raster reprojected will be stored
    '''
    #Opening the raster using gdal
    raster_initial=gdal.Open(raster_original)   
    
    #Changing projection to the same of the shapefile "EPSG:4326"
    gdal.Warp(raster_reprojected, raster_initial, dstSRS="EPSG:4326")
    
    return raster_reprojected

def raster_stats(raster_reprojected, polygonFile, statsFile):
    '''
    "NOTE: FUNCTION UNDER DEVELOP""
    The function takes three arguments:
    1) raster_reprojected: Path to the raster of interest.
    2) polygonFile: path to the shapefile (vector) of interest (field boundary).
    3) statsFile: name and path of the text file in which the stats will be stored.
    '''
    
    #Loading individual bands
    with rasterio.open(raster_reprojected, mode="r") as raster:
        blue=raster.read(1)        
        
    with rasterio.open(raster_reprojected, mode="r") as raster:
        green=raster.read(2)
        
    with rasterio.open(raster_reprojected, mode="r") as raster:
        red=raster.read(3)        
        
    with rasterio.open(raster_reprojected, mode="r") as raster:
        nir=raster.read(4) 
                
    #Getting stats for each band    
    blue_stats = zonal_stats(polygonFile, blue,
                stats="count min mean max median")
    
    green_stats = zonal_stats(polygonFile, green,
                stats="count min mean max median")

    red_stats = zonal_stats(polygonFile, red,
                stats="count min mean max median")

    nir_stats = zonal_stats(polygonFile, nir, stats="count min mean max median")

    #Storing results to a textfile
    with open(statsFile, 'w') as outpu:        
            print("The stats for {} correspoding the the area of interest are the following:".format(raster_reprojected), file=outpu)
            print("Blue: {}".format(blue_stats), file=outpu)
            print("Green: {}".format(green_stats), file=outpu)
            print("Red: {}".format(red_stats), file=outpu)
            print("Nir: {}".format(nir_stats), file=outpu)
            outpu.close()    
    
def clip_fromRaster(raster_reprojected, polygonFile, rasterOut):
    '''
    The function takes three arguments:
    1) raster_reprojected= path to raster from which we are going to clip (MUST BE REPROJECTED to the same CRS as shapefile).
    2) polygonFile= path to shapefile/polygon that we want to clip.
    3) rasterOut= name and path to the new raster created from this clipping process.
    '''
    #Opening raster
    raster=gdal.Open(raster_reprojected)

    #Clipping
    gdal.Warp(rasterOut, raster, cutlineDSName=polygonFile,
                          cropToCutline=True, dstNodata=np.nan)
    
    return rasterOut

def calculate_NDVI(rasterIN, ndviOUT):
    '''
    The function takes two arguments:
    1) rasterIN: path to stacked layer that was clipped for the field of interest.
    2) ndviOUT= name and path to the new ndvi raster created.
    '''
   #Loading RED and NIR bands
    with rasterio.open(rasterIN) as raster:
        red=raster.read(3)        
        
    with rasterio.open(rasterIN) as raster:
        nir=raster.read(4)
    
    
    #Performing NDVI calculation
    ndvi = (nir.astype(float) - red.astype(float)) / (nir + red)
    
    
    #Saving NDVI image
    
    # Set spatial characteristics of the output object to mirror the input
    kwargs = raster.meta
    kwargs.update(
        dtype=rasterio.float32,
        count = 1)

    # Create the file
    with rasterio.open(ndviOUT, 'w', **kwargs) as dst:
        dst.write_band(1, ndvi.astype(rasterio.float32))
        dst.close()

    
def raster_hist(rasterToPlot, histOUT):
    '''
    The function takes two arguments:
    1) rasterToPlot: path to raster of interest.
    2) histOUT= name and path to the figure created.
    '''
    #Opening the raster of interest
    raster=rasterio.open(rasterToPlot)
    
    fig, axhist = plt.subplots(1, 1)
    show_hist(raster, bins=50, lw=0.0, stacked=False, alpha=0.3, 
              histtype='stepfilled', ax=axhist)
    axhist.set_xlabel('NDVI values')
    axhist.set_ylabel('Frequency')
    axhist.set_title('Histogram')
    axhist.get_legend().remove()
    plt.savefig(histOUT)


def raster_plot(rasterToPlot, plotOUT):
    '''
    The function takes two arguments:
    1) rasterToPlot: path to raster of interest.
    2) histOUT= name and path to the figure created.
    '''
   
    ndvi = rasterio.open(rasterToPlot, mode='r', rmasked=True)
    
    fig,ax = plt.subplots(1,1, figsize=(10,4)) #Platform for creating the figure
    show(ndvi, ax=ax, title="Spatial distribution of NDVI", cmap='RdYlGn')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    plt.savefig(plotOUT)
    plt.show()
        

def raster_hist_plot(rasterToPlot, figureOUT):
    '''
    The function takes two arguments:
    1) rasterToPlot: path to raster of interest.
    2) histOUT= name and path to the figure created.
    '''
   
    ndvi = rasterio.open(rasterToPlot, mode='r', rmasked=True)
    
    fig,(ax1, ax2) = plt.subplots(1,2, figsize=(10,4)) #Platform for creating the figure
    show(ndvi, ax=ax1, title="Spatial distribution of NDVI", cmap='RdYlGn')
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    show_hist(ndvi, bins=50, lw=0.0, stacked=False, alpha=0.3, 
              histtype='stepfilled', title="Histogram of NDVI", ax=ax2)
    ax2.set_xlabel('NDVI values')
    ax2.set_ylabel('Frequency')
    ax2.get_legend().remove()
    plt.savefig(figureOUT)
    plt.show()
    

        

# the following condition checks whether we are
# running as a script, in which case run the test code,
# or being imported, in which case don't.

if __name__ == '__main__':

    ''''Defining arguments for functions'''
    
    #Specific arguments
    folderPath= 'L1C_T16TDK_20160513T164811_20160513T213820/'
    rootname= 'S2A_OPER_MSI_L1C_TL_MTI__20160513T213820_A004656_T16TDK'
    polygonFile='Doug_cc_boundary/Doug_cc_boundary.shp'
    statsFile= folderPath+'bgrn_clip_doug_stats.txt'

    #Default
    stackBands= folderPath+'bgrn.tiff' 
    raster_reprojected= folderPath+'bgrn_reprojected.tiff'    
    rasterOut= folderPath+'bgrn_clip_doug.tiff'
    ndviOUT= folderPath+'bgrn_clip_doug_ndvi.tiff'
    histOUT= folderPath+'bgrn_clip_doug_ndvi_hist.png'    
    plotOUT= folderPath+'bgrn_clip_doug_ndvi_plot.png'    
    figureOUT= folderPath+'bgrn_clip_doug_ndvi_figure.png'    
    
    
    ''''Iniciating process'''

    #STEP 1. Stacking bands         
    stack_bands(folderPath, rootname, stackBands)
    
    #STEP 2. Information raster
    info_raster(stackBands)
    
    #STEP 3. Information vector
    info_vector(polygonFile)
    
    #STEP 4: Change coordinate reference system for raster        
    raster_change_crs(stackBands, raster_reprojected)
    
    # STEP PENDING (FUNCTION IS STILL UNDER DEVELOP): Get stats from raster.
    # raster_stats(raster_reprojected, polygonFile, statsFile)
    
    #STEP 5: Clipping field boundary   
    clip_fromRaster(raster_reprojected, polygonFile, rasterOut)
    
    #STEP 6. Calculate NDVI for clipped area
    calculate_NDVI(rasterOut, ndviOUT)
    
    #STEP 7. Create histogram and save to an image
    raster_hist(ndviOUT, histOUT)
    
    #STEP 8. Create plot and save to an image       
    raster_plot(ndviOUT, plotOUT)
    
    #STEP 9. Create plot and save to an image       
    raster_hist_plot(ndviOUT, figureOUT)
    
    
    
    