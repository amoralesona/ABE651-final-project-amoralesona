#!/usr/bin/env python
"""
ABE651 (Spring 2021)

Script created for handling and processing data for the project:
Assessment of Spatial Variability of Large-Scale Fields During the Crop Growing Season 
Based on Vegetative Indices Derived From Sentinel-2 Images

Author: Ana Morales
Last time updated: Friday, May 7, 2021

"""

# import required modules
import rasterio
import matplotlib.pyplot as plt
import geopandas as gpd
from osgeo import gdal
import numpy as np
from rasterio.plot import show
from matplotlib.ticker import FormatStrFormatter #For the 2 decimals in plots
import pandas as pd



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
    
    Source: https://developers.planet.com/planetschool/calculate-an-ndvi-in-python/
    '''
   #Loading RED and NIR bands
    with rasterio.open(rasterIN) as raster:
        red=raster.read(3)               
        nir=raster.read(4)
    
    # Allow division by zero
    np.seterr(divide='ignore', invalid='ignore')
    
    #Performing calculation
    ndvi = (nir.astype(float) - red.astype(float)) / (nir + red)
  
    
    # Set spatial characteristics of the output object to mirror the input
    kwargs = raster.meta
    kwargs.update(
        dtype=rasterio.float32,
        count = 1)

    # Create the file
    with rasterio.open(ndviOUT, 'w', **kwargs) as final_ndvi:
        final_ndvi.write_band(1, ndvi.astype(rasterio.float32))

    return final_ndvi

    
def raster_min_max_mean_stdev(rasterIn):
    '''
    The function takes one argument, which is the raster of interest.
    The output are the min, max, mean, and standart deviation corresponding to the raster.
    
    '''    
    
    #Open raster of interest and assing to a variable "raster"
    raster = gdal.Open(rasterIn)
    
    #Count number of bands in raster and assign to a variable "count"
    count = raster.RasterCount
    
    #Create an empty dictionary
    response = {}
    
    #Loop to get statistics per band in case raster has more than one
    for counter in range(1, count+1):
        stats = raster.GetRasterBand(counter).GetStatistics(0, 1)
        response["band_{}".format(counter)] = "%.3f, %.3f, %.3f, %.3f" % (stats[0], stats[1], stats[2], stats[3])
    print(response)  
    
 
def rasterToArray(rasterIN):
    '''
    The function takes a single raster in (e.g. NDVI map, elevation map, etc)
    and return a 1D array.
    '''    
    
    raster=gdal.Open(rasterIN) #Open raster
    rasterBand=raster.GetRasterBand(1) #Select band of interest, in this case default 1 because of single raster
    rasterArray=rasterBand.ReadAsArray() #Convert to numpy array
    rasterArrayFlat=rasterArray.flatten() #Convert to 1D array
    
    return rasterArrayFlat 


def calculate_GNDVI(rasterIN, gndviOUT):
    '''
    The function takes two arguments:
    1) rasterIN: path to stacked layer that was clipped for the field of interest.
    2) gndviOUT= name and path to the new ndvi raster created.
    '''
   #Loading NIR and GREEN bands
    with rasterio.open(rasterIN) as raster:
        nir=raster.read(4)               
        green=raster.read(2)
    
    # Allow division by zero
    np.seterr(divide='ignore', invalid='ignore')
    
    #Performing calculation
    gndvi = (nir.astype(float) - green.astype(float)) / (nir + green)
  
    
    # Set spatial characteristics of the output object to mirror the input
    kwargs = raster.meta
    kwargs.update(
        dtype=rasterio.float32,
        count = 1)

    # Create the file
    with rasterio.open(gndviOUT, 'w', **kwargs) as final_gndvi:
        final_gndvi.write_band(1, gndvi.astype(rasterio.float32))

    return final_gndvi
      

def calculate_PPRB(rasterIN, pprbOUT):
    '''
    The function takes two arguments:
    1) rasterIN: path to stacked layer that was clipped for the field of interest.
    2) pprbOUT= name and path to the new ndvi raster created.
    '''
   #Loading BLUE and GREEN bands
    with rasterio.open(rasterIN) as raster:
        blue=raster.read(1)               
        green=raster.read(2)
    
    # Allow division by zero
    np.seterr(divide='ignore', invalid='ignore')
    
    #Performing calculation
    pprb = (green.astype(float) - blue.astype(float)) / (green + blue)
  
    
    # Set spatial characteristics of the output object to mirror the input
    kwargs = raster.meta
    kwargs.update(
        dtype=rasterio.float32,
        count = 1)

    # Create the file
    with rasterio.open(pprbOUT, 'w', **kwargs) as final_pprb:
        final_pprb.write_band(1, pprb.astype(rasterio.float32))

    return final_pprb


def raster_individual_plot(rasterToPlot, plotOUT, figureTitle):
    '''
    The function takes two arguments:
    1) rasterToPlot: path to raster of interest.
    2) histOUT= name and path to the figure created.
    3) title= title for the plot (e.g. NDVI)
    ''' 
    
    raster = rasterio.open(rasterToPlot, mode='r', rmasked=True)
    
    fig,ax = plt.subplots(1,1, figsize=(10,4)) #Platform for creating the figure
    show(raster, 
            ax=ax, 
           cmap='RdYlGn', 
           vmin=-1, vmax=+1,
            title=figureTitle)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    plt.savefig(plotOUT)
    plt.show()  


def raster_multiple_plot(ndviOUT, gndviOUT, pprbOUT, xmin, xmax, mulplotOUT):
    '''
    The function takes six arguments:
    1) ndviOUT: path to raster 1 of interest.
    2) gndviOUT= path to raster 2 of interest.
    3) pprbOUT= path to raster 3 of interest.
    4) xmin= min value for scale of color bar
    5) xmax= max value for scale of color bar
    4) mulplotOUT = name and path to the figure created.
    ''' 
    
    ndvi = rasterio.open(ndviOUT, mode='r', rmasked=True)
    gndvi = rasterio.open(gndviOUT, mode='r', rmasked=True)
    pprb = rasterio.open(pprbOUT, mode='r', rmasked=True)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(24,8)) #Platform for creating the figure

    show(ndvi, ax=ax1, cmap='RdYlGn', vmin=xmin, vmax=xmax, title="NDVI")
    show(gndvi, ax=ax2, cmap='RdYlGn', vmin=xmin, vmax=xmax, title="GNDVI")
    show(pprb, ax=ax3, cmap='RdYlGn', vmin=xmin, vmax=xmax, title="PPRB")
    
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
          
    plt.savefig(mulplotOUT)  
    plt.show()    


def raster_multiple_vi(ndviOUT, gndviOUT, pprbOUT, xmin, xmax, mulplotOUT):
    '''
    This function focuses on creating a plot that compares three rasters in 
    three scenarios: BARE SOIL, VEGETATIVE PERIOD, and REPRODUCTIVE PERIOD.
    
    The function takes five arguments:
    1) ndviOUT: path to raster 1 of interest.
    2) gndviOUT= path to raster 2 of interest.
    3) pprbOUT= path to raster 3 of interest.
    4) xmin= min value for scale of color bar
    5) xmax= max value for scale of color bar    
    6) mulplotOUT = name and path to the figure created.
    ''' 
    
    ndvi = rasterio.open(ndviOUT, mode='r', rmasked=True)
    gndvi = rasterio.open(gndviOUT, mode='r', rmasked=True)
    pprb = rasterio.open(pprbOUT, mode='r', rmasked=True)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(24,8)) #Platform for creating the figure

    show(ndvi, ax=ax1, cmap='RdYlGn', vmin=xmin, vmax=xmax, title="BARE SOIL")
    show(gndvi, ax=ax2, cmap='RdYlGn', vmin=xmin, vmax=xmax, title="VEGETATIVE PERIOD")
    show(pprb, ax=ax3, cmap='RdYlGn', vmin=xmin, vmax=xmax, title="REPRODUCTIVE PERIOD")
    
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
          
    plt.savefig(mulplotOUT)  
    plt.show()   
        

# the following condition checks whether we are
# running as a script, in which case run the test code,
# or being imported, in which case don't.

if __name__ == '__main__':

    ''''DEFINING ARGUMENTS FOR FUNCTIONS'''
    
    #Specific arguments
    folderPath= 'Sentinel_Images/L1C_T16TEK_A016425_20180814T165018/IMG_DATA/'
    rootname= 'T16TEK_20180814T163901'
    polygonFile='Field_boundaries/Doug_cs_boundary/Doug_cs_boundary.shp'
    
    outputFolder= 'Outputs/'
    fieldFolder= outputFolder+"Doug/"
    
    imageDate= "D20180814_"
    
    #Default   
    stackBands= folderPath+'bgrn.tiff'     
    raster_reprojected= folderPath+'bgrn_reprojected.tiff'
    rasterOut= fieldFolder+imageDate+'bgrn_clip.tiff'
    
    ndviOUT= fieldFolder+imageDate+'_ndvi.tiff'
    gndviOUT= fieldFolder+imageDate+'_gndvi.tiff'
    pprbOUT= fieldFolder+imageDate+'_pprb.tiff' 
    
    mulplotOUT= fieldFolder+imageDate+'_VI_mulplot.png'
    mulBoxPlotOUT= fieldFolder+imageDate+'_VI_mulBoxPlot.png'
    
    ''''DATA PREPARATION AND QUALITY CHECK'''
    
    #STEP 1. Stacking bands         
    stack_bands(folderPath, rootname, stackBands)
    
    # STEP 2. Information raster
    info_raster(stackBands)
    
    # STEP 3. Information vector
    info_vector(polygonFile)
    
    #STEP 4: Change coordinate reference system for raster        
    raster_change_crs(stackBands, raster_reprojected)
       
    #STEP 5: Clipping field boundary   
    clip_fromRaster(raster_reprojected, polygonFile, rasterOut)
    
    #STEP 6. Calculate VI: NDVI for clipped area
    calculate_NDVI(rasterOut, ndviOUT) 

    #STEP 7. Calculate stats from NDVI for data quality check and clouds   
    raster_min_max_mean_stdev(ndviOUT)
    
    
    
    ''''VEGETATIVE INDICES CALCULATION FOR IMAGES WITH NO CLOUDS'''
        
    #STEP 8. Calculate VI: GNDVI for clipped area
    calculate_GNDVI(rasterOut, gndviOUT)

    #STEP 9. Calculate VI: PPRB for clipped area
    calculate_PPRB(rasterOut, pprbOUT) 



    ''''CALCULATION OF STATS FOR VEGETATIVE INDICESS'''

    #STEP 10. Calculate stats from GNDVI   
    raster_min_max_mean_stdev(gndviOUT)
    
    #STEP 11. Calculate stats from PPRB   
    raster_min_max_mean_stdev(pprbOUT)
     
     
     
    ''''GENERATION OF FIGURE WITH THREE VI (NDVI, GNDVI, PPRB)''' 

    #Plot vegetative indices
    raster_multiple_plot(ndviOUT, gndviOUT, pprbOUT, -1, +1, mulplotOUT)

    
    '''CREATION OF BOX-WHISKER PLOTS OF THREE VI (NDVI, GNDVI, PPRB) 
    CORRESPONDING TO THE SAME DATE'''
    
    #Define scenario for which image will be plotted
    Escenario = "Scenario 3: Reproductive period"

    #Create empthy dataframe 
    df= pd.DataFrame()

    #Creat 1D arrays for ndvi, gndvi, and pprb (using function created)    
    ndvi_array=rasterToArray(ndviOUT)
    gndvi_array=rasterToArray(gndviOUT)
    pprb_array=rasterToArray(pprbOUT)    
    
    #Add arrays to dataframe
    df["NDVI"]=ndvi_array.tolist()
    df["GNDVI"]=gndvi_array.tolist()
    df["PPRB"]=pprb_array.tolist()        

    #Plot boxplot       
    df.boxplot(column=['NDVI', 'GNDVI', 'PPRB'])
    title_boxplot = Escenario
    plt.title( title_boxplot )
    plt.savefig(mulBoxPlotOUT)
    plt.show()    
       
  
    '''PERSONALIZED PLOTS FOR SCENARIOS FOR COMPARISON OF VI (NDVI, GNDVI, PPRB) 
    FROM ****DIFFERENT DATES***'''
    
    #NDVI PLOT OF THREE SCENARIOS---------------------------------------------
    #Define arguments
    ndvi1= fieldFolder+"D20190401_ndvi.tiff"
    ndvi2= fieldFolder+"D20190712_ndvi.tiff"
    ndvi3= fieldFolder+"D20190829_ndvi.tiff"
    ndviplot= fieldFolder+"D2019_ndvi.png"
    
    #call function
    raster_multiple_vi(ndvi1, ndvi2, ndvi3, -1, +1, ndviplot)
    
    #GNDVI PLOT OF THREE SCENARIOS---------------------------------------------
    #Define arguments
    gndvi1= fieldFolder+"D20190401_gndvi.tiff"
    gndvi2= fieldFolder+"D20190712_gndvi.tiff"
    gndvi3= fieldFolder+"D20190829_gndvi.tiff"
    gndviplot= fieldFolder+"D2019_gndvi.png"
    #call function
    raster_multiple_vi(gndvi1, gndvi2, gndvi3, -1, +1, gndviplot)   
    
    #GNDVI PLOT OF THREE SCENARIOS---------------------------------------------
    #Define arguments
    pprb1= fieldFolder+"D20190401_pprb.tiff"
    pprb2= fieldFolder+"D20190712_pprb.tiff"
    pprb3= fieldFolder+"D20190829_pprb.tiff"
    pprbplot= fieldFolder+"D2019_pprb.png"
    #call function
    raster_multiple_vi(pprb1, pprb2, pprb3, -1, +1, pprbplot) 


    #PLOTS PER SCENARIO
    #Define arguments
    escenario1_plot=fieldFolder+"D2019_BareSoil_mulplot.png"
    escenario2_plot=fieldFolder+"D2019_VegP_mulplot.png"
    escenario3_plot=fieldFolder+"D2019_RepP_mulplot.png"
    #Call function
    raster_multiple_plot(ndvi1, gndvi1, pprb1, None, None, escenario1_plot)
    raster_multiple_plot(ndvi2, gndvi2, pprb2, None, None, escenario2_plot)
    raster_multiple_plot(ndvi3, gndvi3, pprb3, None, None, escenario3_plot)        


    
    

    

    
    