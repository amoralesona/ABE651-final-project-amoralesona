# ABE651-final-project-amoralesona
By: Ana Morales
Last time updated: Friday, May 7th, 2021

Script created for handling and processing inputs and outputs generated during vegetative indices calculation.

The program has the following functions for handling and processing raster and vector data:

1) stack_bands(folderPath, rootname, stackBands):
	The function takes three arguments:
	1) folderPath: path in which the individual layers are available
	2) rootname: each individual layer has a rootname that is common to all the bands
	3) stackBands: path and name of the output stacked layer

2) info_raster(rasterPath):
	The function takes one argument, which is the path to the raster of interest.

3) info_vector(polygonFile):
	The function takes one argument, which is the path to the shapefile (vector) of interest.
    
4) raster_change_crs(raster_original, raster_reprojected):
	The function takes two arguments:
	1) raster_original: Path to the raster of interest.
	2) raster_reprojected: Path where the new raster reprojected will be stored
   
5) clip_fromRaster(raster_reprojected, polygonFile, rasterOut):
	The function takes three arguments:
    1) raster_reprojected= path to raster from which we are going to clip (MUST BE REPROJECTED to the same CRS as shapefile).
    2) polygonFile= path to shapefile/polygon that we want to clip.
    3) rasterOut= name and path to the new raster created from this clipping process.

6) calculate_NDVI(rasterIN, ndviOUT):
	The function takes two arguments:
    1) rasterIN: path to stacked layer that was clipped for the field of interest.
    2) ndviOUT= name and path to the new ndvi raster created.    
    Source: https://developers.planet.com/planetschool/calculate-an-ndvi-in-python/

7) raster_min_max_mean_stdev(rasterIn):
    The function takes one argument, which is the raster of interest.
    The output are the min, max, mean, and standart deviation corresponding to the raster.
    
8) rasterToArray(rasterIN):
    The function takes a single raster in (e.g. NDVI map, elevation map, etc)
    and return a 1D array.

9) calculate_GNDVI(rasterIN, gndviOUT):
	The function takes two arguments:
	1) rasterIN: path to stacked layer that was clipped for the field of interest.
	2) gndviOUT= name and path to the new ndvi raster created.
	     
10) calculate_PPRB(rasterIN, pprbOUT):
    The function takes two arguments:
    1) rasterIN: path to stacked layer that was clipped for the field of interest.
    2) pprbOUT= name and path to the new ndvi raster created.

11) raster_individual_plot(rasterToPlot, plotOUT, figureTitle):
    The function takes two arguments:
    1) rasterToPlot: path to raster of interest.
    2) histOUT= name and path to the figure created.
    3) title= title for the plot (e.g. NDVI)    

12) raster_multiple_plot(ndviOUT, gndviOUT, pprbOUT, xmin, xmax, mulplotOUT):
    The function takes six arguments:
    1) ndviOUT: path to raster 1 of interest.
    2) gndviOUT= path to raster 2 of interest.
    3) pprbOUT= path to raster 3 of interest.
    4) xmin= min value for scale of color bar
    5) xmax= max value for scale of color bar
    4) mulplotOUT = name and path to the figure created.   

13) raster_multiple_vi(ndviOUT, gndviOUT, pprbOUT, xmin, xmax, mulplotOUT):
    This function focuses on creating a plot that compares three rasters in 
    three scenarios: BARE SOIL, VEGETATIVE PERIOD, and REPRODUCTIVE PERIOD.   
    The function takes five arguments:
    1) ndviOUT: path to raster 1 of interest.
    2) gndviOUT= path to raster 2 of interest.
    3) pprbOUT= path to raster 3 of interest.
    4) xmin= min value for scale of color bar
    5) xmax= max value for scale of color bar    
    6) mulplotOUT = name and path to the figure created.
