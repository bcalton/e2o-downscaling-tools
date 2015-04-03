# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 07:36:06 2014

@author: schelle
"""

from osgeo import gdal, gdalconst
import logging
import logging.handlers
import ConfigParser
import os
from numpy import *
import sys


def readMap(fileName, fileFormat,logger):
    """
    Read geographical file into memory
    """

    #Open file for binary-reading

    mapFormat = gdal.GetDriverByName(fileFormat)
    mapFormat.Register()
    ds = gdal.Open(fileName)
    if ds is None:
        logger.error('Could not open ' + fileName + '. Something went wrong!! Shutting down')
        sys.exit(1)
    # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX    = geotrans[1]
    resY    = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = linspace(originX+resX/2,originX+resX/2+resX*(cols-1),cols)
    y = linspace(originY+resY/2,originY+resY/2+resY*(rows-1),rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1) # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0,0,cols,rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    del ds
    return resX, resY, cols, rows, x, y, data, FillVal



def getmapname(number,prefix):
    """
    generate a pcraster type mapname based on timestep and prefix
    :var number: number of the mape
    :var prefix: prefix for the map

    :return: Name
    """
    print number
    below_thousand = number % 1000
    above_thousand = number / 1000
    mapname = str(prefix + '%0' + str(8-len(prefix)) + '.f.%03.f') % (above_thousand, below_thousand)

    return mapname

def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """ Write geographical data into file"""

    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(fileFormat)

    # Processing
    if verbose:
        print 'Writing to temporary file ' + fileName + '.tif'
    # Create Output filename from (FEWS) product name and date and open for writing
    TempDataset = driver1.Create(fileName + '.tif',data.shape[1],data.shape[0],1,gdal.GDT_Float32)
    # Give georeferences
    xul = x[0]-(x[1]-x[0])/2
    yul = y[0]+(y[0]-y[1])/2

    TempDataset.SetGeoTransform( [ xul, x[1]-x[0], 0, yul, 0, y[1]-y[0] ] )
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data,0,0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(FillVal)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to ' + fileName + '.map'
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print 'Removing temporary file ' + fileName + '.tif'
    os.remove(fileName + '.tif');

    if verbose:
        print 'Writing to ' + fileName + ' is done!'



def configget(log,config,section,var,default):
    """   
    Gets a string from a config file (.ini) and returns a default value if
    the key is not found. If the key is not found it also sets the value 
    with the default in the config-file
    
    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to get
        - default - default string
        
    Returns:
        - string - either the value from the config file or the default value
    """
    

    try:
        ret = config.get(section, var)
    except:
        ret = default
        log.info( "returning default (" + str(default) + ") for " + section + ":" + var)
        configset(config,section,var,str(default), overwrite=False)
    

    return ret       

def configset(config,section,var,value, overwrite=False):
    """   
    Sets a string in the in memory representation of the config object
    Deos NOT overwrite existing values if overwrite is set to False (default)
    
    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to set
        - value - the value to set
        - overwrite (optional, default is False)
   
    Returns:
        - nothing
        
    """
    
    if not config.has_section(section):
        config.add_section(section)
        config.set(section,var,value)
    else:     
        if not config.has_option(section,var):
            config.set(section,var,value)
        else:
            if overwrite:
                config.set(section,var,value)

def iniFileSetUp(configfile):
    """
    Reads .ini file and sets default values if not present
    """
    # TODO: clean up wflwo specific stuff
    #setTheEnv(runId='runId,caseName='caseName)
    # Try and read config file and set default options
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(configfile)
    return config

def setlogger(logfilename,loggername, level=logging.INFO):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    """

    try:
        #create logger
        logger = logging.getLogger(loggername)
        if not isinstance(level, int):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(level)
        ch = logging.FileHandler(logfilename,mode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        #create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        #add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        #add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger
    except IOError:
        print "ERROR: Failed to initialize logger with logfile: " + logfilename
        sys.exit(2)

def closeLogger(logger, ch):
    """
    Closes the logger
    """
    logger.removeHandler(ch)
    ch.flush()
    ch.close()
    return logger, ch



def resample_grid(gridZ_in,Xin,Yin,Xout,Yout,method='nearest'):
    """
    Resample a regular grid be supplyin original and new x, y row,col coordinates
    Missing values is set to 1E31, the PCRaster standard

    :param gridZ_in: datablock of original data (e.g. from readMap)
    :param Xin: X-coordinates of all columns (e.g. from readMap)
    :param Yin: Y-coordinates of all columns (e.g. from readMap)
    :param Xout: X-coordinates of all columns in the new map (e.g. from readMap)
    :param Yout: Y-coordinates of all columns in the new map (e.g. from readMap)
    :param method: linear, nearest, cubic, quintic
    :return: datablock of new grid.
    """
    
    from scipy import interpolate
    # we need to sort the y data (must be ascending)
    # and thus flip the image
    Yin.sort()
    # define interpolator
    if method in 'nearest linear':
        interobj = interpolate.RegularGridInterpolator((Yin,Xin), flipud(gridZ_in), method=method ,bounds_error=False)
        _x, _y = meshgrid(Xout, Yout)
        yx_outpoints = transpose([_y.flatten(), _x.flatten()])

        # interpolate
        res = interobj(yx_outpoints)
    else:
        interobj = interpolate.interp2d(Xin, Yin, gridZ_in,  kind=method, bounds_error=False)
        res = interobj(Xout, Yout)

    # make all combinations of X,Y so we can supply the grid as a list of x,y coordinates


    #reshape to the new grid
    result = reshape(res, (len(Yout), len(Xout)))
    result[isnan(result)] = 1E31
    
    return result


