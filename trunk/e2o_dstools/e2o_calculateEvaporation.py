"""
Determine downscaled reference evaporation from the eartH2Observe WFDEI forcing


usage:

    e2o_calculateEvaporation.py -I inifile

    -I inifile - ini file with settings which data to get

"""

import getopt, sys, os, netCDF4, glob
import osgeo.gdal as gdal
from osgeo.gdalconst import *
from osgeo import gdal, gdalconst
import datetime
from numpy import *
import numpy as np
from e2o_utils import *
import pdb
import pandas as pd
import pcraster as pcr
from scipy import interpolate
import scipy.ndimage
import shutil



def usage(*args):
    """
    Print usage information

    -  *args: command line arguments given
    """
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

class ncdatset():
    """
    Wrapper around the nc object to simplify things
    
    Opens the dataset and determines the number of dimensions
    3 = T, Lat Lon
    4 = T heigth, Lat Lon
    """
    def __init__(self,ncurl,logger):
        self.logger = logger
        try:
            self.nc = netCDF4.Dataset(ncurl)
        except:
            self.logger.error("Failed to open remote file: " + ncurl)
            sys.exit(12)

        self.lat = self.getlat(self.nc)
        if self.lat == None:
            self.logger.error("No lat information found!")
        self.lon = self.getlon(self.nc)
        if self.lon == None:
            self.logger.error("No lon information found!")

        self.heigth = self.getheigth(self.nc)
        self.dimensions = 3 if self.heigth == None else 4
        self.time =   self.gettime(self.nc)
        self.timesteps = self.time.shape[0]
        self.logger.debug(self.nc)


    def __del__(self):
        self.nc.close()
        
        
    def getlat(self,ncdataset):
        """
        """

        for a in ncdataset.variables:
            if  ncdataset.variables[a].standard_name == 'latitude':
                return ncdataset.variables[a]

        return None

    def getvarbyname(self,name):
        """
        """


        for a in self.nc.variables:
            if  self.nc.variables[a].standard_name == name:
                return self.nc.variables[a]

        return None

    def gettime(self,ncdataset):
        """
        """

        for a in ncdataset.variables:
            if  ncdataset.variables[a].standard_name == 'time':
                return ncdataset.variables[a]

        return None
        

    def getlon(self,ncdataset):
        """
        """

        for a in ncdataset.variables:
            if  ncdataset.variables[a].standard_name == 'longitude':
                return ncdataset.variables[a]

        return None


    def getheigth(self,ncdataset):
        """
        """

        for a in ncdataset.variables:
            if  ncdataset.variables[a].standard_name == 'height':
                return ncdataset.variables[a]

        return None



def get_times_daily(startdate,enddate, serverroot, wrrsetroot, filename,logger):
    """
    generate a dictionary with date/times and the NC files in which the data resides
    """

    numdays = enddate - startdate
    dateList = []
    filelist = {}
    for x in range (0, numdays.days + 1):
        dateList.append(startdate + datetime.timedelta(days = x))
    
    for thedate in dateList:
        ncfile = serverroot + wrrsetroot + "%d" % (thedate.year) + "/" + filename + "%d%02d.nc" % (thedate.year,thedate.month)
        filelist[str(thedate)] = ncfile
    
    return filelist, dateList


def get_times(startdate,enddate, serverroot, wrrsetroot, filename, timestepSeconds, logger):
    """
    generate a dictionary with date/times and the NC files in which the data resides for flexible timestep
    """

    numdays = enddate - startdate
    dateList = []
    filelist = {}
    #days
    for x in range (0, numdays.days + 1):
        #selected time-step in seconds
        delta = 0
        while delta < 86400:
            dateList.append(startdate + datetime.timedelta(seconds = delta))
            delta += timestepSeconds
    
    for thedate in dateList:
        ncfile = serverroot + wrrsetroot + "%d" % (thedate.year) + "/" + filename + "%d%02d.nc" % (thedate.year,thedate.month)
        filelist[str(thedate)] = ncfile
    
    return filelist, dateList


class getstepdaily():
    """
    class to get data from a set of NC files
    Initialise with a list of netcdf files and a variable name (standard_name)

    """
    
    def __init__(self,nclist,BB,varname,logger):
        """
        """
        self.o_nc_files = []
        self.list = nclist
        self.varname = varname
        self.BB = BB
        self.latidx = []
        self.lonidx =[]
        self.lat=[]
        self.lon=[]
        self.data = []
        self.logger = logger

    def getdate(self,thedate):

        datestr = str(thedate)
        lat = None
        lon = None
        window = None
        
        if datestr in self.list.keys():
            self.dset = ncdatset(self.list[datestr],self.logger)
            data = self.dset.getvarbyname(self.varname)
            lat = flipud(self.dset.lat[:])
            
            lon = self.dset.lon[:]
            (latidx,) = logical_and(lat >= self.BB['lat'][0], lat < self.BB['lat'][1]).nonzero()
            (lonidx,) = logical_and(lon >= self.BB['lon'][0], lon < self.BB['lon'][1]).nonzero()

            time = self.dset.time
            timeObj = netCDF4.num2date(time[:], units=time.units, calendar=time.calendar)
            
            dpos = thedate.day -1
            
            if self.dset.dimensions ==3:
                window = data[dpos,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[dpos,0,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            
            self.lat = lat[latidx]
            self.lon = lon[lonidx]
    
        else:
            self.logger.error( "cannot find: " + datestr)
            
        return self.lat, self.lon, window

    def getdates(self,alldates):
        """
        Does not work yet
        """
        lat = None
        lon = None
        ret = []
        
        lastnc = None
        
        # here we loop over nc files for speed reasons
        
        for theone in  unique(self.list.values()):
            self.dset = ncdatset(theone,self.logger)
            
            time = self.dset.time
            tar = time[:]
            timeObj = netCDF4.num2date(tar, units=time.units, calendar=time.calendar)
            #print timeObj
            spos = nonzero(timeObj == alldates[0])[0]
            if len(spos) != 1:
                spos = 0
            else:
                spos = int(spos)
 
            epos = nonzero(timeObj == alldates[-1])[0]
            if len(epos) != 1:
                epos = len(tar)
            else:
                epos = int(epos + 1)
            
            
            self.logger.info("Processing url: " + theone )
            
            data = self.dset.getvarbyname(self.varname)
            
            if data == None:
                self.logger.error("dataset with standard_name " + self.varname + " not found" )
                
            lat = self.dset.lat[:]
            lon = self.dset.lon[:]

            (self.latidx,) = logical_and(lat >= self.BB['lat'][0], lat <= self.BB['lat'][1]).nonzero()
            (self.lonidx,) = logical_and(lon >= self.BB['lon'][0], lon <= self.BB['lon'][1]).nonzero()


            if self.dset.dimensions ==3:
                window = data[spos:epos,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[spos:epos,0,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]

            
            self.lat = lat[self.latidx]
            self.lon = lon[self.lonidx]
            
            if len(ret) == 0:
                ret = window.copy()
            else:
                ret = vstack((ret,window))
        
        return ret
        
    def getdates_seconds(self,alldates):
        """
        Does not work yet
        """
        lat = None
        lon = None
        ret = []
        
        lastnc = None
        
        # here we loop over nc files fro speed reasons
        
        for theone in  unique(self.list.values()):
            self.dset = ncdatset(theone,self.logger)
            
            time = self.dset.time
            tar = time[:]
            timeObj = netCDF4.num2date(tar, units=time.units, calendar=time.calendar)
            #print timeObj
            spos = nonzero(timeObj == alldates[0])[0]
            if len(spos) != 1:
                spos = 0
            else:
                spos = int(spos)
 
            epos = nonzero(timeObj == alldates[-1])[0]
            if len(epos) != 1:
                epos = len(tar)
            else:
                epos = int(epos + 1)
            
            
            self.logger.info("Processing url: " + theone )
            
            data = self.dset.getvarbyname(self.varname)
            
            if data == None:
                self.logger.error("dataset with standard_name " + self.varname + " not found" )
                
            lat = self.dset.lat[:]
            lon = self.dset.lon[:]

            (self.latidx,) = logical_and(lat >= self.BB['lat'][0], lat <= self.BB['lat'][1]).nonzero()
            (self.lonidx,) = logical_and(lon >= self.BB['lon'][0], lon <= self.BB['lon'][1]).nonzero()


            if self.dset.dimensions ==3:
                window = data[spos:epos,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[spos:epos,0,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]

            
            self.lat = lat[self.latidx]
            self.lon = lon[self.lonidx]
            
            if len(ret) == 0:
                ret = window.copy()
            else:
                ret = vstack((ret,window))
        
        return ret
        
class getstep():
    """
    class to get data from a set of NC files for user defined timestep in seconds
    Initialise with a list of netcdf files and a variable name (standard_name)

    """
    
    def __init__(self,nclist,BB,varname,timestepSeconds,logger):
        """
        """
        self.o_nc_files = []
        self.list = nclist
        self.varname = varname
        self.BB = BB
        self.latidx = []
        self.lonidx =[]
        self.lat=[]
        self.lon=[]
        self.data = []
        self.logger = logger

    def getdate(self,thedate):

        datestr = str(thedate)
        lat = None
        lon = None
        window = None
        
        if datestr in self.list.keys():
            self.dset = ncdatset(self.list[datestr],self.logger)
            data = self.dset.getvarbyname(self.varname)
            lat = flipud(self.dset.lat[:])
            
            lon = self.dset.lon[:]
            (latidx,) = logical_and(lat >= self.BB['lat'][0], lat < self.BB['lat'][1]).nonzero()
            (lonidx,) = logical_and(lon >= self.BB['lon'][0], lon < self.BB['lon'][1]).nonzero()

            time = self.dset.time
            timeObj = netCDF4.num2date(time[:], units=time.units, calendar=time.calendar)
            
            if timestepSecond < 3600:
                dpos = thedate.second -1
            elif timestepSecond < 86400:
                dpos = thedate.hour -1
            else:
                dpos = thedate.day -1
            
            if self.dset.dimensions ==3:
                window = data[dpos,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[dpos,0,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            
            self.lat = lat[latidx]
            self.lon = lon[lonidx]
    
        else:
            self.logger.error( "cannot find: " + datestr)
            
        return self.lat, self.lon, window

    def getdates(self,alldates):
        """
        Does not work yet
        """
        lat = None
        lon = None
        ret = []
        
        lastnc = None
        
        # here we loop over nc files fro speed reasons
        
        for theone in  unique(self.list.values()):
            self.dset = ncdatset(theone,self.logger)
            
            time = self.dset.time
            tar = time[:]
            timeObj = netCDF4.num2date(tar, units=time.units, calendar=time.calendar)
            #print timeObj
            spos = nonzero(timeObj == alldates[0])[0]
            if len(spos) != 1:
                spos = 0
            else:
                spos = int(spos)
 
            epos = nonzero(timeObj == alldates[-1])[0]
            if len(epos) != 1:
                epos = len(tar)
            else:
                epos = int(epos + 1)
            
            
            self.logger.info("Processing url: " + theone )
            
            data = self.dset.getvarbyname(self.varname)
            
            if data == None:
                self.logger.error("dataset with standard_name " + self.varname + " not found" )
                
            lat = self.dset.lat[:]
            lon = self.dset.lon[:]

            (self.latidx,) = logical_and(lat >= self.BB['lat'][0], lat <= self.BB['lat'][1]).nonzero()
            (self.lonidx,) = logical_and(lon >= self.BB['lon'][0], lon <= self.BB['lon'][1]).nonzero()


            if self.dset.dimensions ==3:
                window = data[spos:epos,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[spos:epos,0,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]

            
            self.lat = lat[self.latidx]
            self.lon = lon[self.lonidx]
            
            if len(ret) == 0:
                ret = window.copy()
            else:
                ret = vstack((ret,window))
        
        return ret


def getmapname(number,prefix):
    """
    generate a pcraster type mapname based on timestep and prefix
    :var number: number of the mape
    :var prefix: prefix for the map

    :return: Name
    """
    below_thousand = number % 1000
    above_thousand = number / 1000
    mapname = str(prefix + '%0' + str(8-len(prefix)) + '.f.%03.f') % (above_thousand, below_thousand)

    return mapname

def save_as_mapsstack(lat,lon,data,times,directory,prefix="E2O",oformat="PCRaster"):        
    
    cnt = 0
    if not os.path.exists(directory):
        os.mkdir(directory)
    for a in times:
            mapname = getmapname(cnt,prefix)
            #print "saving map: " + os.path.join(directory,mapname)
            #writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[cnt,:,:]),-999.0)
            writeMap(os.path.join(directory,mapname),oformat,lon,lat,data[cnt,:,:],-999.0)
            cnt = cnt + 1    

def save_as_mapsstack_per_day(lat,lon,data,ncnt,directory,prefix="E2O",oformat="PCRaster",FillVal=1E31):
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    mapname = getmapname(ncnt,prefix)
    #print "saving map: " + os.path.join(directory,mapname)
    #writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[:,:]),-999.0)
    writeMap(os.path.join(directory,mapname),oformat,lon,lat,data[:,:],FillVal)

def save_as_gtiff(lat,lon,data,ncnt,directory,prefix,oformat='GTiff'):        
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    mapname = prefix + '.tif'
    #print "saving map: " + os.path.join(directory,mapname)
    writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[:,:]),-999.0)


def resample(highResdemname,prefix,ncnt,logger):
    """

    :param highResdemname:
    :param prefix:
    :param ncnt: map to resample
    :param logger:
    :return: resampled map
    """

    #create resample dir
    try:
        os.stat('resampled')
    except:
        os.mkdir('resampled')

    tif_mapname         = prefix+'.tif'
    pcraster_mapname    = getmapname(ncnt,prefix)

    tif_filename        = os.path.join('temp',tif_mapname)
    pcraster_filename   = os.path.join('temp',pcraster_mapname)
    pcraster_resFilename   = os.path.join('resampled',pcraster_mapname)

    command= 'gdal_translate -of %s %s %s' % ('GTiff',pcraster_filename,tif_filename)
    os.system(command)

    # Source
    src_filename    = os.path.join('temp',tif_mapname)
    src             = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
    src_proj        = src.GetProjection()
    src_geotrans    = src.GetGeoTransform()

    # We want a section of source that matches this:
    match_filename  = highResdemname
    match_ds        = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
    match_proj      = match_ds.GetProjection()
    match_geotrans  = match_ds.GetGeoTransform()
    wide            = match_ds.RasterXSize
    high            = match_ds.RasterYSize

    # Output / destination
    dst_filename = os.path.join('resampled',tif_mapname)
    dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
    dst.SetGeoTransform( match_geotrans )
    dst.SetProjection( match_proj)

    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)

    del dst # Flush

    resX, resY, cols, rows, x, y, data, FillVal = readMap(dst_filename,'GTiff',logger)

    # nodig?? data    = np.flipud(dataUD)

    return data

def resampleDEM(nameHighResDEM, nameLowResDEM,logger):
    """

    :param nameHighResDEM:
    :param nameLowResDEM:
    :param logger:
    :return: elevationcorrection, hiresdem, upscaled_lowresdem
    """
    
    #create temp dir
    try:
        os.stat('temp')
    except:
        os.mkdir('temp')
    
    # Source
    src_filename    = nameLowResDEM
    src             = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
    src_proj        = src.GetProjection()
    src_geotrans    = src.GetGeoTransform()

    # We want a section of source that matches this:
    match_filename  = nameHighResDEM
    match_ds        = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
    match_proj      = match_ds.GetProjection()
    match_geotrans  = match_ds.GetGeoTransform()
    wide            = match_ds.RasterXSize
    high            = match_ds.RasterYSize

    # Output / destination
    dst_filename = os.path.join('temp','DEM.tif')
    dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
    dst.SetGeoTransform( match_geotrans )
    dst.SetProjection( match_proj)

    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)

    del dst # Flush

    resX, resY, cols, rows, XI, YI, resLowResDEM, FillVal = readMap(dst_filename,'GTiff',logger)
    resX, resY, cols, rows, XI, YI, highResDEM, FillVal = readMap(match_filename,'GTiff',logger)

    highResDEM  = np.maximum(0,highResDEM)
    elevationCorrection = highResDEM - resLowResDEM
       
    return elevationCorrection, highResDEM, resLowResDEM
    
def resample(highResdemname,prefix,ncnt,logger):
    """

    :param highResdemname:
    :param prefix:
    :param ncnt: map to resample
    :param logger:
    :return: resampled map
    """
    
    #create resample dir
    try:
        os.stat('resampled')
    except:
        os.mkdir('resampled')    
    
    tif_mapname         = prefix+'.tif'
    pcraster_mapname    = getmapname(ncnt,prefix)
    
    tif_filename        = os.path.join('temp',tif_mapname)
    pcraster_filename   = os.path.join('temp',pcraster_mapname)
    pcraster_resFilename   = os.path.join('resampled',pcraster_mapname)
    
    command= 'gdal_translate -of %s %s %s' % ('GTiff',pcraster_filename,tif_filename)
    os.system(command)
      
    # Source
    src_filename    = os.path.join('temp',tif_mapname)
    src             = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
    src_proj        = src.GetProjection()
    src_geotrans    = src.GetGeoTransform()

    # We want a section of source that matches this:
    match_filename  = highResdemname
    match_ds        = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
    match_proj      = match_ds.GetProjection()
    match_geotrans  = match_ds.GetGeoTransform()
    wide            = match_ds.RasterXSize
    high            = match_ds.RasterYSize

    # Output / destination
    dst_filename = os.path.join('resampled',tif_mapname)
    dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
    dst.SetGeoTransform( match_geotrans )
    dst.SetProjection( match_proj)

    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)

    del dst # Flush
       
    resX, resY, cols, rows, x, y, data, FillVal = readMap(dst_filename,'GTiff',logger)

    # nodig?? data    = np.flipud(dataUD)
    
    return data

def correctTemp(Temp,elevationCorrection):

    lapse_rate = 0.006 # [ K m-1 ]


    #apply elevation correction
    Temp_cor   = Temp - lapse_rate * elevationCorrection
    
    return Temp_cor 

def correctRsin(Rsin,currentdate,radiationCorDir,logger):
    """
    Corrects incoming radiation using the information from the e2o_radiation module

    :param Rsin:
    :param currentdate:
    :param radiationCorDir:
    :param logger:
    :return:  corrected incoming radiation
    """
    #get day of year
    logger.info("Correcting incoming radiation with DEM...")
    tt  = currentdate.timetuple()
    JULDAY = tt.tm_yday
    #read data from radiation correction files
    mapname     = getmapname(JULDAY,'FLAT')
    fileName    = os.path.join(radiationCorDir,mapname)
    resX, resY, cols, rows, x, y, data, FillVal          = readMap(fileName,'PCRaster',logger)
    
    resX, resY, cols, rows, x, y, flat, FillVal            = readMap((os.path.join(radiationCorDir,(getmapname(JULDAY,'FLAT')))),'PCRaster',logger)   
    resX, resY, cols, rows, x, y, flatdir, FillVal         = readMap((os.path.join(radiationCorDir,(getmapname(JULDAY,'FLATDIR')))),'PCRaster',logger)
    resX, resY, cols, rows, x, y, cor, FillVal             = readMap((os.path.join(radiationCorDir,(getmapname(JULDAY,'COR')))),'PCRaster',logger)
    resX, resY, cols, rows, x, y, cordir, FillVal          = readMap((os.path.join(radiationCorDir,(getmapname(JULDAY,'CORDIR')))),'PCRaster',logger)
    #ratio direct - diffuse
    missmask = cor == FillVal
    flat[flat == 0.0] = 0.00001
    flatdir[flatdir == 0.0] = 0.00001
    ratio           = flatdir / flat
    Rsin_dir        = ratio * Rsin
    #corrected Rsin direct for elevation and slope
    Rsin_dir_cor    = (cordir/flatdir)*Rsin_dir
    Rsin_cor        = Rsin_dir_cor + (Rsin - Rsin_dir)
    Rsin_cor[missmask] = FillVal
    Rsin_cor[Rsin_cor < 0.0] = FillVal
    
    return Rsin_cor

def correctPres(relevantDataFields, Pressure, highResDEM, resLowResDEM,FillVal=1E31):
    """
    Correction of air pressure for DEM based altitude correction based on barometric formula

    :param relevantDataFields:
    :param Pressure:
    :param highResDEM:
    :param resLowResDEM:
    :return: corrected pressure

    relevantDataFields : ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','SurfaceIncidentShortwaveRadiation','NearSurfaceWindSpeed']
    """
    

    Tmean   =  relevantDataFields[0]
              


    g            = 9.81         # gravitational constant [m s-2]
    R_air        = 8.3144621    # specific gas constant for dry air [J mol-1 K-1]
    Mo           = 0.0289644    # molecular weight of gas [g / mol]
    lapse_rate   = 0.006        # lapse rate [K m-1]

    # Why is this, switched off for now...
    #highResDEM  = np.maximum(0,highResDEM)

    Pres_cor = zeros_like(Tmean)
    Pres_cor    = Pressure *( (Tmean/ ( Tmean + lapse_rate * (highResDEM - resLowResDEM))) ** (g * Mo / (R_air * lapse_rate)))

    Pres_cor[isnan(Pres_cor)] = FillVal

    return Pres_cor

def PenmanMonteith(lat, currentdate, relevantDataFields, Tmax, Tmin):
    """

    :param lat:
    :param currentdate:
    :param relevantDataFields:
    :param Tmax:
    :param Tmin:
    :return:

    relevantDataFields : ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','SurfaceIncidentShortwaveRadiation','NearSurfaceWindSpeed']
    """

    Tmean   =  relevantDataFields[0]
    Rlin    =  relevantDataFields[1]
    Pres    =  relevantDataFields[2]
    Q       =  relevantDataFields[3]
    Rsin    =  relevantDataFields[4]
    Wsp     =  relevantDataFields[5]
              
    """
    Computes Penman-Monteith reference evaporation
    Inputs:
        Rsin:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) incoming shortwave radiation [W m-2]
        Rlin:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) incoming longwave radiation [W m-2]
        Tmean:       netCDF obj or NumPy array   -- 3D array (time, lat, lon) daily mean temp [K]
        Tmax:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) daily max. temp [K]
        Tmin:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) daily min. temp [K]
        Wsp:         netCDF obj or NumPy array   -- 3D array (time, lat, lon) wind speed [m s-2]
        Q:           netCDF obj or NumPy array   -- 3D array (time, lat, lon) specific humididy [kg kg-1]
        Pres:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) Surface air pressure [Pa]
        Pet:         netCDF obj or NumPy array   -- 3D array (time, lat, lon) for target data
    Outputs:
        trg_var:    netCDF obj or NumPy array   -- 3D array (time, lat, lon) for target data, updated with computed values
    """
    cp           = 1013         # specific heat of air 1013 [J kg-1 K-1]
    TimeStepSecs = 86400        # timestep in seconds
    karman       = 0.41         # von Karman constant [-]
    vegh         = 0.12         # vegetation height [m] 
    refh         = 2            # reference height where wind speed is measured
    alpha        = 0.23         # albedo, 0.23 [-]
    rs           = 70           # surface resistance, 70 [s m-1]
    R            = 287.058      # Universal gas constant [J kg-1 K-1]
    convmm       = 1000*TimeStepSecs # conversion from meters to millimeters
    sigma        = 4.903e-9     # stephan boltzmann [W m-2 K-4]
    eps          = 0.622        # ratio of water vapour/dry air molecular weights [-]
    g            = 9.81         # gravitational constant [m s-2]
    R_air        = 8.3144621    # specific gas constant for dry air [J mol-1 K-1]
    Mo           = 0.0289644    # molecular weight of gas [g / mol]
    lapse_rate   = 0.006        # lapse rate [K m-1]
    
    #CALCULATE EXTRATERRESTRIAL RADIATION
    #get day of year
    tt  = currentdate.timetuple()
    JULDAY = tt.tm_yday
#    #Latitude radians
    LatRad= lat*np.pi/180.0
    test = np.tan(LatRad)
#    ### water euivalent extraterrestial radiation ###    
#    # declination (rad)
    declin = 0.4093*(np.sin(((2.0*pi*JULDAY)/365.0)-1.405))
#    # sunset hour angle
    arccosInput = (-(np.tan(LatRad))*(np.tan(declin)))
#    
    arccosInput = np.minimum(1,arccosInput)
    arccosInput = np.maximum(-1,arccosInput)
    sunangle = np.arccos(arccosInput)
#    # distance of earth to sun
    distsun = 1+0.033*(np.cos((2*pi*JULDAY)/365.0))
    # Ra = water equivalent extra terrestiral radiation in MJ day-1
    Ra = ((24 * 60 * 0.082) / 3.14) * distsun * (sunangle*(np.sin(LatRad))*(np.sin(declin))+(np.cos(LatRad))*(np.cos(declin))*(np.sin(sunangle)))
    
    #CALCULATE ACTUAL VAPOR PRESSURE
    # saturation vapour pressure [Pa]
    es = lambda T:610.8*np.exp((17.27*(Tmean-273.15))/((Tmean-273.15)+237.3))
    es_min  = es(Tmin)
    es_max  = es(Tmax)
    es_mean = (es_min+es_max)/2.

    # actual vapour pressure
    ea = lambda Pres, Q, eps: -(Q*Pres)/((eps-1)*Q-eps)
    ea_mean = ea(Pres, Q, eps)
    ea_mean_kPa = ea_mean / 1000
    
    #clear sky solar radiation MJ d-1
    Rso = np.maximum(0.1,((0.75+(2*0.00005)) * Ra))
    
    Rsin_MJ = 0.086400 * Rsin   
    
    Rlnet_MJ = - sigma * ((Tmax**4+Tmin**4)/2) * (0.34 - 0.14 * np.sqrt(np.maximum(0,(ea_mean_kPa)))) * (1.35*np.minimum(1,(Rsin_MJ / Rso))-0.35)
    
    Rlnet_Watt = Rlnet_MJ / 0.086400
    
    Rnet  = np.maximum(0,((1-alpha)*Rsin + Rlnet_Watt))
    
    # vapour pressure deficit
    vpd = np.maximum(es_mean - ea_mean, 0.)
    
    # density of air [kg m-3]
    rho = Pres/(Tmean*R)
    
    # Latent heat [J kg-1]
    Lheat = (2.501-(0.002361*(Tmean-273.15)))*1e6

    # slope of vapour pressure [Pa K-1]
    deltop  = 4098. *(610.8*np.exp((17.27*(Tmean-273.15))/((Tmean-273.15)+237.3)))
    delbase = ((Tmean-273.15)+237.3)**2
    delta   = deltop/delbase

    # psychrometric constant
    gamma   = cp*Pres/(eps*Lheat)
    
    # aerodynamic resistance
    z = 10 # height of wind speed variable (10 meters above surface)
    Wsp_2 = Wsp*4.87/(np.log(67.8*z-5.42))
    ra = 208./Wsp_2
    
    PETtop  = np.maximum((delta*Rnet + rho*cp*vpd/ra),1)
    PETbase = np.maximum((delta + gamma*(1+rs/ra)),1)
    PET     = np.maximum(PETtop/PETbase, 0)
    PETmm   = np.maximum((PET/Lheat*TimeStepSecs),0)

    if PET.any() == float("inf"):
        sys.exit("Value infinity found")
    else:
        pass
    
    return PETmm
    
def PriestleyTaylor(lat, currentdate, relevantDataFields, Tmax, Tmin):
    
    """
    relevantDataFields : ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','SurfaceIncidentShortwaveRadiation','NearSurfaceWindSpeed']
    """
    
    Tmean   =  relevantDataFields[0]
    Rlin    =  relevantDataFields[1]
    Pres    =  relevantDataFields[2]
    Q       =  relevantDataFields[3]
    Rsin    =  relevantDataFields[4]             

    alpha        = 0.23         # albedo, 0.23 [-]
    sigma        = 4.903e-9     # stephan boltzmann [MJ K-4 m-2 day-1]
    cp_pt        = 0.001013     # specific heat of air 1013 [MJ kg-1 C-1]
    a            = 1.26         # Priestley-Taylor coefficient [-]
    eps          = 0.622        # ratio of molecular weight of water to dry air [-]
    g            = 9.81         # gravitational constant [m s-2]
    R_air        = 8.3144621    # specific gas constant for dry air [J mol-1 K-1]
    Mo           = 0.0289644    # molecular weight of gas [g / mol]
    lapse_rate   = 0.006        # lapse rate [K m-1]
    
    """ http://agsys.cra-cin.it/tools/evapotranspiration/help/Priestley-Taylor.html """
    
    #CALCULATE EXTRATERRESTRIAL RADIATION
    #get day of year
    tt  = currentdate.timetuple()
    JULDAY = tt.tm_yday
#    #Latitude radians
    LatRad= lat*np.pi/180.0
    test = np.tan(LatRad)
#    ### water euivalent extraterrestial radiation ###    
#    # declination (rad)
    declin = 0.4093*(np.sin(((2.0*pi*JULDAY)/365.0)-1.405))
#    # sunset hour angle
    arccosInput = (-(np.tan(LatRad))*(np.tan(declin)))
#    
    arccosInput = np.minimum(1,arccosInput)
    arccosInput = np.maximum(-1,arccosInput)
    sunangle = np.arccos(arccosInput)
#    # distance of earth to sun
    distsun = 1+0.033*(np.cos((2*pi*JULDAY)/365.0))
    # Ra = water equivalent extra terrestiral radiation in MJ day-1
    Ra = ((24 * 60 * 0.082) / 3.14) * distsun * (sunangle*(np.sin(LatRad))*(np.sin(declin))+(np.cos(LatRad))*(np.cos(declin))*(np.sin(sunangle)))
    
    #CALCULATE ACTUAL VAPOR PRESSURE
    # saturation vapour pressure [Pa]
    es = lambda T:610.8*np.exp((17.27*(Tmean-273.15))/((Tmean-273.15)+237.3))
    es_min  = es(Tmin)
    es_max  = es(Tmax)
    es_mean = (es_min+es_max)/2.

    # actual vapour pressure
    ea = lambda Pres, Q, eps: -(Q*Pres)/((eps-1)*Q-eps)
    ea_mean = ea(Pres, Q, eps)
    ea_mean_kPa = ea_mean / 1000
        
    #clear sky solar radiation MJ d-1
    Rso = np.maximum(0.1,((0.75+(2*0.00005)) * Ra))
    
    Rsin_MJ = 0.086400 * Rsin 
    
    Rlnet_MJ = - sigma * ((Tmax**4+Tmin**4)/2) * (0.34 - 0.14 * np.sqrt(np.maximum(0,(ea_mean_kPa)))) * (1.35*np.minimum(1,(Rsin_MJ / Rso))-0.35)
    
    Rlnet_Watt = Rlnet_MJ / 0.086400
        
    preskPa     = Pres / 1000 
    latentHeat  = 2.501 - ( 0.002361 * ( Tmean - 273.15 ) ) # latent heat of vaporization (MJ kg-1)

    slope_exp   = (17.27*(Tmean - 273.15)) / ((Tmean - 273.15) + 237.3)
    slope_div   = ((Tmean - 273.15) + 237.3)**2
    slope       = 4098 * (0.6108 * (np.exp(slope_exp))) / slope_div
    
    psychConst  = cp_pt * ( preskPa ) / (latentHeat * eps ) # psychrometric constant (kPa degreesC-1)
          
    # net radiation [MJ m-2]
    Rnet  = np.maximum(0,((1-alpha) *Rsin_MJ + Rlnet_MJ))
    
    PETmm = a * (1 / latentHeat) *( (slope * Rnet ) /  ( slope + psychConst ) )
    
    return PETmm
    
def hargreaves(lat, currentdate, relevantDataFields, Tmax, Tmin):
    
    Tmean = relevantDataFields[0]

    #get day of year
    tt  = currentdate.timetuple()
    JULDAY = tt.tm_yday
#    #Latitude radians
    LatRad= lat*np.pi/180.0
    test = np.tan(LatRad)
#    ### water euivalent extraterrestial radiation ###    
#    # declination (rad)
    declin = 0.4093*(np.sin(((2.0*pi*JULDAY)/365.0)-1.405))
#    # sunset hour angle
    arccosInput = (-(np.tan(LatRad))*(np.tan(declin)))
#    
    arccosInput = np.minimum(1,arccosInput)
    arccosInput = np.maximum(-1,arccosInput)
    sunangle = np.arccos(arccosInput)
#    # distance of earth to sun
    distsun = 1+0.033*(np.cos((2*pi*JULDAY)/365.0))
#    # SO = water equivalent extra terrestiral radiation in mm/day
    Ra = 15.392*distsun*(sunangle*(np.sin(LatRad))*(np.sin(declin))+(np.cos(LatRad))*(np.cos(declin))*(np.sin(sunangle)))
    strDay       = str(JULDAY)

    airT = relevantDataFields[0]
    PETmm = 0.0023*Ra*((np.maximum(0,(Tmean-273.0))) + 17.8)*sqrt(np.maximum(0,(Tmax-Tmin)))
 
    return Tmean

#### MAIN ####

def main(argv=None):

    # Set all sorts of defaults.....
    serverroot = "http://wci.earth2observe.eu/thredds/dodsC/"
    wrrsetroot = "ecmwf/met_forcing_v0/"
    
    #available variables with corresponding file names and standard_names as in NC files
    variables = ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','Rainfall','SurfaceIncidentShortwaveRadiation','SnowfallRate','NearSurfaceWindSpeed']
    filenames = ["Tair_daily_E2OBS_","LWdown_daily_E2OBS_","PSurf_daily_E2OBS_","Qair_daily_E2OBS_",\
                    "Rainf_daily_E2OBS_","SWdown_daily_E2OBS_","Snowf_daily_E2OBS_","Wind_daily_E2OBS_"]
    standard_names = ['air_temperature','surface_downwelling_longwave_flux_in_air','surface_air_pressure','specific_humidity',\
                        'rainfal_flux','surface_downwelling_shortwave_flux_in_air','snowfall_flux','wind_speed']
    prefixes = ["Tmean","LWdown","PSurf","Qair",\
                    "Rainf","SWdown","Snowf","Wind"]
    
    #tempdir

    #defaults in absence of ini file
    filename = "Tair_daily_E2OBS_"
    standard_name ='air_temperature'
    startyear = 1980
    endyear= 1980
    startmonth = 1
    endmonth = 12
    latmin = 51.25
    latmax = 51.75
    lonmin = 5.25
    lonmax = 5.75
    startday = 1
    endday = 1
    getDataForVar = True
    calculateEvap = False
    evapMethod = None
    downscaling = None
    nrcalls = 0
    loglevel=logging.INFO

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            exit()   
    try:
        opts, args = getopt.getopt(argv, 'I:l:')
    except getopt.error, msg:
        usage(msg)
    
    for o, a in opts:
        if o == '-I': inifile = a
    
    logger = setlogger("e2o_calculateEvap.log","e2o_calculateEvaporation",level=loglevel)
    #logger, ch = setlogger("e2o_getvar.log","e2o_getvar",level=loglevel)
    logger.info("Reading settings from ini: " + inifile)
    theconf = iniFileSetUp(a)
    
    # Read period from file
    startyear = int(configget(logger,theconf,"selection","startyear",str(startyear)))
    endyear = int(configget(logger,theconf,"selection","endyear",str(endyear)))
    endmonth = int(configget(logger,theconf,"selection","endmonth",str(endmonth)))
    startmonth = int(configget(logger,theconf,"selection","startmonth",str(startmonth)))
    endday = int(configget(logger,theconf,"selection","endday",str(endday)))
    startday = int(configget(logger,theconf,"selection","startday",str(startday)))
    start = datetime.datetime(startyear,startmonth,startday)
    end = datetime.datetime(endyear,endmonth,endday)
        
    #read remaining settings from in file        
    lonmax = float(configget(logger,theconf,"selection","lonmax",str(lonmax)))
    lonmin = float(configget(logger,theconf,"selection","lonmin",str(lonmin)))
    latmax = float(configget(logger,theconf,"selection","latmax",str(latmax)))
    latmin = float(configget(logger,theconf,"selection","latmin",str(latmin)))
    BB = dict(
           lon=[ lonmin, lonmax],
           lat= [ latmin, latmax]
           )
    serverroot = configget(logger,theconf,"url","serverroot",serverroot)
    wrrsetroot = configget(logger,theconf,"url","wrrsetroot",wrrsetroot)
    oformat = configget(logger,theconf,"output","format","PCRaster")
    odir = configget(logger,theconf,"output","directory","output/")
    oprefix = configget(logger,theconf,"output","prefix","E2O")
    radcordir = configget(logger,theconf,"downscaling","radiationcordir","output_rad")
    FNhighResDEM = configget(logger,theconf,"downscaling","highResDEM","downscaledem.map")
    FNlowResDEM = configget(logger,theconf,"downscaling","lowResDEM","origdem.map")
    saveAllData = int(configget(logger,theconf,"output","saveall","0"))

    # Check whether downscaling should be applied
    downscaling   = configget(logger,theconf,"selection","downscaling",downscaling)
    if downscaling == 'True':
        #get grid info
        resX, resY, cols, rows, highResLon, highResLat, highResDEM, FillVal = readMap(FNhighResDEM,'GTiff',logger)
        LresX, LresY, Lcols, Lrows, lowResLon, lowResLat, lowResDEM, FillVal = readMap(FNlowResDEM,'GTiff',logger)
        #writeMap("DM.MAP","PCRaster",highResLon,highResLat,highResDEM,FillVal)
        #elevationCorrection, highResDEM, resLowResDEM = resampleDEM(FNhighResDEM,FNlowResDEM,logger)
        demmask=highResDEM != FillVal
        mismask=highResDEM == FillVal
        Ldemmask=lowResDEM != FillVal
        Lmismask=lowResDEM == FillVal
        # Fille gaps in high res DEM with Zeros for ineterpolation purposes
        lowResDEM[Lmismask] = 0.0
        resLowResDEM = resample_grid(lowResDEM,lowResLon, lowResLat,highResLon, highResLat,method='linear',FillVal=0.0)
        lowResDEM[Lmismask] = FillVal
        elevationCorrection = highResDEM - resLowResDEM

        #writeMap("1.MAP","PCRaster",highResLon,highResLat,resLowResDEM,FillVal)
        #writeMap("2.MAP","PCRaster",highResLon,highResLat,elevationCorrection,FillVal)
        #writeMap("3.MAP","PCRaster",highResLon,highResLat,highResDEM,FillVal)


    #Check whether evaporation should be calculated
    calculateEvap   = configget(logger,theconf,"selection","calculateEvap",calculateEvap)

    if calculateEvap == 'True':
        evapMethod      = configget(logger,theconf,"selection","evapMethod",evapMethod)
                
    if evapMethod == 'PenmanMonteith':
        relevantVars = ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','SurfaceIncidentShortwaveRadiation','NearSurfaceWindSpeed']        
    elif evapMethod == 'Hargreaves':
        relevantVars = ['Temperature']
    elif evapMethod == 'PriestleyTaylor':
        relevantVars = ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','SurfaceIncidentShortwaveRadiation']
        
    currentdate = start
    ncnt = 0
    while currentdate <= end:
        ncnt += 1
        if ncnt > 0:
            # Get all daily datafields needed and aad to list
            relevantDataFields = []
            for i in range (0,len(variables)):
                if variables[i] in relevantVars:
                    mapname = os.path.join(odir,getmapname(ncnt,oprefix))
                    if os.path.exists(mapname):
                        logger.info("Skipping map: " + mapname)
                    else:
                        logger.info("Getting data field: " + filename)
                        filename = filenames[i]
                        standard_name = standard_names[i]
                        logger.info("Get file list..")
                        tlist, timelist = get_times_daily(currentdate,currentdate,serverroot, wrrsetroot, filename,logger)
                        logger.info("Get dates..")
    
                        ncstepobj = getstepdaily(tlist,BB,standard_name,logger)
    
                        logger.info("Get data...: " + str(timelist))
                        mstack = ncstepobj.getdates(timelist)
                        mean_as_map = mstack.mean(axis=0)
                        logger.info("Get data body...")
                        if downscaling == 'True':
                            logger.info("Downscaling..." + variables[i])
                            #save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,mean_as_map,int(ncnt),'temp',prefixes[i],oformat='GTiff')
                            #mean_as_map = resample(FNhighResDEM,prefixes[i],int(ncnt),logger)
                            mean_as_map = resample_grid(mean_as_map,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method='linear',FillVal=FillVal)
                            mean_as_map = flipud(mean_as_map)
                            mean_as_map[mismask] = FillVal
                            if variables[i]     == 'Temperature':
                                mean_as_map     = correctTemp(mean_as_map,elevationCorrection)
                            if variables[i]     == 'SurfaceIncidentShortwaveRadiation':
                                mean_as_map     = correctRsin(mean_as_map,currentdate,radcordir,logger)
                            if variables[i]     == 'SurfaceAtmosphericPressure':
                                mean_as_map     = correctPres(relevantDataFields, mean_as_map, highResDEM, resLowResDEM,FillVal=FillVal)
                            mean_as_map[mismask] = FillVal

                        relevantDataFields.append(mean_as_map)
                                            
                        #only needed once
                        if nrcalls ==0:
                            nrcalls = nrcalls + 1
                            latitude = ncstepobj.lat[:]
                            #assuming a resolution of 0.5 degrees
                            LATITUDE = np.ones(((2*(latmax-latmin)),(2*(lonmax-lonmin))))
                            for i in range (0,int((2*(lonmax-lonmin)))):
                                LATITUDE[:,i]=LATITUDE[:,i]*latitude
                            if downscaling == 'True':    
                                #save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,LATITUDE,int(ncnt),'temp','lat',oformat=oformat)
                                #LATITUDE = resample(FNhighResDEM,'lat',int(ncnt),logger)
                                LATITUDE = zeros_like(highResDEM)
                                for i in range(0,LATITUDE.shape[1]):
                                    LATITUDE[:,i] = highResLat


                            #assign longitudes and lattitudes grids
                            if downscaling == 'True':  
                                lons = highResLon
                                lats = highResLat
                            else:
                                lons = ncstepobj.lon
                                lats = ncstepobj.lat
    
            if evapMethod == 'PenmanMonteith':
    
                mapname = os.path.join(odir,getmapname(ncnt,oprefix))
                if os.path.exists(mapname):
                    logger.info("Skipping map: " + mapname)
                else:
                    # retrieve 3 hourly Temperature and calculate max and min Temperature
                    filename = 'Tair_E2OBS_'
                    standard_name = 'air_temperature'
                    timestepSeconds = 10800
                
                    tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, filename,timestepSeconds,logger )
                    ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)
                    mstack = ncstepobj.getdates(timelist)
                    tmin = flipud(mstack.min(axis=0))
                    tmax = flipud(mstack.max(axis=0))
                    if downscaling == 'True':
                        tmin = resample_grid(tmin,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method='linear',FillVal=FillVal)
                        tmax = resample_grid(tmax,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method='linear',FillVal=FillVal)
                        tmin = correctTemp(tmin,elevationCorrection)
                        tmax = correctTemp(tmax,elevationCorrection)
                        tmax[mismask] = FillVal
                        tmin[mismask] = FillVal


                    PETmm = PenmanMonteith(LATITUDE, currentdate, relevantDataFields, tmax, tmin)
                    # FIll out unrealistic values
                    PETmm[mismask] = FillVal
                    PETmm[PETmm < -10.0] = FillVal
                    PETmm[PETmm > 135.0] = FillVal

                    logger.info("Saving PM PET data for: " +str(currentdate))
                    save_as_mapsstack_per_day(lats,lons,PETmm,int(ncnt),odir,prefix=oprefix,oformat=oformat,FillVal=FillVal)
                    if saveAllData:
                        save_as_mapsstack_per_day(lats,lons,tmin,int(ncnt),odir,prefix='TMIN',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,tmax,int(ncnt),odir,prefix='TMAX',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[1],int(ncnt),odir,prefix='RLIN',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[2],int(ncnt),odir,prefix='PRESS',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[3],int(ncnt),odir,prefix='REL',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[4],int(ncnt),odir,prefix='RSIN',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[5],int(ncnt),odir,prefix='WIN',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[0],int(ncnt),odir,prefix='TEMP',oformat=oformat,FillVal=FillVal)


                  
            if evapMethod == 'PriestleyTaylor':
                mapname = os.path.join(odir,getmapname(ncnt,oprefix))
                if os.path.exists(mapname):
                    logger.info("Skipping map: " + mapname)
                else:
                    # retrieve 3 hourly Temperature and calculate min Temperature
                    filename = 'Tair_E2OBS_'
                    standard_name = 'air_temperature'
                    timestepSeconds = 10800
                
                    tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, filename,timestepSeconds,logger )
                    ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)
                    mstack = ncstepobj.getdates(timelist)
                    tmin = mstack.min(axis=0)
                    tmax = mstack.max(axis=0)
                    if downscaling == 'True':
                        tmin = resample_grid(tmin,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method='linear',FillVal=FillVal)
                        tmax = resample_grid(tmax,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method='linear',FillVal=FillVal)
                        tmin = correctTemp(tmin,elevationCorrection)
                        tmax = correctTemp(tmax,elevationCorrection)
                        tmax[mismask] = FillVal
                        tmin[mismask] = FillVal

                    PETmm = PriestleyTaylor(LATITUDE, currentdate, relevantDataFields, tmax, tmin)
                    PETmm[mismask] = FillVal
                    PETmm[PETmm < -10.0] = FillVal
                    PETmm[PETmm > 135.0] = FillVal
                    logger.info("Saving PriestleyTaylor PET data for: " +str(currentdate))

                    save_as_mapsstack_per_day(lats,lons,PETmm,int(ncnt),odir,prefix=oprefix,oformat=oformat,FillVal=FillVal)
                    if saveAllData:
                        save_as_mapsstack_per_day(lats,lons,tmin,int(ncnt),odir,prefix='TMIN',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,tmax,int(ncnt),odir,prefix='TMAX',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[4],int(ncnt),odir,prefix='RSIN',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[0],int(ncnt),odir,prefix='TEMP',oformat=oformat,FillVal=FillVal)


            if evapMethod == 'Hargreaves':
                mapname = os.path.join(odir,getmapname(ncnt,oprefix))
                if os.path.exists(mapname):
                    logger.info("Skipping map: " + mapname)
                else:
                    # retrieve 3 hourly Temperature and calculate max and min Temperature
                    filename = 'Tair_E2OBS_'
                    standard_name = 'air_temperature'
                    timestepSeconds = 10800
    
                    logger.info("Get times 3 hr data..")
                    tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, filename,timestepSeconds,logger )
                    logger.info("Get actual 3hr data...")
                    ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)
                    mstack = ncstepobj.getdates(timelist)
                    tmin = mstack.min(axis=0)
                    tmax = mstack.max(axis=0)
                    if downscaling == 'True':
                        tmin = resample_grid(tmin,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method='linear',FillVal=FillVal)
                        tmax = resample_grid(tmax,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method='linear',FillVal=FillVal)
                        tmin = correctTemp(tmin,elevationCorrection)
                        tmax = correctTemp(tmax,elevationCorrection)
                        tmax[mismask] = FillVal
                        tmin[mismask] = FillVal
                        
                    logger.info("Start hargreaves..")
                    PETmm = hargreaves(LATITUDE,currentdate,relevantDataFields, tmax, tmin)
                    PETmm[mismask] = FillVal
                    PETmm[PETmm < -10.0] = FillVal
                    PETmm[PETmm > 135.0] = FillVal
    
                    logger.info("Saving Hargreaves PET data for: " +str(currentdate))
                    save_as_mapsstack_per_day(lats,lons,PETmm,int(ncnt),odir,prefix=oprefix,oformat=oformat,FillVal=FillVal)
                    if saveAllData:
                        save_as_mapsstack_per_day(lats,lons,tmin,int(ncnt),odir,prefix='TMIN',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,tmax,int(ncnt),odir,prefix='TMAX',oformat=oformat,FillVal=FillVal)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[0],int(ncnt),odir,prefix='TEMP',oformat=oformat,FillVal=FillVal)

        else:
            pass
        
        currentdate += datetime.timedelta(days=1)       
    
    logger.info("Done.")

if __name__ == "__main__":
    main()



