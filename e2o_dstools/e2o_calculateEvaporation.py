import getopt, sys, os, glob
import osgeo.gdal as gdal
from osgeo.gdalconst import *
from osgeo import gdal, gdalconst
import datetime
from numpy import *
import numpy as np
from e2o_dstools.e2o_utils import *
import pdb
import pandas as pd
import pcraster as pcr
from scipy import interpolate
import scipy.ndimage
import shutil
import e2o_dstools

"""
Determine downscaled reference evaporation from the eartH2Observe WFDEI forcing

usage:

    e2o_calculateEvaporation.py -I inifile [-S start][-E end][-l loglevel]

    -I inifile - ini file with settings which data to get


    -S start - start timestep (default is 1)
    -E end - end timestep default is last one defined in ini file (from date)
    -l loglevel - DEBUG, INFO, WARN, ERROR
"""



def usage(*args):
    """
    Print usage information

    -  *args: command line arguments given
    """
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)


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

def fixmapwithmask(data,mask):
    """
    Fill maps value within mask with average
    :param data:
    :param mask:
    :return:
    """

    #avgfill = data[isfinite(data)].flatten().mean()

    data = fill_data(data)
    data[~mask] = nan

    return data

def save_as_mapsstack_per_day(lat,lon,data,ncnt,directory,prefix="E2O",oformat="PCRaster",FillVal=1E31,gzip=False,miscount=3827895):
    import platform

    if not os.path.exists(directory):
        os.mkdir(directory)
    mapname = getmapname(ncnt,prefix)
    #print "saving map: " + os.path.join(directory,mapname)
    #writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[:,:]),-999.0)

    mymiscount = sum(isnan(data))
    if not miscount ==  mymiscount:
        print mapname + ": oops expected nr missing data: "+ str(miscount) + " got "+ str(mymiscount)

    data[isnan(data)] = FillVal
    writeMap(os.path.join(directory,mapname),oformat,lon,lat,data[:,:],FillVal)
    if gzip:
        if 'Linux' in platform.system():
            os.system('gzip ' + os.path.join(directory,mapname))




def save_as_gtiff(lat,lon,data,ncnt,directory,prefix,oformat='GTiff'):

    if not os.path.exists(directory):
        os.mkdir(directory)
    mapname = prefix + '.tif'
    #print "saving map: " + os.path.join(directory,mapname)
    writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[:,:]),-999.0)

def convertCoordinateSystem(sourcetif,desttif,destmap): 

    #AZ - UTM38N:
    #command= 'gdalwarp -s_srs EPSG:4326 -t_srs EPSG:32638 -te 147000 4151000 986000 4775000 -tr 1000 1000 %s %s' % (sourcetif,desttif)
    #AZ - lat-lon:
    command= 'gdalwarp -s_srs EPSG:4326 -t_srs EPSG:4326 -te 41.0 37.5 50.5 43.0 -tr 0.01 0.01 %s %s' % (sourcetif,desttif)
    os.system(command)     

    command= 'gdal_translate -of PCRaster -ot Float32 %s %s' % (desttif,destmap)
    os.system(command)   

def correctRsin(Rsin,currentdate,radiationCorDir,logger,missmask):
    """
    Corrects incoming radiation using the information from the e2o_radiation module

    :param Rsin:
    :param currentdate:
    :param radiationCorDir:
    :param logger:
    :return:  corrected incoming radiation
    """

    #AtmPcorOrg = power(((288.0-0.0065*AltitudeOrg)/288.0),5.256)
    #AtmPcorDownscale = power(((288.0-0.0065*AltitudeDownscale)/288.0),5.256)
    #AtmCorFac = AtmPcorOrg/AtmPcorDownscale
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
    resX, resY, cols, rows, x, y, optcoradjust, FillVal          = readMap((os.path.join(radiationCorDir,(getmapname(JULDAY,'OPT')))),'PCRaster',logger)

    #ratio direct - diffuse

    flat[flat == 0.0] = 0.00001
    flatdir[flatdir == 0.0] = 0.00001
    #ratio           = flatdir / flat
    # Determine clear sky factor
    Rsin = Rsin * optcoradjust
    Kc = minimum(1.0,Rsin/maximum(0.0001,flat))
    Rsin_dir        = Kc * Rsin
    #corrected Rsin direct for elevation and slope
    Rsin_dir_cor    = (cordir/maximum(0.0001,flatdir))*Rsin_dir
    Rsin_cor        = Rsin_dir_cor + (Rsin - Rsin_dir)
    Rsin_cor        = Rsin_dir_cor + (Rsin - Rsin_dir)
    Rsin_cor[missmask] = nan
    Rsin_cor[Rsin_cor < 0.0] = nan

    return Rsin_cor, Kc


def correctRsinAlt(Rsin,AltitudeOrg,AltitudeDownscale,logger):
    """
    Corrects incoming radiation for DEM only

    :param Rsin:
    :param AltitudeOrg:
    :param AltitudeDownscale:
    :param logger:
    :return radcor:
    """

    #AtmPcorOrg = power(((288.0-0.0065*AltitudeOrg)/288.0),5.256)
    #AtmPcorDownscale = power(((288.0-0.0065*AltitudeDownscale)/288.0),5.256)
    #AtmCorFac = AtmPcorOrg/AtmPcorDownscale
    Rsin_cor = Rsin# * AtmCorFac

    return Rsin_cor


def correctPres(relevantDataFields, Pressure, LapseRate, highResDEM, resLowResDEM,FillVal=1E31):
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

    
    array_thousands         = 1000 * ones_like(highResDEM)
    array_minumThousands    = -1000 * ones_like(highResDEM)
    
    elevationDiff = np.minimum((highResDEM - resLowResDEM),array_thousands)
    elevationDiff = np.maximum(elevationDiff,array_minumThousands)

    Pres_cor = zeros_like(Tmean)
    Pres_cor    = Pressure *( (Tmean/ ( Tmean + LapseRate * elevationDiff)) ** (g * Mo / (R_air * LapseRate)))

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
    alpha        = 0.23         # albedo, 0.23 [-]
    rs           = 70           # surface resistance, 70 [s m-1]
    R            = 287.058      # Universal gas constant [J kg-1 K-1]
    convmm       = 1000*TimeStepSecs # conversion from meters to millimeters
    sigma        = 4.903e-9     # stephan boltzmann [W m-2 K-4]
    eps          = 0.622        # ratio of water vapour/dry air molecular weights [-]
    g            = 9.81         # gravitational constant [m s-2]
    R_air        = 8.3144621    # specific gas constant for dry air [J mol-1 K-1]
    Mo           = 0.0289644    # molecular weight of gas [g / mol]


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

    return PETmm

#### MAIN ####


def main(argv=None):

    # Set all sorts of defaults.....
    serverroot = "http://wci.earth2observe.eu/thredds/dodsC/"
    wrrsetroot = "ecmwf/met_forcing_v0/"

    #available variables with corresponding file names and standard_names as in NC files
    variables = ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','Rainfall','SurfaceIncidentShortwaveRadiation','SnowfallRate','NearSurfaceWindSpeed','LapseRate']
    standard_names = ['air_temperature','surface_downwelling_longwave_flux_in_air','surface_air_pressure','specific_humidity',\
                        'rainfal_flux','surface_downwelling_shortwave_flux_in_air','snowfall_flux','wind_speed','air_temperature_lapse_rate']
    prefixes = ["Tmean","LWdown","PSurf","Qair",\
                    "Rainf","SWdown","Snowf","Wind","lapseR"]

    #tempdir

    #defaults in absence of ini file
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
    FNlowResDEM = 'lowresdem\demWRR1.tif'
    getDataForVar = True
    calculateEvap = False
    evapMethod = None
    downscaling = None
    resampling = None
    useVarLapseRate = None
    StartStep = 1
    EndStep = 0

    nrcalls = 0
    loglevel=logging.INFO

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            exit()
    try:
        opts, args = getopt.getopt(argv, 'I:l:S:E:')
    except getopt.error, msg:
        usage(msg)

    for o, a in opts:
        if o == '-I': inifile = a
        if o == '-E': EndStep = int(a)
        if o == '-S': StartStep = int(a)
        if o == '-l': exec "loglevel = logging." + a

    logger = setlogger("e2o_calculateEvap.log","e2o_calculateEvaporation",level=loglevel)
    #logger, ch = setlogger("e2o_getvar.log","e2o_getvar",level=loglevel)
    logger.info("Reading settings from ini: " + inifile)
    theconf = iniFileSetUp(inifile)


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
    saveAllData = int(configget(logger,theconf,"output","saveall","0"))
    netcdfout = configget(logger,theconf,"output","netcdfout","None")

    # Set low resolution DEM filename for either WRR1 or WRR2:
    if 'met_forcing_v1' in wrrsetroot:
        FNlowResDEM     = e2o_dstools.get_data('DEM-WRR2.tif')
    else:
        FNlowResDEM     = e2o_dstools.get_data('DEM-WRR1.tif')
    # Get filenames
    if 'met_forcing_v1' in wrrsetroot: 
        filenames = \
            ["Tair_daily_EI_025_","LWdown_daily_EI_025_","PSurf_daily_EI_025_","Qair_daily_EI_025_",\
                    "Rainf_daily_EI_025_","SWdown_daily_EI_025_","Snowf_daily_EI_025_","Wind_daily_EI_025_","lapseM_EI_025_"]
        temperatureFile = ["Tair_EI_025_"]
    else:
        filenames = ["Tair_daily_E2OBS_","LWdown_daily_E2OBS_","PSurf_daily_E2OBS_","Qair_daily_E2OBS_",\
                    "Rainf_daily_E2OBS_","SWdown_daily_E2OBS_","Snowf_daily_E2OBS_","Wind_daily_E2OBS_"]
        temperatureFile = ["Tair_E2OBS_"]        
                    
    # Check whether downscaling should be applied
    resamplingtype   = configget(logger,theconf,"downscaling","resamplingtype","linear")
    downscaling   = configget(logger,theconf,"downscaling","downscaling",downscaling)
    resampling   = configget(logger,theconf,"downscaling","resampling",resampling)
    useVarLapseRate   = configget(logger,theconf,"downscaling","useVarLapseRate",useVarLapseRate)

    if downscaling == 'True' or resampling =="True":
        FNhighResDEM = configget(logger,theconf,"downscaling","highResDEM","downscaledem.map")
        # get grid info
        resX, resY, cols, rows, highResLon, highResLat, highResDEM, FillVal = readMap(FNhighResDEM,'GTiff',logger)
        LresX, LresY, Lcols, Lrows, lowResLon, lowResLat, lowResDEM, FillVal = readMap(FNlowResDEM,'GTiff',logger)
        X = highResLon
        Y = highResLat
        #writeMap("DM.MAP","PCRaster",highResLon,highResLat,highResDEM,FillVal)
        #elevationCorrection, highResDEM, resLowResDEM = resampleDEM(FNhighResDEM,FNlowResDEM,logger)
        demmask=isfinite(highResDEM)
        miscount = sum(isnan(highResDEM))
        mismask=isnan(highResDEM)
        Ldemmask=isfinite(lowResDEM)

        Lmismask=isnan(lowResDEM)
        # Fille gaps in high res DEM with Zeros for ineterpolation purposes
        lowResDEM[Lmismask] = 0.0
        resLowResDEM = resample_grid(lowResDEM,lowResLon, lowResLat,highResLon, highResLat,method=resamplingtype,FillVal=0.0)

        lowResDEM[Lmismask] = nan
        elevationCorrection = highResDEM - resLowResDEM

        #set Bounding box used to minimize data retrieval from nc file
        lonmax = math.ceil(highResLon[-1]+1)
        lonmin = math.floor(highResLon[0]-1)
        latmax = math.ceil(highResLat[0]+1)
        latmin = math.floor(highResLat[-1]-1)
        BB = dict(
           lon=[ lonmin, lonmax],
           lat= [ latmin, latmax]
           )
    else:
        FillVal     = nan
        LresX, LresY, Lcols, Lrows, lowResLon, lowResLat, lowResDEM, FillVal = readMap(FNlowResDEM,'GTiff',logger)
        mismask     = lowResDEM == nan
        X = lowResLon
        Y = lowResLat

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
    if EndStep == 0:
        EndStep = (end - start).days + 1

    logger.info("Will start at step: " + str(StartStep) + " date/time: " + str(start + datetime.timedelta(days=StartStep)))
    logger.info("Will stop at step: " + str(EndStep) + " date/time: " + str(start + datetime.timedelta(days=EndStep)))

    if netcdfout != 'None':
        ncout = netcdfoutput(netcdfout,X,Y,logging,start,EndStep - StartStep)
    while currentdate <= end:
        ncnt += 1

        if ncnt > 0 and ncnt >= StartStep and ncnt <= EndStep:
            # Get all daily datafields needed and aad to list
            relevantDataFields = []
            # Get all data for this timestep
            mapname = os.path.join(odir,getmapname(ncnt,oprefix))
            if os.path.exists(mapname) or os.path.exists(mapname + ".gz") or os.path.exists(mapname + ".zip"):
                logger.info("Skipping map: " + mapname)
            else:
                #retrieve lapse_rate if needed
                if useVarLapseRate == 'True' and 'met_forcing_v1' in wrrsetroot:            
                    if downscaling == 'True' or resampling =="True":
                        filename        = filenames[-1]
                        standard_name   = standard_names[-1]
                        tlist, timelist = get_times_daily(currentdate,currentdate,serverroot, wrrsetroot, filename,logger)
                        ncstepobj       = getstepdaily(tlist,BB,standard_name,logger)
                        mstack          = ncstepobj.getdates(timelist)
                        lapse_rate      = flipud(mstack.mean(axis=0))
                        lapse_rate      = resample_grid(lapse_rate,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method=resamplingtype,FillVal=nan)
                else:
                    lapse_rate = -0.006       
                #retrieve all other variables                                                                              
                for i in range (0,len(variables)):
                    if variables[i] in relevantVars:
                        filename = filenames[i]
                        logger.info("Getting data field: " + filename)
                        standard_name = standard_names[i]
                        logger.info("Get file list..")
                        tlist, timelist = get_times_daily(currentdate,currentdate,serverroot, wrrsetroot, filename,logger)
                        logger.info("Get dates..")

                        ncstepobj = getstepdaily(tlist,BB,standard_name,logger)

                        logger.info("Get data...: " + str(timelist))
                        mstack = ncstepobj.getdates(timelist)
                        mean_as_map = flipud(mstack.mean(axis=0))

                        logger.info("Get data body...")
                        if downscaling == 'True' or resampling =="True":
                            logger.info("Downscaling..." + variables[i])
                            #save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,mean_as_map,int(ncnt),'temp',prefixes[i],oformat='GTiff')
                            #mean_as_map = resample(FNhighResDEM,prefixes[i],int(ncnt),logger)
                            mean_as_map = resample_grid(mean_as_map,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method=resamplingtype,FillVal=FillVal)
                            #mean_as_map = flipud(mean_as_map)
                            mean_as_map[mismask] = nan
                            if downscaling == "True":
                                if variables[i]     == 'Temperature':
                                    mean_as_map     = correctTemp(mean_as_map, lapse_rate, elevationCorrection)

                                if variables[i]     == 'SurfaceIncidentShortwaveRadiation':
                                    if radcordir != 'None':
                                        mean_as_map, Kc    = correctRsin(mean_as_map,currentdate,radcordir, logger,mismask)
                                    else:
                                        mean_as_map = correctRsinAlt(mean_as_map,lowResDEM,highResDEM,  logger)

                                if variables[i]     == 'SurfaceAtmosphericPressure':
                                    mean_as_map     = correctPres(relevantDataFields, mean_as_map, lapse_rate, highResDEM, resLowResDEM,FillVal=nan)
                            mean_as_map[mismask] = nan

                        relevantDataFields.append(mean_as_map)

                        #only needed once
                        if nrcalls ==0:
                            nrcalls = nrcalls + 1
                            latitude = ncstepobj.lat[:]
                            longitude   = ncstepobj.lon[:]
                            if downscaling == 'True' or resampling == "True":
                                LATITUDE = zeros_like(highResDEM)
                                for i in range(0,LATITUDE.shape[1]):
                                    LATITUDE[:,i] = highResLat
                            else:
                                LATITUDE = np.zeros((len(latitude),len(longitude)),float)
                                for i in range(0,LATITUDE.shape[1]):
                                    LATITUDE[:,i] = latitude

                            #assign longitudes and lattitudes grids
                            if downscaling == 'True' or resampling == "True":
                                lons = highResLon
                                lats = highResLat
                            else:
                                lons = ncstepobj.lon
                                lats = ncstepobj.lat

                if evapMethod == 'PenmanMonteith':
                    # retrieve 3 hourly Temperature and calculate max and min Temperature
                    standard_name = 'air_temperature'
                    timestepSeconds = 10800
                    temperatureFilename = temperatureFile[0]
                    tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, temperatureFilename,timestepSeconds,logger )
                    ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)
                    mstack = ncstepobj.getdates(timelist)
                    tmin = flipud(mstack.min(axis=0))
                    tmax = flipud(mstack.max(axis=0))
                    if downscaling == 'True' or resampling == "True":
                        tmin = resample_grid(tmin,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method=resamplingtype,FillVal=nan)
                        tmax = resample_grid(tmax,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method=resamplingtype,FillVal=nan)
                        if downscaling == "True":
                            tmin = correctTemp(tmin,lapse_rate,elevationCorrection)
                            tmax = correctTemp(tmax,lapse_rate,elevationCorrection)
                        tmax[mismask] = nan
                        tmin[mismask] = nan


                    PETmm = PenmanMonteith(LATITUDE, currentdate, relevantDataFields, tmax, tmin)
                    # FIll out unrealistic valuea
                    PETmm[mismask] = nan
                    PETmm[PETmm < -10.0] = nan
                    PETmm[PETmm > 135.0] = nan

                    PETmm = fixmapwithmask(PETmm,demmask)
                    tmin = fixmapwithmask(tmin, demmask)
                    tmax = fixmapwithmask(tmax, demmask)
                    for a in [0,1,2,3,4,5]:
                        relevantDataFields[a] = fixmapwithmask(relevantDataFields[a], demmask)


                    logger.info("Saving PM PET data for: " +str(currentdate))

                    if netcdfout != 'None':
                        ncout.savetimestep(ncnt,PETmm,unit="mm",var="PotEvap",name='water_potential_evaporation_amount')
                    else:
                        save_as_mapsstack_per_day(lats, lons, PETmm, int(ncnt), odir, prefix=oprefix, oformat=oformat,
                                                  FillVal=FillVal,miscount=miscount)

                    if saveAllData:
                        if netcdfout != 'None':
                            ncout.savetimestep(ncnt, relevantDataFields[0], unit="degrees", var=relevantVars[0],
                                               name='air_temperature')
                            ncout.savetimestep(ncnt, relevantDataFields[1], unit="W/m^2", var=relevantVars[1],
                                               name='surface_air_pressure')
                        else:
                            save_as_mapsstack_per_day(lats,lons,tmin,int(ncnt),odir,prefix='TMIN',oformat=oformat,FillVal=FillVal,miscount=miscount)
                            save_as_mapsstack_per_day(lats,lons,tmax,int(ncnt),odir,prefix='TMAX',oformat=oformat,FillVal=FillVal,miscount=miscount)
                            save_as_mapsstack_per_day(lats,lons,relevantDataFields[1],int(ncnt),odir,prefix='RLIN',oformat=oformat,FillVal=FillVal,miscount=miscount)
                            save_as_mapsstack_per_day(lats,lons,relevantDataFields[2],int(ncnt),odir,prefix='PRESS',oformat=oformat,FillVal=FillVal,miscount=miscount)
                            save_as_mapsstack_per_day(lats,lons,relevantDataFields[3],int(ncnt),odir,prefix='REL',oformat=oformat,FillVal=FillVal,miscount=miscount)
                            save_as_mapsstack_per_day(lats,lons,relevantDataFields[4],int(ncnt),odir,prefix='RSIN',oformat=oformat,FillVal=FillVal,miscount=miscount)
                            save_as_mapsstack_per_day(lats,lons,relevantDataFields[5],int(ncnt),odir,prefix='WIN',oformat=oformat,FillVal=FillVal,miscount=miscount)
                            save_as_mapsstack_per_day(lats,lons,relevantDataFields[0],int(ncnt),odir,prefix='TEMP',oformat=oformat,FillVal=FillVal,miscount=miscount)
                            #save_as_mapsstack_per_day(lats,lons,Kc,int(ncnt),odir,prefix='KC',oformat=oformat,FillVal=FillVal)


                if evapMethod == 'PriestleyTaylor':
                    # retrieve 3 hourly Temperature and calculate min Temperature
                    temperatureFilename = temperatureFile[0]
                    standard_name = 'air_temperature'
                    timestepSeconds = 10800

                    tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, temperatureFilename,timestepSeconds,logger )
                    ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)
                    mstack = ncstepobj.getdates(timelist)
                    tmin = flipud(mstack.min(axis=0))
                    tmax = flipud(mstack.max(axis=0))
                    if downscaling == 'True' or resampling == "True":
                        tmin = resample_grid(tmin,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method=resamplingtype,FillVal=nan)
                        tmax = resample_grid(tmax,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method=resamplingtype,FillVal=nan)
                        if downscaling == "True":
                            tmin = correctTemp(tmin,lapse_rate,elevationCorrection)
                            tmax = correctTemp(tmax,lapse_rate,elevationCorrection)
                        tmax[mismask] = nan
                        tmin[mismask] = nan

                    PETmm = PriestleyTaylor(LATITUDE, currentdate, relevantDataFields, tmax, tmin)
                    PETmm[mismask] = nan
                    PETmm[PETmm < -10.0] = nan
                    PETmm[PETmm > 135.0] = nan

                    PETmm = fixmapwithmask(PETmm,demmask)
                    tmin = fixmapwithmask(tmin, demmask)
                    tmax = fixmapwithmask(tmax, demmask)
                    for a in [0,4]:
                        relevantDataFields[a] = fixmapwithmask(relevantDataFields[a], demmask)

                    logger.info("Saving PriestleyTaylor PET data for: " +str(currentdate))

                    save_as_mapsstack_per_day(lats,lons,PETmm,int(ncnt),odir,prefix=oprefix,oformat=oformat,FillVal=FillVal,miscount=miscount)
                    if saveAllData:
                        save_as_mapsstack_per_day(lats,lons,tmin,int(ncnt),odir,prefix='TMIN',oformat=oformat,FillVal=FillVal,miscount=miscount)
                        save_as_mapsstack_per_day(lats,lons,tmax,int(ncnt),odir,prefix='TMAX',oformat=oformat,FillVal=FillVal,miscount=miscount)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[4],int(ncnt),odir,prefix='RSIN',oformat=oformat,FillVal=FillVal,miscount=miscount)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[0],int(ncnt),odir,prefix='TEMP',oformat=oformat,FillVal=FillVal,miscount=miscount)


                if evapMethod == 'Hargreaves':
                    # retrieve 3 hourly Temperature and calculate max and min Temperature
                    temperatureFilename = temperatureFile[0]
                    standard_name = 'air_temperature'
                    timestepSeconds = 10800

                    logger.info("Get times 3 hr data..")
                    tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, temperatureFilename,timestepSeconds,logger )
                    logger.info("Get actual 3hr data...")
                    ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)
                    mstack = ncstepobj.getdates(timelist)
                    tmin = flipud(mstack.min(axis=0))
                    tmax = flipud(mstack.max(axis=0))

                    if downscaling == 'True' or resampling == "True":
                        tmin = resample_grid(tmin,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method=resamplingtype,FillVal=nan)
                        tmax = resample_grid(tmax,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method=resamplingtype,FillVal=nan)
                        if resampling == "True":
                            tmin = correctTemp(tmin,lapse_rate,elevationCorrection)
                            tmax = correctTemp(tmax,lapse_rate,elevationCorrection)
                        tmax[mismask] = nan
                        tmin[mismask] = nan

                    logger.info("Start hargreaves..")
                    PETmm = hargreaves(LATITUDE,currentdate,relevantDataFields, tmax, tmin)
                    PETmm[mismask] = nan
                    PETmm[PETmm < -10.0] = nan
                    PETmm[PETmm > 135.0] = nan

                    PETmm = fixmapwithmask(PETmm,demmask)
                    tmin = fixmapwithmask(tmin, demmask)
                    tmax = fixmapwithmask(tmax, demmask)
                    relevantDataFields[0] = fixmapwithmask(relevantDataFields[0], demmask)

                    logger.info("Saving Hargreaves PET data for: " +str(currentdate))

                    save_as_mapsstack_per_day(lats,lons,PETmm,int(ncnt),odir,prefix=oprefix,oformat=oformat,FillVal=FillVal,miscount=miscount)
                    if saveAllData:
                        save_as_mapsstack_per_day(lats,lons,tmin,int(ncnt),odir,prefix='TMIN',oformat=oformat,FillVal=FillVal,miscount=miscount)
                        save_as_mapsstack_per_day(lats,lons,tmax,int(ncnt),odir,prefix='TMAX',oformat=oformat,FillVal=FillVal,miscount=miscount)
                        save_as_mapsstack_per_day(lats,lons,relevantDataFields[0],int(ncnt),odir,prefix='TEMP',oformat=oformat,FillVal=FillVal,miscount=miscount)

        else:
            pass

        currentdate += datetime.timedelta(days=1)

    logger.info("Done.")

if __name__ == "__main__":
    main()