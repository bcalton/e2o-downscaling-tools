#!/usr/bin/python

# e2o_dstools is Free software, see below:
#
# Copyright (c) Deltares 2005-2014
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



# $Rev:: 904           $:  Revision of last commit
# $Author:: schelle    $:  Author of last commit
# $Date:: 2014-01-13 1#$:  Date of last commit

"""
Determine clear sky radiation over a Digital Elevation Model.

Usage::

    e2o_radiation -D DEM [-O outputdir][-S start day][-E end day]
              [-M][-x lon][-y lat][-h][-l loglevel][-T minutes]

    -D DEM Filename of the digital elevation model
    -O outputdir (default is . )
    -S Startday - Start day of the simulation (1 Jan is 1)
    -E EndDay - End day of the simulation
    -T minutes - timeresolution in minutes (60 default is 1 hour)
    -M The DEM xy units are in metres (instead of lat/lon)
    -x longitute of the map left (if map xy in metres)
    -y lattitude of the map bottom (if map xy in metres)
    -l loglevel Set loglevel (DEBUG, INFO, WARNING,ERROR)
    -s start hour (per day) of the calculations (default =1)
    -e end hour (per day) of the calculations (default = 23)
    -h This information
    -f output format as a gdal type string (http://www.gdal.org/formats_list.html) Default is PCRaster
    -p add the file format string as a postfix to all filename (default is False)


The program produces the following map stacks, one for each day of
the year:

::

    COR00000.??? - Total clear sky radiation on DEM
    SUN00000.??? - Nr of time intervals a pixel was in the sun
    FLAT0000.??? - Total clear sky radiation on a flat surface
    CORDIR00.??? - Direct clear sky radiation on DEM
    FLATDIR0.??? - Direct clear sky radiation on a flat surface

"""


from pcraster import *
import sys
import logging
import e2o_utils
import getopt
import numpy as np


def lattometres(lat):
    """
    Determines the length of one degree lat/long at a given latitude (in meter).
    Code taken from http:www.nga.mil/MSISiteContent/StaticFiles/Calculators/degree.html

    :param lat:  map with lattitude values for each cell
    :return: length of a cell lat, length of a cell long
    """

    #radlat = spatial(lat * ((2.0 * math.pi)/360.0))
    #radlat = lat * (2.0 * math.pi)/360.0
    radlat = spatial(lat) # pcraster cos/sin work in degrees!


    m1 = 111132.92        # latitude calculation term 1
    m2 = -559.82        # latitude calculation term 2
    m3 = 1.175            # latitude calculation term 3
    m4 = -0.0023        # latitude calculation term 4
    p1 = 111412.84        # longitude calculation term 1
    p2 = -93.5            # longitude calculation term 2
    p3 = 0.118            # longitude calculation term 3
    # # Calculate the length of a degree of latitude and longitude in meters

    latlen = m1 + (m2 * cos(2.0 * radlat)) + (m3 * cos(4.0 * radlat)) + (m4 * cos(6.0 * radlat))
    longlen = (p1 * cos(radlat)) + (p2 * cos(3.0 * radlat)) + (p3 * cos(5.0 * radlat))

    return latlen, longlen

def detRealCellLength(ZeroMap,sizeinmetres):
    """
    Determine cellength. Always returns the length
    in meters.

    :param ZeroMap: pcraster map object used as clone
    :param sizeinmetres: is set to one the currect cell size is assumed to be in metres, otherwise in lat,log
    :return:
    """

    if sizeinmetres:
            reallength = celllength()
            xl = celllength()
            yl = celllength()
    else:
        aa = ycoordinate(boolean(cover(ZeroMap + 1,1)))
        yl, xl = lattometres(aa)

        xl = xl * celllength()
        yl = yl * celllength()
        # Average length for surface area calculations.

        reallength = (xl + yl) * 0.5

    return xl,yl,reallength


def correctrad(Day,Hour,Lat,Lon,Slope,Aspect,Altitude,Altitude_UnitLatLon):
    """ 
    Determines radiation over a DEM assuming clear sky for a specified hour of
    a day
    
    :var Day: Day of the year (1-366)
    :var Hour: Hour of the day (0-23)
    :var Lat: map with latitudes for each grid cell
    :var Lon: map with longitudes for each grid cell
    :var Slope: Slope in degrees
    :var Aspect: Aspect in degrees relative to north for each cell
    :var Altitude: Elevation in metres
    :var Altitude_Degree: Elevation in degrees. If the actual pcraster maps
                          are in lat lon this maps should hold the Altitude converted
                          to degrees. If the maps are in metres this maps should also
                          be in metres

    :return Stot: Total radiation on the dem, shadows not taken into account
    :return StotCor: Total radiation on the dem taking shadows into acount
    :return StotFlat: Total radiation on the dem assuming a flat surface
    :return SUN: Map with shade (0) or no shade (1) pixels
    """

    Sc  = 1367.0          # Solar constant (Gates, 1980) [W/m2]
    Trans   = 0.6             # Transmissivity tau (Gates, 1980)    
    pi = 3.1416
    a = pow(100,5.256)
    #report(a,"zz.map")

    AtmPcor = pow(((288.0-0.0065*Altitude)/288.0),5.256)
    #Lat = Lat * pi/180
    ##########################################################################
    # Calculate Solar Angle and correct radiation ############################
    ##########################################################################
    # Solar geometry
    # ----------------------------
    # SolDec  :declination sun per day  between +23 & -23 [deg]
    # HourAng :hour angle [-] of sun during day
    # SolAlt  :solar altitude [deg], height of sun above horizon
    # SolDec  = -23.4*cos(360*(Day+10)/365);
    # Now added a new function that should work on all latitudes! 
    #theta    =(Day-1)*2 * pi/365  # day expressed in radians
    theta    =(Day-1)*360.0/365.0  # day expressed in degrees
     
    SolDec =180/pi * (0.006918-0.399912 * cos(theta)+0.070257 * sin(theta) - 0.006758 * cos(2*theta)+0.000907 * sin(2*theta) - 0.002697 * cos(3*theta)+0.001480 * sin(3*theta))
    
    #HourAng = 180/pi * 15*(Hour-12.01)
    HourAng = 15.0*(Hour-12.01) 
    SolAlt  = scalar(asin(scalar(sin(Lat)*sin(SolDec)+cos(Lat)*cos(SolDec)*cos(HourAng))))
    
    # Solar azimuth                    
    # ----------------------------
    # SolAzi  :angle solar beams to N-S axes earth [deg]
    SolAzi = scalar(acos((sin(SolDec)*cos(Lat)-cos(SolDec)* sin(Lat)*cos(HourAng))/cos(SolAlt)))
    SolAzi = ifthenelse(Hour <= 12, SolAzi, 360 - SolAzi)
    

    # Surface azimuth
    # ----------------------------
    # cosIncident :cosine of angle of incident; angle solar beams to angle surface
    cosIncident = sin(SolAlt)*cos(Slope)+cos(SolAlt)*sin(Slope)*cos(SolAzi-Aspect)
    # Fro flat surface..  
    FlatLine = spatial(scalar(0.00001))
    FlatSpect = spatial(scalar(0.0000))
    cosIncidentFlat = sin(SolAlt)*cos(FlatLine)+cos(SolAlt)*sin(FlatLine)*cos(SolAzi-FlatSpect)
    # Fro flat surface..    
    #cosIncident = sin(SolAlt) + cos(SolAzi-Aspect)

    # Critical angle sun
    # ----------------------------
    # HoriAng  :tan maximum angle over DEM in direction sun, 0 if neg 
    # CritSun  :tan of maximum angle in direction solar beams
    # Shade    :cell in sun 1, in shade 0
    # NOTE: for a changing DEM in time use following 3 statements and put a #
    #       for the 4th CritSun statement
    HoriAng   = cover(horizontan(Altitude_UnitLatLon,directional(SolAzi)),0)
    #HoriAng   = horizontan(Altitude,directional(SolAzi))
    HoriAng   = ifthenelse(HoriAng < 0, scalar(0), HoriAng)
    CritSun   = ifthenelse(SolAlt > 90, scalar(0), scalar(atan(HoriAng)))
    Shade   = SolAlt > CritSun
    #Shade = spatial(boolean(1))
    # Radiation outer atmosphere
    # ----------------------------
    #report(HoriAng,"hor.map")

    OpCorr = Trans**((sqrt(1229+(614*sin(SolAlt))**2) -614*sin(SolAlt))*AtmPcor)    # correction for air masses [-] 
    Sout   = Sc*(1+0.034*cos(360*Day/365.0)) # radiation outer atmosphere [W/m2]
    Snor   = Sout*OpCorr                   # rad on surface normal to the beam [W/m2]

    # Radiation at DEM
    # ----------------------------
    # Sdir   :direct sunlight on dem surface [W/m2] if no shade
    # Sdiff  :diffuse light [W/m2] for shade and no shade
    # Stot   :total incomming light Sdir+Sdiff [W/m2] at Hour
    # Radiation :avg of Stot(Hour) and Stot(Hour-HourStep)
    # NOTE: PradM only valid for HourStep & DayStep = 1

    
    SdirCor   = ifthenelse(Snor*cosIncident*scalar(Shade)<0,0.0,Snor*cosIncident*scalar(Shade))
    Sdir   = ifthenelse(Snor*cosIncident<0,0.0,Snor*cosIncident)
    SdirFlat   = ifthenelse(Snor*cosIncidentFlat<0,0.0,Snor*cosIncidentFlat)
    Sdiff  = ifthenelse(Sout*(0.271-0.294*OpCorr)*sin(SolAlt)<0, 0.0, Sout*(0.271-0.294*OpCorr)*sin(SolAlt))
    #AtmosDiffFrac = ifthenelse(Sdir > 0, Sdiff/Sdir, 1)
    Shade = ifthenelse(Sdir <=0, 0,Shade)

    # Stot   = cover(Sdir+Sdiff,windowaverage(Sdir+Sdiff,3));     # Rad [W/m2]
    Stot   = Sdir + Sdiff                                             # Rad [W/m2]
    StotCor   = SdirCor + Sdiff                                   # Rad [W/m2]
    StotFlat = SdirFlat + Sdiff
    
    
     
    return StotCor, StotFlat, Shade, SdirCor, SdirFlat


def GenRadMaps(SaveDir, Lat, Lon, Slope, Aspect, Altitude, DegreeDem, logje, start=1, end=2, interval=60, shour=1, ehour=23, outformat='PCRaster',Addpostfix=False):
    """
    Generate radiation masp for a number of days

    :param SaveDir: When to save the maps
    :param Lat: LAttitude for each pixel
    :param Lon: Longitude for each pixel
    :param Slope: Slope for each pixel
    :param Aspect: Aspect
    :param Altitude: Altitude
    :param DegreeDem: Altitude with elevation rescaled to degree lat/lon
    :param logje: Logger
    :param start: Start day of the year
    :param end: End day of the year
    :param interval: Interval for calculations in minutes
    :param shour: Start hour of the calculations
    :param outformat: File format of results as a gdal string
    :return: Nothing
    """

    Intperday = 1440./interval
    Starthour = shour
    EndHour = ehour
    Calcsteps = Intperday/24 * 24
    calchours = np.arange(Starthour,EndHour,24/Intperday)

    for Day in range(start,end+1):
        nr = "%0.3d" % Day
        # check if step already existst
        ckfile = SaveDir + "/COR00000." + nr
        if os.path.exists(ckfile):
            logging.warn(ckfile + " exists, skipping this step...")
        else:
            avgrad = 0.0 * Altitude
            _flat = 0.0 * Altitude
            avshade = 0.0 * Altitude

            cordir = 0.0 * Altitude
            flatdir = 0.0 * Altitude
            id = 1
            logje.info("Calulations for day: " + str(Day))
            for Hour in calchours:
                logje.debug("Hour: " + str(Hour))
                crad,  flat, shade, craddir, craddirflat = correctrad(Day,float(Hour),Lat,Lon,Slope,Aspect,Altitude,DegreeDem)
                avgrad=avgrad + crad
                _flat = _flat + flat
                avshade=avshade + scalar(shade)
                cordir = cordir + craddir
                flatdir = flatdir + craddirflat
                id = id + 1


            miss = float(1E31)
            if Addpostfix:
                postfix = "." + outformat
            else:
                postfix = ""

            x = pcr2numpy(xcoordinate(boolean(cover(1.0))),1E31)[0,:]
            y = pcr2numpy(ycoordinate(boolean(cover(1.0))),1E31)[:,0]

            data = pcr2numpy(avgrad/Calcsteps,1E31)
            e2o_utils.writeMap(os.path.join(SaveDir,"COR00000." + nr + postfix),outformat,x,y,data, 1E31)
            data = pcr2numpy(avshade,1E31)
            e2o_utils.writeMap(os.path.join(SaveDir,"SUN00000." + nr+ postfix),outformat,x,y,data, 1E31)
            data = pcr2numpy(_flat/Calcsteps,1E31)
            e2o_utils.writeMap(os.path.join(SaveDir,"FLAT0000." + nr+ postfix),outformat,x,y,data, 1E31)
            data = pcr2numpy(cordir/Calcsteps,1E31)
            e2o_utils.writeMap(os.path.join(SaveDir,"CORDIR00." + nr+ postfix),outformat,x,y,data, 1E31)
            data = pcr2numpy(flatdir/Calcsteps,1E31)
            e2o_utils.writeMap(os.path.join(SaveDir,"FLATDIR0." + nr+ postfix),outformat,x,y,data, 1E31)


def usage(*args):
    """

    :param args:
    :return:
    """

    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)


def main(argv=None):
    """

    :param argv: See ussage
    :return:
    """

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return


    try:
        opts, args = getopt.getopt(argv, 'hD:Mx:y:l:O:S:E:T:s:e:f:p')
    except getopt.error, msg:
        usage(msg)


    thedem = "mydem.map"
    xymetres = False
    lat = 52
    lon = 10
    loglevel = logging.INFO
    outputdir="output_rad"
    startday = 1
    endday = 2
    calc_interval = 60
    shour=1
    ehour=23
    oformat ='PCRaster'
    postfix =False


    for o, a in opts:
        if o == '-h': usage()
        if o == '-O': outputdir = a
        if o == '-D': thedem = a
        if o == '-M': xymetres = true
        if o == '-x': lat = int(a)
        if o == '-y': lon = int(a)
        if o == '-S': startday = int(a)
        if o == '-E': endday = int(a)
        if o == '-T': calc_interval = int(a)
        if o == '-l': exec "loglevel = logging." + a
        if o == '-s': shour = int(a)
        if o == '-e': ehour = int(a)
        if o == '-f': oformat = a
        if o == '-p': postfix = True


    logger = e2o_utils.setlogger("e2o_radiation.log","e2o_radiation",level=loglevel)
    if not os.path.exists(thedem):
        logger.error("Cannot find dem: " + thedem + " exiting.")
        sys.exit(1)
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    logger.debug("Reading dem: " " thedem")
    setclone(thedem)
    dem = readmap(thedem)

    logger.debug("Calculating slope and aspect...")
    if xymetres:
        LAT = spatial(scalar(lat))
        LON= spatial(scalar(lon))
        Slope = max(0.00001,slope(dem))
        DEMxyUnits = dem
    else:
        LAT= ycoordinate(boolean(dem))
        LON = xcoordinate(boolean(dem))
        Slope = slope(dem)
        xl, yl, reallength = detRealCellLength(dem * 0.0, 0)
        Slope = max(0.00001, Slope * celllength() / reallength)
        DEMxyUnits = dem * celllength() / reallength

    # Get slope in degrees
    Slope = scalar(atan(Slope))
    Aspect = cover(scalar(aspect(dem)),0.0)

    GenRadMaps(outputdir,LAT,LON,Slope,Aspect,dem,DEMxyUnits,logger,start=startday,end=endday,interval=calc_interval,shour=shour,ehour=ehour,outformat=oformat,Addpostfix=postfix)


if __name__ == "__main__":
    main()
