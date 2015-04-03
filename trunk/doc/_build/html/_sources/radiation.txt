Radiation correction on a digital elevation model
=================================================


Introduction
------------
The radiation module can be used to create celar sky radiation estimates
using a (high resolution) digital elevation model. In turn, these results can be used
by the evaporation module to better estimate local evaporation.


Usage
=====

The folowing command takes the digital elevation model wflow_dem.map and calcilate the
radiation components for day 1 to 366 and save these as a pcraster mapstack. The internal
timesteps is set to 60 minutes and the calculations are done for 5 hr in the morning to 22 hour at night.

::

    e2o_radiation.py -S 1 -E 366 -D wflow_dem.map -O output_rad -s 5 -e 22 -T 60 -f PCRaster



Implementation
==============

.. automodule:: e2o_dstools.e2o_radiation
    :members:




Description
-----------

Adjusted after van Dam 2000

This section gives a short description. Another description can
be found at http://re.jrc.cec.eu.int/pvgis/pv/solres/solres.htm.

Potential solar radiation is the radiation of an unobstructed or
cloudless sky. The magnitude of the potential solar radiation depends
on the position of the sun the solar altitude or solar angle during the day,
the inclination of the solar rays with the earth’s surface, the amount of
radiation at the outer layer of the earth’s atmosphere, the transmissivity
of the sky and the altitude of the earth’s surface.

Solar declination is the annual fluctuation of the sun between the two
tropics and varies between –23 and +23 degrees latitude. Solar declination is
calculated per Day (Julian day number):


.. math::

    \delta = -23.4 cos(360 (Day + 10 )/365)

The hour angle describes the movement of the sun around the earth in 24 hours,
which  equals 15 degrees longitude per hour (360 :math:`deg` /24h). The hour angle n
is calculated for each Hour (whole hour of the day):

.. math::

    n = ( 15 \dot Hour - 12)


The position or height of the sun above the horizon is called the solar altitude
or solar angle. Solar altitude :math:`\alpha` (deg) is calculated for each location,
determined by the location’s latitude :math:`\psi` (deg), declination and hour angle:


.. math::

    sin(\alpha) = sin(\psi) sin(\delta) + cos(\psi) cos(\delta) cos(n)

Solar azimuth is the angle between the solar rays and the North-South axis of the
earth. Solar azimuth :math:`{\beta}_s` (deg) is calculated by:

.. math::

    &cos({\beta}_s) = (sin(\delta) cos(\psi) - cos(\delta) sin(\psi) cos(n))/cos(\alpha) \\
    &for Hour \le 12: {\beta}_s = {\beta}_s \\
    &for Hour > 12: {\beta}_s = 360 - {\beta}_s

Surface azimuth or aspect :math:`{\beta}_1` (deg) is the orientation of the land
surface or slope to the North-South axis of the sun. Slope :math:`\varphi` (deg) is
the maximum rate of change in elevation.

The angle of incidence is the angle between the perpendicular plane of the
incoming solar rays and the surface on which they are projected, defined by the
aspect and slope of that surface. The angle of incidence :math:`\vartheta` (deg) is
calculated with the solar angle :math:`\alpha` (deg), the slope of the land
surface :math:`\varphi` (deg), the azimuth of the sun :math:`{\beta}_s` (deg) and
azimuth of the land surface  :math:`{\beta}_1` (deg):

.. math::

    cos(\vartheta) = cos(\alpha) sin(\varphi) cos({\beta}_s - {\beta}_1) + sin(\alpha) cos(\varphi)

The second section of the radiation module calculates the potential solar energy. The
amount of solar radiation that reaches the outer atmosphere is decreased by the
travelling distance of the solar rays through the sky to the surface, the transmissivity
of the sky and the cloud factor.

Solar energy at the outer layer of the atmosphere :math:`Sout (Wm^2)` is
calculated by (Kreider & Kreith 1975):

.. math::

    S_{out} = S_c (1 + 0.034 cos(360 Day/365))

where :math:`S_c  (Wm^2)` is the solar constant of 1367 :math:`Wm^2` (Duffie & Beckman 1991).
The solar ‘constant’ is subject to much discussion. Gates (1980) gives a value
of 1360 :math:`Wm^2`. The NASA reports a value of 1353 :math:`Wm^2` (Jansen 1985),
while Duncan et al. (1982) give a value of 1367 :math:`Wm^2`. Monteith and Unsworth (1990)
measured the highest value of 1373 W.m-2. The World Radiation Centre uses a
value of 1367 :math:`Wm^2` (Duffie & Beckman 1991) and this value is also used in this study.

The solar radiation energy that reaches the earth’s surface is decreased due to the
length of the air mass it has to pass through and the transmissivity :math:`\tau`
(% or fraction) of the sky. The radiation flux through a hypothetical plane
normal to the beam  (:math:`S_{nor} Wm^2`) is given by (Gates 1980):

.. math::

    S_{nor} = S_{out} \tau^{Mh}

in which Mh (% or fraction) is the relative path length of the optical air mass at
altitude h (m). Transmissivity (:math:`\tau`) is usually between 0.5 and 0.8, but can be as low
as 0.4 in the tropics (Whitmore et al. 1993), but mostly a value of 0.6 is used
(Gates 1980). To calculate the relative path length of an optical air mass at
altitude h (m), the relative path length of an optical air mass at sea level M0
(% or fraction) is corrected for the atmospheric pressure at altitude h. Mh
(% or fraction) is calculated using (Kreider & Kreith 1975):

.. math::

    Mh = M_0 P_h/P_0

in which :math:`P_h / P_0` (mbar.mbar-1) is an atmospheric pressure correction. The relative path length
of the optical air mass at sea level M0 is obtained by (Kreider & Kreith 1975):

.. math::

    M_0 = \sqrt(1299 + (614 sin(\alpha))^2) - 614 sin(\alpha)

The atmospheric pressure correction :math:`P_h / P_0` is written as (List 1984):


.. math::

    P_h / P_0 = ((288 0.0065h) / 288)^5.256


The incoming radiation normal to the beam Snor must be corrected by the orientation and
slope of the surface, defined by the angle of incidence :math:`\vartheta`, to calculate the incoming radiation
Sdir (:math:`Wm^2`) on the earth’s surface:

.. math::

    S_dir = S_nor cos(\vartheta)

Direct light is scattered in the atmosphere. This daylight scattering or diffuse radiation is
approximately 15% of direct radiation (Gates 1980). A more accurate empirical estimation for diffuse
radiation Sdif (:math:`Wm^2`) in a clear not dust-free sky reads as (Liu and Jordan in Gates 1980):


.. math::

    S_dif = S_out (0.271 - 0.294 \vartheta^{Mh} sin(\alpha)

During daylight when the sun is above the horizon, it is assumed that all cells receive the same amount
of diffuse radiation. Total incoming radiation Sin (:math:`Wm^2`) is the sum of direct and diffuse radiation:

.. math::

    S_in = S_dir + S_dif


Total incoming radiation Sin as calculated with the above is actually a radiation flux for that moment. In the procedure
given above, radiation is calculated per time step. If this amount of radiation is used in a water balance model, the amount
of radiation and therewith the amount of evapotranspiration will be overestimated or under estimated, depending on the
time of the day and the position of the sun.


Most of the work done for the shading is implemented in the pcraster horizontan function.



Example
=======

The image below show the difference between the average daily clear sky radiation on
a flat surface compared to the actual radiation received by each grid cell for
julian day 180. IN the steep terrain of the Snowy Mountains of the mUrumbidgee catchment
differences can be up to 30 :math:`W/m^2` on this 1x1km DEM.

.. figure:: _static/murumbidgee.png
    :width: 640px
    :align: center

    Difference in :math:`W/m^2` between horizontal surface radiation and
    inclined surface radiation



