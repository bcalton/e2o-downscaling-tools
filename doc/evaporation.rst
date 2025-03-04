

Downscaling
===========

Radiation
---------
The paragraph below is adapted from the r.sun grass manual:

The real-sky irradiance/irradiation are calculated from clear-sky raster maps by the
application of a factor parameterizing the attenuation of cloud cover. Examples of explicit
calculations of this parameter can be found in Becker (2001), Kitler and Mikler (1986). However, the cloudiness
observation by a meteorological service routine is usually prone to subjective errors and does
not describe sufficiently the physical nature and dynamic spatial-temporal pattern of different
types of cloud cover. Therefore, a simpler parameter has to be used. The solutions for horizontal
and inclined surfaces are slightly different. For the assessment of global irradiance/irradiation
on a horizontal surface under overcast conditions Gh the clear-sky values Ghc are multiplied by
clear-sky index kc (Beyer et al 1996, Hammer et al 1998, Rigollier et al. 2001):

.. math::

	Gh = Ghc kc

The index kc represents the atmospheric transmission expressed as a ratio between horizontal
global radiation under overcast and clear-sky conditions. For a set of ground meteorological
stations the clear-sky index can be calculated from measured global radiation Ghs and
computed values of clear-sky global radiation Ghc:

.. math::

	kc = Ghs/Ghc

As an alternative the kc can be derived also from other climatologic data
(e.g. cloudiness, cf. Kasten and Czeplak 1980). The raster maps of kc must be
then derived by spatial interpolation. The kc can be calculated directly as a raster map from
short-wave surface irradiance measured by satellites. This method is based on the complementarity
between the planetary albedo recorded by the radiometer and the surface radiant flux
(Cano et al 1986, Beyer et al 1996, Hammer et al 1998).
To compute the overcast global irradiance/irradiation for inclined surfaces, Gi
the diffuse Dh and beam Bh components of overcast global radiation and of the clear-sky index kc
have to be treated separately as follows from the following equations:

.. math::

	Dh = Dhc kdc

	Bh = Bhc kbc

The ratio of diffuse to the global radiation Dh/Gh for clear and overcast skies changes
according to the cloudiness. In Europe the Dh/Gh values are typically in interval 0.3-1.0
(Kasten and Czeplak 1980). The underlying physical processes are quite complicated and computationally
represented only by empirical equations (cf. Scharmer and Greif, 2000, Kasten and Czeplak 1980, Hrvoľ 1991).
However, for many meteorological stations, besides the global horizontal radiation Ghs, the diffuse component
Dhs is either measured or calculated from cloudiness, sunshine or other climatologic data.
The raster map of Dhs/Ghs can be derived from the point values by spatial interpolation.
Consecutively, the raster maps of diffuse and beam components of the clear sky index can be computed:

.. math::

	Dh = Gh Dhs/Ghs

	Bh = Gh – Dh

	kdc = Dh/Dh

	kbc = Bh/Bhc


where subscript s is meant to distinguish data measured on meteorological stations Bhs
nd Dhs from the estimated values Bh, and Dh.

Temperature
-----------

Temperature is downscaled using a laps-rate which is fixed at 0.006 K m-1. The down-scaling is based on the
difference in elevation between (1) the global low-resolution DEM belonging to the eartH2Observe dataset and (2) the
local high resolution DEM. The down-scaled temperature is calculated for each grid cell of the high resolution DEM.
Min and max temperature: Daily minimum and maximum temperature are calculated by taking the daily minimum and maximum
 values from the 3-hourly temperature time-series.


Pressure
--------

Pressure is down-scaled with the barometric formula. The barometric formula gives the pressure in the atmosphere as a
function of height. Since temperature and the composition of the atmosphere are complicated functions of height, and
because gravitation is an inverse function of the distance to the centre of the planet, three simplifications have
been made within the version of the equation we apply:

- temperature, gravitation, and composition are assumed constant throughout the atmosphere;
- the atmosphere is assumed an ideal gas;
- temperature is assumed to decrease linearly with height with a slope L, the atmospheric lapse rate.

The downscaling is based on the difference (ΔH [m]) between the low resolution global and the high resolution local
DEM according to the following equation:

.. math::

        P_{cor} = P \left(\frac{\bar{T}}{\bar{T} + L \Delta H}^\frac{g M_0}{L R} \right)



Where P = uncorrected pressure from the eartH2Observe dataset [Pa], $\bar{T}$  = the daily average temperature from the
eartH2Observe dataset [K], g = 9.81 - gravitational constant [m s-2], R = 8.3144621 - specific gas constant for dry
air [J mol-1 K-1], Mo = 0.0289644 - molecular weight of gas [g / mol] and L = 0.006 - lapse rate [K m-1].



Evaporation methods
===================

The following reference evaporation equations have been implemented:

- Hargreaves - temperature based
- Priestley-Taylor - temperature and radiation based
- Penman-Monteith – fully physically based


These are defined as follows:

Hargreaves:

.. math::

    E{T_o} = 0.0023 \cdot {R_a} \cdot (\overline T  + 17.8) \cdot {(TR)^{0.50}}

Priestley-Taylor:

.. math::

    E{T_o} = \alpha \frac{{\Delta {R_n}}}{{{\lambda _v}(\Delta  + \gamma )}}


Penman-Monteith (at z= 10m):

.. math::

    E{T_o} = \frac{{\Delta ({R_n} - G) + {\rho _a}{c_p}\frac{{({e_s} - {e_a})}}{{{r_a}}}}}{{{\lambda _v}\Delta  + \gamma (1 + \frac{{{r_s}}}{{{r_a}}})}}




Where λv = Latent heat of vaporization (J/g), Δ = the slope
of the saturation vapour pressure temperature relationship (Pa K-1), Rn = Net radiation (W m-2), G = Soil heat flux
(W m-2), cp = specific heat of the air (J kg-1 K-1), ρa = mean air density at constant pressure (kg m-3), es - ea =
vapor pressure deficit (Pa), rs = surface resistances (m s-1), ra = aerodynamic resistances (m s-1), γ =
Psychrometric constant (66 Pa K-1), Ra = extraterrestrial radiation (MJ m-2 day-1), = mean daily temperature (°C),
TR= temperature range (°C), α = empirical multiplier (-;1.26),




Implementation
==============

ini file configuration
----------------------

The .ini file below shows the available options

.. literalinclude:: _download/e2o_calculateEvaporation.ini

The file can be downloaded here: :download:`here. <_download/e2o_calculateEvaporation.ini>`



Code description
----------------

.. automodule:: e2o_dstools.e2o_calculateEvaporation
    :members:

