Potential evaporation
=====================



Introduction
------------


In general evaporation amounts are determined for about 90% by
radiation input. Radiation at the earth’s surface is determined by
the potential solar radiation at the edge of the earths atmosphere
and the filtering within the atmosphere. The first component can be
easily determined from equations. The reduction due to clouds etc
can be estimated by incorporation short wave measurements but if
these are not available cloud cover estimates can also be used. By
combining this with a DEM radiation at the earths surface can be
determined including the effects of aspect and shading.

The figure below shows the steps used to generated down-scale reference evaporation.

.. digraph:: steps

    "Clear sky radiation maps" [shape=box];
    "Reference evaporation" [shape=box];
    "e2o_radiation.py" -> "Clear sky radiation maps" [label =" Correct for aspect and slope with DEM"];
     "e2o_calculateEvaporation.py" -> "Reference evaporation" [label =" Downscale using DEM and clear-sky maps"]
    "Clear sky radiation maps" -> "e2o_calculateEvaporation.py"

    dpi=69;

ini file configuration
----------------------

The .ini file below shows the available options

.. literalinclude:: _download/e2o_calculateEvaporation.ini

The file can be downloaded here: :download:`here. <_download/e2o_calculateEvaporation.ini>`


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

Temperature is downscaled using a fixed laps-rate.

Pressure
--------

Pressure is downscaled .... PM


Implementation
==============





.. automodule:: e2o_dstools.e2o_calculateEvaporation
    :members:

