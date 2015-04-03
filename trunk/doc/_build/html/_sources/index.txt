===============================
eartH2Observe Downscaling tools
===============================

.. note::
      This documentation is for version |version| of e2o_dstools, release |release|
      This documentation was generated |today|

Introduction
============

e2o_downscaling-tools consists of a number of python programs and procedures that facilitate local
application of the earth2observe global water resources reanalysis. The tools
can connect directly to the project's data server and save (resampled) data to a local computer
for further analysis or direct application. The current first versions
of the tool focusses on downscaling the global forcing dataset used in
the project :cite:`weedonwfdei2014`.



Usage
=====

PET Determination and downscaling
---------------------------------

The tools assume you have a digital elevation model fo3 your area. This file should be in
a GDAL supported format (preferably GTiff).

+ Optionally, first run the e2o_radiation script. This will generate Clear-Sky radiation maps for each day of the
  year (four maps per day). These maps can be used by the e2o_calculateEvaportion script to downscalle
  reference ET
+ Next run the e2o_calculateEvaporation script. This will calculate downscaled ET based on a local DEM for
  the priod you specify in the .ini file

See the documentation per module for more information.


Retrieving variables
--------------------
The e2o_getvar script allows you to retrieve single or multiple variable from the e2o server fro
a specified region and timespan.

The radiation module
====================
.. toctree::
   :maxdepth: 2

   radiation

The evaporation module
======================
.. toctree::
   :maxdepth: 2

   evaporation


The e2o_getvar script
=====================
.. toctree::
   :maxdepth: 2

   e2o_getvar

Examples and tests
==================
.. toctree::
   :maxdepth: 2
   

FAQ
===
.. toctree::
   :maxdepth: 2

   faq

   
Release notes
=============
.. toctree::
   :maxdepth: 2

   release-notes


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



TODO
====

.. todolist::


.. only:: html

   .. rubric:: References

.. bibliography:: e2o.bib
   :style: alpha
