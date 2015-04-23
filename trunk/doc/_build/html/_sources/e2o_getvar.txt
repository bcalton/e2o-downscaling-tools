e2o_getvar - retrieve variables
===============================

e2o_getval gett variables from the server for a specific region (specified in the ini file)
and optionally resamples those to a specified Digital Elevation Model (DEM).

Usage:
::

    e2o_getvar.py -I the_ini_file


ini file configuration
----------------------
The .ini file below shows the available options

.. literalinclude:: _download/e2o_getvar.ini


An example ini file be found :download:`here. <_download/e2o_getvar.ini>`


Implementation
--------------



.. automodule:: e2o_dstools.e2o_getvar
    :members:

