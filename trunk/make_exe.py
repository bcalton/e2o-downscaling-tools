from distutils.core import setup
from bbfreeze import Freezer
from _version import *
import ctypes

def dependencies_for_freeezing():
	import netCDF4_utils 

nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)



f = Freezer("e2o_dstools_" + str(nrbits) +"_bit")
f.addScript("e2o_dstools/e2o_radiation.py")
f.addScript("e2o_dstools/e2o_getvar.py")
f.addScript("e2o_dstools/e2o_calculateEvaporation.py")
f()    # starts the freezing process
