import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from setuptools import find_packages
here = os.path.abspath(os.path.dirname(__file__))

README = open(os.path.join(here, 'README.txt')).read()

requires = [
    'netCDF4',
    'numpy',
    'matplotlib',
    'pcraster', 'osgeo','scipy']



# Source dist
setup(name='e2o_dstoools',
      version= "0.1",
      packages=['e2o_dstools'],
      package_dir={'e2o_dstools': 'e2o_dstools'},
      author='J. Schellekens/F Sperna Weiland',
      author_email='jaap.schellekens@deltares.nl',
      url='http://www.earth2observe.eu',
      license = "GPL",
      scripts=['e2o_dstools/e2o_getvar.py','e2o_dstools/e2o_calculateEvaporation.py'],
      description='Download and downscaling tools for the earth2observe datasets',
      )

