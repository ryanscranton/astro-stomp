#!/usr/bin/env python

"""
setup.py file for STOMP library
"""

import sys, os
import subprocess

from distutils.core import setup, Extension

# make sure we compile in the python-specific
# methods to the stomp classes

addflags='-DWITH_PYTHON'
if 'CPPFLAGS' in os.environ:
    os.environ['CPPFLAGS'] +=  ' ' +addflags
else:
    os.environ['CPPFLAGS'] = addflags

# also try to add numpy support
try:
    import numpy
    include_dirs=[numpy.get_include()]
    addflags='-DWITH_NUMPY'
    os.environ['CPPFLAGS'] +=  ' ' +addflags
    depends=['NumpyVector.h']
except:
    sys.stdout.write("Numpy not found, not building with numpy support\n")
    include_dirs=[]
    depends=[]



# create the ups table
pyvers='%s.%s' % sys.version_info[0:2]
d1='lib/python%s/site-packages' % pyvers
d2='lib64/python%s/site-packages' % pyvers

if not os.path.exists('ups'):
    os.mkdir('ups')
tablefile=open('ups/stomp.table','w')
tab="""
setupOptional("gflags")
envPrepend(PATH,${PRODUCT_DIR}/bin)
envPrepend(LD_LIBRARY_PATH,${PRODUCT_DIR}/lib)
envPrepend(LIBRARY_PATH,${PRODUCT_DIR}/lib)
envPrepend(C_INCLUDE_PATH,${PRODUCT_DIR}/include)
envPrepend(CPATH,${PRODUCT_DIR}/include)
envPrepend(PYTHONPATH,${PRODUCT_DIR}/%s)
envPrepend(PYTHONPATH,${PRODUCT_DIR}/%s)
""" % (d1,d2)
tablefile.write(tab)
tablefile.close()



stomp_module = Extension("_stomp",
                         depends=depends,
                         sources=["../stomp/stomp_core.cc",
                                  "../stomp/stomp_angular_bin.cc",
                                  "../stomp/stomp_angular_coordinate.cc",
                                  "../stomp/stomp_angular_correlation.cc",
                                  "../stomp/stomp_pixel.cc",
                                  "../stomp/stomp_scalar_pixel.cc",
                                  "../stomp/stomp_tree_pixel.cc",
                                  "../stomp/stomp_itree_pixel.cc",
                                  "../stomp/stomp_base_map.cc",
                                  "../stomp/stomp_map.cc",
                                  "../stomp/stomp_scalar_map.cc",
                                  "../stomp/stomp_tree_map.cc",
                                  "../stomp/stomp_itree_map.cc",
                                  "../stomp/stomp_geometry.cc",
                                  "../stomp/stomp_util.cc",
                                  "stomp_wrap.cxx"],
                         )

setup (name = "stomp",
       version = "0.1",
       author      = "Ryan Scranton",
       description = """
STOMP is a set of libraries for doing astrostatistical analysis on the
celestial sphere.  The goal is to enable descriptions of arbitrary regions
on the sky which may or may not encode futher spatial information (galaxy
density, CMB temperature, observational depth, etc.) and to do so in such
a way as to make the analysis of that data as algorithmically efficient as
possible.
""",
       data_files=[('ups',['ups/stomp.table'])],
       ext_modules = [stomp_module],
       py_modules = ["stomp"],
       include_dirs=include_dirs)
