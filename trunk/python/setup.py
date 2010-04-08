#!/usr/bin/env python

"""
setup.py file for STOMP library
"""

from distutils.core import setup, Extension


stomp_module = Extension("_stomp",
                           sources=["../stomp/stomp_core.cc",
                                    "../stomp/stomp_angular_bin.cc",
                                    "../stomp/stomp_angular_coordinate.cc",
                                    "../stomp/stomp_angular_correlation.cc",
                                    "../stomp/stomp_pixel.cc",
                                    "../stomp/stomp_scalar_pixel.cc",
                                    "../stomp/stomp_tree_pixel.cc",
                                    "../stomp/stomp_base_map.cc",
                                    "../stomp/stomp_map.cc",
                                    "../stomp/stomp_scalar_map.cc",
                                    "../stomp/stomp_tree_map.cc",
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
       ext_modules = [stomp_module],
       py_modules = ["stomp"],
       )
