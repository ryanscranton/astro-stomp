#!/usr/bin/env python

# Unit testing module for STOMP core class
# Copyright (c) 2010, Ryan Scranton
# 
# All rights reserved.

"""
STOMP is a set of libraries for doing astrostatistical analysis on the
celestial sphere.  The goal is to enable descriptions of arbitrary regions
on the sky which may or may not encode futher spatial information (galaxy
density, CMB temperature, observational depth, etc.) and to do so in such
a way as to make the analysis of that data as algorithmically efficient as
possible.

This module tests the basic constants necessary to make the library function.
"""

__author__ = "Ryan Scranton (ryan.scranton@gmail.com)"
__copyright__ = "Copyright 2010, Ryan Scranton"
__license__ = "BSD"
__version__ = "1.0"

import stomp
import math
import unittest

class TestStompCore(unittest.TestCase):

    """
    Unit testing class for core STOMP variables.
    """

    def testConstants(self):
        self.assertEqual(stomp.HPixResolution, 4)
        self.assertEqual(stomp.MaxPixelResolution, 32768)
        self.assertEqual(stomp.HPixLevel, 2)
        self.assertEqual(stomp.MaxPixelLevel, 15)
        self.assertEqual(stomp.MaxSuperpixnum, 7488)
        self.assertAlmostEqual(stomp.Pi, math.pi)
        self.assertAlmostEqual(stomp.DegToRad, math.pi/180.0)
        self.assertAlmostEqual(stomp.RadToDeg, 180.0/math.pi)

    def testMostSignifcantBit(self):
        self.assertEqual(stomp.MostSignificantBit(4), 2)       
        self.assertEqual(stomp.MostSignificantBit(32768), 15)
        self.assertEqual(stomp.MostSignificantBit(32769), 15)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStompCore)
    unittest.TextTestRunner(verbosity=2).run(suite)

