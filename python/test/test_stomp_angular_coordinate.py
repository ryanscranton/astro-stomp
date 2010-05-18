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

This module tests the wrapping of the AngularCoordinate class and its
derived classes.
"""

__author__ = "Ryan Scranton (ryan.scranton@gmail.com)"
__copyright__ = "Copyright 2010, Ryan Scranton"
__license__ = "BSD"
__version__ = "1.0"

import stomp
import math
import unittest

class TestStompAngularCoordinate(unittest.TestCase):

    """
    Unit testing class for the Stomp::AngularCoordinate class and derivatives.
    """

    def testAngularCoordinateBasic(self):
        """Test the basic constructors and coordinate transforms."""
        ang = stomp.AngularCoordinate()
        ang.SetSurveyCoordinates(20.0, 60.0)

        tmp_ang = stomp.AngularCoordinate(20.0, 60.0,
                                          stomp.AngularCoordinate.Survey)

        self.assertAlmostEqual(ang.RA(), tmp_ang.RA())
        self.assertAlmostEqual(ang.DEC(), tmp_ang.DEC())
        self.assertAlmostEqual(ang.GalLon(), tmp_ang.GalLon())
        self.assertAlmostEqual(ang.GalLat(), tmp_ang.GalLat())

    def testAngularCoordinateVectorMath(self):
        """Test our dot & cross product methods."""
        ang_a = stomp.AngularCoordinate(0.0, 0.0,
                                        stomp.AngularCoordinate.Survey)
        ang_b = stomp.AngularCoordinate(45.0, 0.0,
                                        stomp.AngularCoordinate.Survey)
        ang_c = stomp.AngularCoordinate(60.0, 0.01,
                                        stomp.AngularCoordinate.Survey)
        ang_d = stomp.AngularCoordinate(0.0, 60.0,
                                        stomp.AngularCoordinate.Survey)

        self.assertAlmostEqual(ang_a.DotProduct(ang_b), ang_b.DotProduct(ang_a))
        # TODO(ryan.scranton) Need some clean cross product tests.

    def testAngularCoordinatePositionAngle(self):
        """Test our position angle methods."""
        ang_a = stomp.AngularCoordinate(0.0, 0.0,
                                        stomp.AngularCoordinate.Equatorial)
        ang_b = stomp.AngularCoordinate(1.0, 0.0,
                                        stomp.AngularCoordinate.Equatorial)
        self.assertAlmostEqual(ang_a.PositionAngle(ang_b), 90.0)
        self.assertAlmostEqual(ang_b.PositionAngle(ang_a), 270.0)

        ang_b.SetEquatorialCoordinates(0.0, 1.0)
        self.assertAlmostEqual(ang_a.PositionAngle(ang_b), 0.0)
        self.assertAlmostEqual(ang_b.PositionAngle(ang_a), 180.0)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(
        TestStompAngularCoordinate)
    unittest.TextTestRunner(verbosity=2).run(suite)

