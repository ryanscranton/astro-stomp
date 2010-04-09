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

This module tests the wrapping of the Pixel class.
"""

__author__ = "Ryan Scranton (ryan.scranton@gmail.com)"
__copyright__ = "Copyright 2010, Ryan Scranton"
__license__ = "BSD"
__version__ = "1.0"

import stomp
import math
import unittest

class TestStompPixel(unittest.TestCase):

    """
    Unit testing class for the Stomp::Pixel class.
    """

    def testPixelBasic(self):
        """Test the basic constructors and hierarchy."""
        ang = stomp.AngularCoordinate(60.0, 20.0,
                                      stomp.AngularCoordinate.Survey)
        pix = stomp.Pixel(ang, stomp.HPixResolution)

        for level in range(stomp.HPixLevel, stomp.MaxPixelLevel+1):
            resolution = stomp.Pixel.LevelToResolution(level)
            tmp_pix = stomp.Pixel(ang, resolution)

            self.assertEqual(pix.Pixnum(), tmp_pix.Superpixnum())
            self.assertEqual(pix.Superpixnum(), tmp_pix.Superpixnum())


    def testPixelAnnulusIntersection(self):
        """Test our method for annulus intersection checking."""
        # Initialize a pixel at a given location
        ang = stomp.AngularCoordinate(60.0, 20.0,
                                      stomp.AngularCoordinate.Survey)
        pix = stomp.Pixel(ang, stomp.HPixResolution)
        
        # now find the rough angular scale of the pixel
        pix_center = stomp.AngularCoordinate()
        pix.Ang(pix_center)
        ang.SetSurveyCoordinates(pix.Lambda(), pix.EtaMin())
        pixel_radius = pix_center.AngularDistance(ang)

        # Start with an annulus that definitely contains the pixel
        self.assertEqual(pix.IntersectsAnnulus(pix_center, 0.0,
                                               10.0*pixel_radius), 1)

        # Now shrink to a radius that is inside, but doesn't contain the pixel
        self.assertEqual(pix.IntersectsAnnulus(pix_center, 0.0,
                                               0.1*pixel_radius), -1)

        # Inscribed, but make it an annulus
        self.assertEqual(pix.IntersectsAnnulus(pix_center, 0.8*pixel_radius,
                                               1.25*pixel_radius), -1)

        # Circle outside the pixel
        ang.SetSurveyCoordinates(pix.Lambda(), pix.Eta()+10.0*pixel_radius)
        self.assertEqual(pix.IntersectsAnnulus(ang, 0.0, 0.1*pixel_radius), 0)

        # Annulus outside the pixel
        ang.SetSurveyCoordinates(pix.Lambda(), pix.Eta()+10.0*pixel_radius)
        self.assertEqual(pix.IntersectsAnnulus(ang, 0.1*pixel_radius,
                                               0.2*pixel_radius), 0)

        # Annulus with center outside, but should contain the pixel
        ang.SetSurveyCoordinates(pix.Lambda(), pix.Eta()+10.0*pixel_radius)
        self.assertEqual(pix.IntersectsAnnulus(ang, 0.1*pixel_radius,
                                               20.0*pixel_radius), 1)

        # Annulus with center outside, but should intersect the pixel
        ang.SetSurveyCoordinates(pix.Lambda(), pix.Eta()+2.0*pixel_radius)
        self.assertEqual(pix.IntersectsAnnulus(ang, 1.25*pixel_radius,
                                               2.5*pixel_radius), -1)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStompPixel)
    unittest.TextTestRunner(verbosity=2).run(suite)

