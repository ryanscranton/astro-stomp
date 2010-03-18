// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the basic constants necessary to make the library
// function and should be included in any program using the library.


#ifndef STOMP_CORE_H
#define STOMP_CORE_H

#include <stdint.h>

namespace Stomp {

// Here we're going to define the basic constants and methods that will be
// used throughout all of the other library classes.

// First some trigonometric values.
extern const double Pi;
extern const double DegToRad;
extern const double RadToDeg;
extern const double StradToDeg;

// For historical reasons, coordinate system is built around the SDSS
// survey coordinates rather than traditional equatorial RA-DEC coordinates.
// To switch to those coordinates, the next five functions would need to be
// modified so that EtaOffSet, Node and EtaPole all return zero.
extern const double EtaOffSet;
extern const double SurveyCenterRA;
extern const double SurveyCenterDEC;
extern const double Node;
extern const double EtaPole;

// For the purposes of rapid localization, we set a basic level of
// pixelization that divides the sphere into 7488 superpixels (the value is
// chosen such that the width of one pixel matches the fiducial width of a
// stripe in the SDSS survey coordinates.
//
// Pixels are addressed hierarchically, so a given pixel is refined into 4
// sub-pixels and is joined with 3 nearby pixels to form a superpixel.  The
// level of refinement is encoded in the "resolution" level, with the
// following two functions defining the limits on acceptable values (basically
// limited by the number of pixels that can be addressed in a single
// superpixel with 32-bit integers).  The current limits allow for about
// half a terapixel on the sphere, which corresponds to roughly 2 arcsecond
// resolution.  Valid resolution values are all powers of 2; refining the
// pixel scale increases resolution by a factor of 2 at each level and
// coarsening the pixel scale reduces the resolution by a factor of 2 at each
// level.
extern const uint32_t Nx0;
extern const uint32_t Ny0;
extern const uint8_t HPixLevel;
extern const uint8_t MaxPixelLevel;
extern const uint32_t HPixResolution;
extern const uint32_t MaxPixelResolution;
extern const uint8_t ResolutionLevels;
extern const double HPixArea;
extern const uint32_t MaxPixnum;
extern const uint32_t MaxSuperpixnum;

// Some methods to deal with comparisons between doubles.
bool DoubleLT(double a, double b);
bool DoubleLE(double a, double b);
bool DoubleGT(double a, double b);
bool DoubleGE(double a, double b);
bool DoubleEQ(double a, double b);

// And finally, a method to let us translate between resolution and
// resolution level.
uint8_t MostSignificantBit(uint32_t input_int);

} // end namespace Stomp

#endif
