// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file loads the header files for the various pieces of the STOMP
// library.  This is intended to be a convenience so that users can have access
// to the library objects without having to specify which pieces of the library
// they want to access explicity.

#ifndef STOMP_H
#define STOMP_H

#include <stomp/stomp_core.h>
#include <stomp/stomp_angular_coordinate.h>
#include <stomp/stomp_angular_bin.h>
#include <stomp/stomp_angular_correlation.h>
#include <stomp/stomp_pixel.h>
#include <stomp/stomp_scalar_pixel.h>
#include <stomp/stomp_tree_pixel.h>
#include <stomp/stomp_base_map.h>
#include <stomp/stomp_map.h>
#include <stomp/stomp_scalar_map.h>
#include <stomp/stomp_tree_map.h>
#include <stomp/stomp_geometry.h>
#include <stomp/stomp_util.h>

#endif
