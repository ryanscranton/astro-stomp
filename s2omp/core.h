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


#ifndef CORE_H_
#define CORE_H_

#include <algorithm>
#include <ctime>
#include <math.h>
#include <stdint.h>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include <s2/s2.h>
#include <base/integral_types.h>

namespace s2omp {

// Here we're going to define the basic constants and methods that will be
// used throughout all of the other library classes.

// As these are the two core classes of the library we define their existence
// now as well vectors of these quantities.
class pixel;
class point;

class bound_interface;
class circle_bound;
class annulus_bound;

class angular_bin;

class coverer;
class region_map;

typedef std::vector<pixel> pixel_vector;
typedef pixel_vector::const_iterator pixel_iterator;

typedef std::set<pixel> pixel_set;
typedef pixel_set::const_iterator pixel_set_iterator;

typedef std::vector<point> point_vector;
typedef point_vector::const_iterator point_iterator;

typedef std::vector<circle_bound *> circle_ptr_vector;
typedef circle_ptr_vector::const_iterator circle_ptr_iterator;

typedef std::vector<angular_bin> theta_vector;
typedef theta_vector::iterator theta_iterator;
typedef theta_vector::const_iterator theta_const_iterator;
typedef std::pair<theta_iterator, theta_iterator> theta_pair;
typedef std::vector<angular_bin *> theta_ptr_vector;
typedef theta_ptr_vector::const_iterator theta_ptr_iterator;

// First some trigonometric values.
static double const PI = 2.0 * asin(1.0);
static double const DEG_TO_RAD = PI / 180.0;
static double const RAD_TO_DEG = 180.0 / PI;
static double const STRAD_TO_DEG2 = 180.0 * 180.0 / (PI * PI);

// Angle of obliquity for converting Equatorial to Ecliptic coordinates:
// 23° 26′ 21.406″ from JPL (2000) in J2000.
static double const OBLIQUITY_RAD = 23.439279444444445 * DEG_TO_RAD;
static double const COS_OBLIQUITY = cos(OBLIQUITY_RAD);
static double const SIN_OBLIQUITY = sin(OBLIQUITY_RAD);

// Defaults from S2 for the maximum pixel level and the default number of
// pixels in an exterior covering.
static int const MAX_LEVEL = S2::kMaxCellLevel;
static long const DEFAULT_COVERING_PIXELS = 8;

// It can be useful to define a default resolution level that should be useful
// for a rough pixelization.  Level = 6 corresponds to roughly 24.5k pixels on
// the sky with an average size of about 1.27 degrees on a side.
static int const DEFAULT_LEVEL = 6;

// For the tree_pixel and tree_union classes, it's also useful to define a
// default carrying capacity for each node in the tree structure.
static uint const DEFAULT_NODE_CAPACITY = 200;

// From S2: Multiply a positive number by this constant to ensure that the
// result of a floating point operation is at least as large as the true
// infinite-precision result.
static double const FLOAT_ROUND_UP = 1.0 + 1.0 / (uint64(1) << 52);

// Default region value to return if we can't find an input pixel or point in
// a region_map.
static int const INVALID_REGION_VALUE = -1;

// For generating random numbers, we define a few static classes.
static std::mt19937_64 MT_GENERATOR(std::time(0));
static std::uniform_int_distribution<uint64> UNIFORM_UINT64;
static std::uniform_real_distribution<double> UNIFORM_DOUBLE(0.0, 1.0);

inline bool double_lt(double a, double b) {
  return a < b * FLOAT_ROUND_UP;
}

inline bool double_le(double a, double b) {
  return a <= b * FLOAT_ROUND_UP;
}

inline bool double_gt(double a, double b) {
  return a > b / FLOAT_ROUND_UP;
}

inline bool double_ge(double a, double b) {
  return a >= b / FLOAT_ROUND_UP;
}

inline bool double_eq(double a, double b) {
  return double_le(a, b) && double_ge(a, b);
}

inline int random_int(int max_value) {
  return static_cast<int>(max_value * UNIFORM_DOUBLE(MT_GENERATOR));
}

inline double random_double() {
  return UNIFORM_DOUBLE(MT_GENERATOR);
}

inline uint64 random_uint64() {
  return UNIFORM_UINT64(MT_GENERATOR);
}

inline uint64 random_uint64(uint64 max_value) {
  return static_cast<uint64>(max_value * UNIFORM_DOUBLE(MT_GENERATOR));
}

} // end namespace s2omp

#endif /* CORE_H_ */
