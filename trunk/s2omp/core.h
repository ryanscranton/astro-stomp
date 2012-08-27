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
#include <stdint.h>
#include <math.h>
#include <vector>
#include <base/integral_types.h>

namespace s2omp {

// Here we're going to define the basic constants and methods that will be
// used throughout all of the other library classes.

// As these are the two core classes of the library we define their existence
// now as well vectors of these quantities.
class pixel;
class point;

typedef std::vector<pixel> pixel_vector;
typedef pixel_vector::iterator pixel_iterator;

typedef std::vector<point> point_vector;
typedef point_vector::iterator point_iterator;

// First some trigonometric values.
extern const double PI = 2.0*asin(1.0);
extern const double DEG_TO_RAD = PI/180.0;
extern const double RAD_TO_DEG = 180.0/PI;
extern const double STRAD_TO_DEG2 = 180.0*180.0/(PI*PI);

extern const int MAX_LEVEL = 30;

inline bool double_lt(double a, double b) {
	return a < b - 1.0e-15 ? true : false;
}

inline bool double_le(double a, double b) {
	return a <= b - 1.0e-15 ? true : false;
}

inline bool double_gt(double a, double b) {
	return a > b + 1.0e-15 ? true : false;
}

inline bool double_ge(double a, double b) {
  return a >= b - 1.0e-15 ? true : false;
}

inline bool double_eq(double a, double b) {
	return double_le(a, b) && double_ge(a, b) ? true : false;
}

} // end namespace s2omp

#endif /* CORE_H_ */
