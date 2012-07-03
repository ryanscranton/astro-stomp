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

#include <stdint.h>

namespace s2omp {

// Here we're going to define the basic constants and methods that will be
// used throughout all of the other library classes.

// First some trigonometric values.
extern const double PI;
extern const double DEG_TO_RAD;
extern const double RAD_TO_DEG;
extern const double STRAD_TO_DEG2;

extern const int MAX_LEVEL;

bool DoubleLT(double a, double b);
bool DoubleLE(double a, double b);
bool DoubleGT(double a, double b);
bool DoubleGE(double a, double b);
bool DoubleEQ(double a, double b);

} // end namespace s2omp

#endif /* CORE_H_ */
