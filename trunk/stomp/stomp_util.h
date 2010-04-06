// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains classes which are extraneous to the core library
// functions, but which are useful utility classes to for some applications.

#ifndef STOMP_UTIL_H
#define STOMP_UTIL_H

#include <sys/time.h>

namespace Stomp {

class Cosmology;
class StompWatch;

class Cosmology {
  // Some cosmological parameters.  We want to be able to convert our angular
  // scales into co-moving angular distances, so we need a few cosmological
  // parameters.  This needs to be done in an analytic manner, so we're limited
  // to a flat, LambdaCDM model where we can use a fitting function for
  // the co-moving distance as a function of redshift (M. Lampton 2005).  We'll
  // default to values from the WMAP5 cosmology, but these can be modified
  // if necessary.
 public:
  static double omega_m;
  static double h;
  static double a_;
  static double b_;
  static const double AA_;
  static const double BB_;

  // First, some basic setters and getters for our cosmological parameters.
  static double OmegaM();
  static double HubbleConstant();
  static double HubbleDistance();
  static double OmegaL();
  static void SetOmegaM(double omega_m);
  static void SetHubbleConstant(double hubble);
  static void SetOmegaL(double omega_lambda);

  // Now the derived parameters.  First, the basic distances as a function
  // of redshift.  The returned values are in Mpc/h.
  static double ComovingDistance(double z);
  static double AngularDiameterDistance(double z);
  static double LuminosityDistance(double z);

  // Finally, two methods that allow us to translate between angular and
  // physical scales at a given redshift.  As before, the distance scales are
  // in Mpc/h and the angles are in degrees.
  static double ProjectedDistance(double z, double theta);
  static double ProjectedAngle(double z, double radius);
};

class StompWatch {
  // A simple class to handle time keepig inside of a program.  The return
  // values are in seconds, so the primary goal here is not high precision so
  // much as providing a tool for tuning algorithm performance on a fixed data
  // set.

 public:
  StompWatch();

  // Methods to start and stop the internal timer in the class.
  void StartTimer();
  void StopTimer();

  // Return the elapsed time (in seconds) between the last StartTimer() call
  // and the latest StopTimer() call.
  double ElapsedTime();

 private:
  timeval start;
  timeval stop;
};

// A simple utility function for breaking up a string into its component items,
// generally when the string is a list of filenames or such (separated by the
// input delimiter).
void Tokenize(const std::string& str,
	       std::vector<std::string>& tokens,
	       const std::string& delimiters = " ");

} // end namespace Stomp

#endif
