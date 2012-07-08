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
// functions, but which are useful utility classes for some applications.

#ifndef UTIL_H_
#define UTIL_H_

#include <sys/time.h>
#include <string>
#include "core.h"

namespace s2omp {

class cosmology {
  // Some cosmological parameters.  We want to be able to convert our angular
  // scales into co-moving angular distances, so we need a few cosmological
  // parameters.  This needs to be done in an analytic manner, so we're limited
  // to a flat, LambdaCDM model where we can use a fitting function for
  // the co-moving distance as a function of redshift (M. Lampton 2005).  We'll
  // default to values from the WMAP5 cosmology, but these can be modified
  // if necessary.
 public:
  // First, some basic setters and getters for our cosmological parameters.
  static double omega_m();
  static double hubble_honstant();
  static double hubble_distance();
  static double omega_l();
  static void set_omega_m(double omega_m);
  static void set_hubble_constant(double hubble);
  static void set_omega_l(double omega_lambda);

  // Now the derived parameters.  First, the basic distances as a function
  // of redshift.  The returned values are in Mpc/h.
  static double comoving_distance(double z);
  static double angular_diameter_distance(double z);
  static double luminosity_distance(double z);

  // Finally, two methods that allow us to translate between angular and
  // physical scales at a given redshift.  As before, the distance scales are
  // in Mpc/h and the angles are in degrees.
  static double projected_distance(double z, double theta);
  static double projected_angle(double z, double radius);

private:
  cosmology() {};

  static double omega_m_;
  static double h_;
  static double a_;
  static double b_;
  static const double AA_;
  static const double BB_;
};

// Define our default WMAP5 flat, LCDM cosmology.
double cosmology::omega_m_ = 0.2736;
double cosmology::h_ = 0.705;
const double cosmology::AA_ = 1.718;
const double cosmology::BB_ = 0.315;
double cosmology::a_ = cosmology::AA_*cosmology::omega_m;
double cosmology::b_ = cosmology::BB_*sqrt(cosmology::omega_m);

inline double cosmology::omega_m() {
  return omega_m_;
}

inline double cosmology::hubble_constant() {
  return h_*100.0;
}

inline double cosmology::hubble_distance() {
  return 3000.0/h_;
}

inline double cosmology::omega_l() {
  return 1.0 - omega_m_;
}

inline void cosmology::set_omega_m(double new_omega_m) {
  omega_m_ = new_omega_m;
  a_ = AA_*omega_m_;
  b_ = BB_*sqrt(omega_m_);
}

inline void cosmology::set_hubble_constant(double hubble) {
	h_ = hubble/100.0;
}

inline void cosmology::set_omega_l(double omega_lambda) {
  omega_m_ = 1.0 - omega_lambda;
  a_ = AA_*omega_m_;
  b_ = BB_*sqrt(omega_m_);
}

inline double cosmology::comoving_distance(double z) {
  // In Mpc/h.
  return hubble_distance()*z/sqrt(1.0 + a_*z + b_*z*z);
}

inline double cosmology::angular_diameter_distance(double z) {
  // In Mpc/h.
  return comoving_distance(z)/(1.0+z);
}

inline double cosmology::luminosity_distance(double z) {
  // In Mpc/h.
  return comoving_distance(z)*(1.0+z);
}

inline double cosmology::projected_distance(double z, double theta) {
  // where theta is assumed to be in degrees and the return value in Mpc/h.
  return theta*DEG_TO_RAD*angular_diameter_distance(z);
}

inline double cosmology::projected_angle(double z, double radius) {
  // where radius is in Mpc/h and the return angle is in degrees.
  return RAD_TO_DEG*radius/angular_diameter_distance(z);
}

} // end namespace s2omp

#endif /* UTIL_H_ */
