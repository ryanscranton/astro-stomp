// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains a single class for dealing with projected radial
// bins on the sky.  It's functionally similar to the AngularBin class, but
// the angular limits are determined by the combination of radial scale and
// redshift.

#ifndef STOMP_RADIAL_BIN_H
#define STOMP_RADIAL_BIN_H

#include <vector>
#include <math.h>

namespace Stomp {

class AngularBin;
class RadialBin;

typedef std::vector<RadialBin> RadialVector;
typedef RadialVector::iterator RadialIterator;
typedef std::pair<RadialIterator, RadialIterator> RadialPair;
typedef std::vector<RadialBin *> RadialPtrVector;
typedef RadialPtrVector::iterator RadialPtrIterator;

class RadialBin : public AngularBin {
  // Class object for holding the data associated with a single projected
  // radial annulus.  Each instance of the class contains a lower and upper
  // limit that defines the annulus as well as methods for testing against
  // those limits and data fields that are used for calculating angular
  // auto-correlations and cross-correlations via the AngularCorrelation
  // class described in stomp_angular_correlation.h.

 public:
  RadialBin();
  virtual ~RadialBin();

  // The simplest RadialBin object needs only a minimum and maximum radius
  // range (generally denoted by r_p in the literature) and a redshift, which
  // is necessary to translate the projected radii into angular scales.  The
  // assumed units for the radius are Mpc/h, as in the Cosmology class.
  RadialBin(double r_min, double r_max, double redshift);

  // A common method for calculating the error on angular correlations is to
  // divide the survey area up into equal area regions and use jack-knife
  // methods to estimate the variance on the correlation function.  This
  // constructor sets up the RadialBin object for that sort of operation.
  RadialBin(double r_min, double r_max, double redshift, int16_t n_regions);

  // Seting the radial minima, maxima and mid-point.
  //
  // Depending on whether we're using linear or logarithmic angular binning,
  // we'll need to set the mid-point of the bin by hand.  We also need to
  // over-ride the methods to set the angular limits of the bin so that we can
  // ensure that all of our bounds stay in sync.
  void SetRadiusMin(double r_min);
  void SetRadiusMax(double r_max);
  void SetRadius(double r);
  void SetRedshift(double z);

  // Some simple methods to complement those in AngularBin
  bool WithinRadialBounds(double r);
  double Radius();
  double RadiusMin();
  double RadiusMax();
  double Redshift();

  // Finally, some static methods which the AngularCorrelation method will use
  // to order its vectors of RadialBin objects.
  static bool RadialOrder(RadialBin r_a, RadialBin r_b);
  static bool ReverseResolutionOrder(RadialBin r_a,
				     RadialBin r_b);

 private:
  double r_min_, r_max_, redshift_, r_;
};

} // end namespace Stomp

#endif
