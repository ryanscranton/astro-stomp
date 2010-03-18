// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the class for calculating angular correlations on
// the sphere.  In general, different methods are more efficient on small vs.
// large angular scales, so this class draws on nearly the entire breadth of
// the STOMP library.

#ifndef STOMP_ANGULAR_CORRELATION_H
#define STOMP_ANGULAR_CORRELATION_H

#include <vector>
#include "stomp_angular_coordinate.h"
#include "stomp_angular_bin.h"

namespace Stomp {

class Map;        // class declaration in stomp_map.h
class ScalarMap;  // class declaration in stomp_scalar_map.h
class TreeMap;    // class declaration in stomp_tree_map.h
class AngularCorrelation;

typedef std::vector<AngularCorrelation> WThetaVector;
typedef WThetaVector::iterator WThetaIterator;

class AngularCorrelation {
  // Class object for calculating auto-correlations and cross-correlations
  // given a set of objects and a Map.  Broadly speaking, this is a
  // container class for a set of AngularBin objects which collectively
  // span some range of angular scales.  Accordingly, the methods are generally
  // intended to package the machinery of the auto-correlation and
  // cross-correlation calculations into simple, one-line calls.

 public:
  // The first constructor takes an angular minimum and maximum (in degrees)
  // and constructs a logrithmic binning scheme using the specified number
  // of bins per decade (which can be a non-integer value, obviously).  The
  // bins are such that the minimum angular scale of the first bin will be
  // theta_min and the maximum angular scale of the last bin with be
  // theta_max.  The last boolean argument controls whether or not an
  // pixel resolution will be assigned to the bins.  If it is false, then
  // the resolution values will all be -1.
  AngularCorrelation(double theta_min, double theta_max,
		     double bins_per_decade, bool assign_resolutions = true);

  // The alternate constructor is used for a linear binning scheme.  The
  // relationship between theta_min and theta_max remains the same and the
  // spacing of the bins is determined based on the requested number of bins.
  AngularCorrelation(uint32_t n_bins, double theta_min, double theta_max,
		     bool assign_resolutions = true);
  ~AngularCorrelation() {
    thetabin_.clear();
  };

  // Find the resolution we would use to calculate correlation functions for
  // each of the bins.  If this method is not called, then the resolution
  // for each bin is set to -1, which would indicate that any correlation
  // calculation with that bin should be done using a pair-based estimator.
  void AssignBinResolutions(double lammin = -70.0, double lammax = 70.0,
			    uint32_t min_resolution = MaxPixelResolution);

  // For small angular scales, it's usually faster and more memory
  // efficient to use a pair-based estimator.  To set this scale, we choose
  // a maximum resolution scale we're willing to use our pixel-based estimator
  // on and modify all smaller angular bins to use the pair-based estimator.
  void SetMaxResolution(uint32_t resolution);

  // Additionally, if we are using regions to calculate correlation functions,
  // we need to set the minimum resolution to match the resolution used to
  // divide the total survey area.
  void SetMinResolution(uint32_t resolution);

  // Some wrapper methods for find the auto-correlation and cross-correlations
  void FindAutoCorrelation(Map& stomp_map,
			   WAngularVector& galaxy,
			   uint8_t random_iterations = 1);
  void FindCrossCorrelation(Map& stomp_map,
			    WAngularVector& galaxy_a,
			    WAngularVector& galaxy_b,
			    uint8_t random_iterations = 1);

  // Variation on the wrapper methods that use regions to calculate the
  // cosmic variance on the correlation functions.  If you don't specify the
  // number of regions to use, the code will default to twice the number of
  // angular bins.
  void FindAutoCorrelationWithRegions(Map& stomp_map,
				      WAngularVector& galaxy,
				      uint8_t random_iterations = 1,
				      uint16_t n_regions = 0);
  void FindCrossCorrelationWithRegions(Map& stomp_map,
				       WAngularVector& galaxy_a,
				       WAngularVector& galaxy_b,
				       uint8_t random_iterations = 1,
				       uint16_t n_regions = 0);

  // In general, the code will use a pair-based method for small angular
  // scanes and a pixel-based method for large angular scales.  In the above
  // methods, this happens automatically.  If you want to run these processes
  // separately, these methods allow you to do this.  If the Map used
  // to call these methods has initialized regions, then the estimators will
  // use the region-based methods.
  void FindPixelAutoCorrelation(Map& stomp_map, WAngularVector& galaxy);
  void FindPixelAutoCorrelation(ScalarMap& stomp_map);
  void FindPixelCrossCorrelation(Map& stomp_map, WAngularVector& galaxy_a,
				 WAngularVector& galaxy_b);
  void FindPixelCrossCorrelation(ScalarMap& stomp_map_a,
				 ScalarMap& stomp_map_b);
  void FindPairAutoCorrelation(Map& stomp_map, WAngularVector& galaxy,
			       uint8_t random_iterations = 1);
  void FindPairCrossCorrelation(Map& stomp_map,
				WAngularVector& galaxy_a,
				WAngularVector& galaxy_b,
				uint8_t random_iterations = 1);

  // Now, some accessor methods for finding the angular range of the bins
  // with a given resolution attached to them (the default value returns the
  // results for all angular bins; for pair-based bins, resolution = -1).
  double ThetaMin(uint32_t resolution = 1);
  double ThetaMax(uint32_t resolution = 1);
  double Sin2ThetaMin(uint32_t resolution = 1);
  double Sin2ThetaMax(uint32_t resolution = 1);
  ThetaIterator Begin(uint32_t resolution = 1);
  ThetaIterator End(uint32_t resolution = 1);
  ThetaIterator Find(ThetaIterator begin, ThetaIterator end,
		     double sin2theta);
  ThetaIterator BinIterator(uint8_t bin_idx = 0);
  uint32_t NBins();
  uint32_t MinResolution();
  uint32_t MaxResolution();

 private:
  ThetaVector thetabin_;
  ThetaIterator theta_pixel_begin_;
  ThetaIterator theta_pair_begin_, theta_pair_end_;
  double theta_min_, theta_max_, sin2theta_min_, sin2theta_max_;
  uint32_t min_resolution_, max_resolution_;
};

} // end namespace Stomp

#endif
