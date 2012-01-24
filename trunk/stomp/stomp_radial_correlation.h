// Copyright 2011  All Rights Reserved.
// Authors: morrison.chrisb@gmail.com (Chris Morrison), 
//          ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible. For the most part this class a takes much of its code from 
// AngularCoorelation with the cavet that the RadialBins no longer have a set
// angular projection or resolution on the sky.
//
// This header file contains the class for calculating projected radial 
// correlations on the sphere.  In general, different methods are more 
// efficient on small vs.large angular scales, so this class draws on nearly 
// the entire breadth of the STOMP library.

#ifndef STOMP_RADIAL_CORRELATION_H
#define STOMP_RADIAL_CORRELATION_H

#include <vector>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_angular_bin.h"
#include "stomp_radial_bin.h"

namespace Stomp {

class Map;        // class declaration in stomp_map.h
class ScalarMap;  // class declaration in stomp_scalar_map.h
class TreeMap;    // class declaration in stomp_tree_map.h
class AngularCorrelation;
class RadialCorrelation;

typedef std::vector<RadialCorrelation> RWThetaVector;
typedef RWThetaVector::iterator RWThetaIterator;

class RadialCorrelation : public AngularCorrelation {
  // Class object for calculating auto-correlations and cross-correlations
  // given a set of objects and a Map in projected physical distand. 
  // Broadly speaking, this is a container class for a set of RadialBin objects
  // which collectively span some range of Radial scales.  Accordingly, the 
  // methods are generall intended to package the machinery of the 
  // auto-correlation and cross-correlation calculations into simple, one-line 
  // calls.

  // The first constructor takes an radial minimum and maximum (in Mpc/h)
  // and constructs a logrithmic binning scheme using the specified number
  // of bins per decade (which can be a non-integer value, obviously).  The
  // bins are such that the minimum radial scale of the first bin will be
  // radius_min and the maximum radial scale of the last bin with be
  // radius_max. The last boolean argument controls whether or not an
  // pixel resolution will be assigned to the bins.  If it is false, then
  // the resolution values will all be -1.
 public:
  RadialCorrelation(double radius_min, double radius_max,
		    double bins_per_decade);

  // The alternate constructor is used for a linear binning scheme.  The
  // relationship between radius_min and radius_max remains the same and the
  // spacing of the bins is determined based on the requested number of bins.
  RadialCorrelation(uint32_t n_bins, double radius_min, double radius_max);

  ~RadialCorrelation() {
    radialbin_.clear();
  };

  // If we're going to use regions to find jack-knife errors, then we need
  // to initialize the AngularBins to handle this state of affairs or possibly
  // clear out previous calculations.
  void InitializeRegions(int16_t n_regions);
  void ClearRegions();
  int16_t NRegion();

  // Some wrapper methods for find the auto-correlation and cross-correlations
  void FindAutoCorrelation(Map& stomp_map,
			   CosmoVector& galaxy,
			   uint8_t random_iterations = 1);
  void FindAutoCorrelationWithRegions(Map& stomp_map,
				      CosmoVector& galaxy,
				      uint8_t random_iterations = 1,
				      uint16_t n_regions = 0);

  // In general, the code will use a pair-based method for small angular
  // scanes and a pixel-based method for large angular scales.  In the above
  // methods, this happens automatically.  If you want to run these processes
  // separately, these methods allow you to do this.  If the Map used
  // to call these methods has initialized regions, then the estimators will
  // use the region-based methods.
  //void FindPixelAutoCorrelation(Map& stomp_map, CosmoVector& galaxy);
  //void FindPixelAutoCorrelation(ScalarMap& stomp_map);
  //void FindPixelCrossCorrelation(Map& stomp_map, CosmoVector& galaxy_a,
  //				 CosmoVector& galaxy_b);
  //void FindPixelCrossCorrelation(ScalarMap& stomp_map_a,
  //				 ScalarMap& stomp_map_b);
  void FindPairAutoCorrelation(Map& stomp_map, CosmoVector& galaxy,
			       uint8_t random_iterations = 1);
  void FindPairCrossCorrelation(Map& stomp_map,
				CosmoVector& galaxy_a,
				CosmoVector& galaxy_b,
				uint8_t random_iterations = 1);

  // Once we're done calculating our correlation function, we can write it out
  // to an ASCII file.  The output format will be
  //
  //   RADIUS_P   W(RADIUS_P)   dW(RADIUS_P)
  //
  // where RADIUS_P is the radial scale in Mpc/h and dW(RADIUS_P) is the jack-knife
  // error based on regionating the data.  If the radial correlation has been
  // calculated without regions, then this column will be omitted.
  bool Write(const std::string& output_file_name);

  double Covariance(uint8_t bin_idx_a, uint8_t bin_idx_b);
  bool WriteCovariance(const std::string& output_file_name);

  // In some cases, we want to default to using either the pair-based or
  // pixel-based estimator for all of our bins, regardless of angular scale.
  // These methods allow us to over-ride the default behavior of the
  // correlation code.
  void UseOnlyPairs();

  // Now, some accessor methods for finding the angular range of the bins
  // with a given resolution attached to them (the default value returns the
  // results for all angular bins; for pair-based bins, resolution = -1).
  double RadiusMin(uint32_t resolution = 1);
  double RadiusMax(uint32_t resolution = 1);
  RadialIterator Begin(uint32_t resolution = 1);
  RadialIterator End(uint32_t resolution = 1);
  RadialIterator BinIterator(uint8_t bin_idx = 0);
  uint32_t NBins();

 private:
  RadialVector radialbin_;
  RadialIterator radial_pixel_begin_, radial_pixel_end_;
  RadialIterator radial_pair_begin_, radial_pair_end_;
  double r_min_, r_max_;
  uint32_t min_resolution_, max_resolution_, regionation_resolution_;
  int16_t n_region_;
  bool manual_resolution_break_;
};

} // end namespace Stomp

#endif
