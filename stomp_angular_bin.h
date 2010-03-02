// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains a single class for dealing with angular annuli
// on the sky.  There is no heavy coding here (indeed, all of the class code is
// in this header), but the functionality here is necessary for the angular
// correlation operations in stomp_correlation.h as well as the pair finding
// in the TreePixel and TreeMap classes.

#ifndef STOMP_ANGULAR_BIN_H
#define STOMP_ANGULAR_BIN_H

#include <vector>
#include <math.h>

namespace Stomp {

class AngularBin;

typedef std::vector<AngularBin> ThetaVector;
typedef ThetaVector::iterator ThetaIterator;
typedef std::pair<ThetaIterator,ThetaIterator> ThetaPair;
typedef std::vector<AngularBin *> ThetaPtrVector;
typedef ThetaPtrVector::iterator ThetaPtrIterator;

class AngularBin {
  // Class object for holding the data associated with a single angular
  // annulus.  Each instance of the class contains a lower and upper angular
  // limit that defines the annulus as well as methods for testing against
  // those limits and data fields that are used for calculating angular
  // auto-correlations and cross-correlations via the AngularCorrelation
  // class described in stomp_angular_correlation.h.

 public:
  AngularBin();
  ~AngularBin();

  // The simplest AngularBin object needs only a minimum and maximum angular
  // range (generally denoted by _theta_ in the literature).  Theta is taken to
  // be in degrees.
  AngularBin(double theta_min, double theta_max);

  // A common method for calculating the error on angular correlations is to
  // divide the survey area up into equal area regions and use jack-knife
  // methods to estimate the variance on the correlation function.  This
  // constructor sets up the AngularBin object for that sort of operation.
  AngularBin(double theta_min, double theta_max, int16_t n_regions);

  // These two methods do the work of clearing and initializing the variables
  // used to store the region-related data.
  void ClearRegions();
  void InitializeRegions(int16_t n_regions);

  // There are two different methods for calculating the angular correlation
  // function, w(theta).  One is based on counting pairs separated by a given
  // angular distance.  The other pixelizes the survey area and sums the product
  // of over-densities for pixels separated by a given angular distance.  To
  // maximize the efficiency of the latter method, the resolution of the pixel
  // map needs to be matched to the angular scale of interest.  By storing this
  // resolution in the AngularBin, we can tell the AngularCorrelation object
  // which maps to use for this bin.  Alternatively, by setting this value to
  // an illegal resolution, we can signal that this angular scale is better
  // suited to the pair-based method.
  void SetResolution(uint16_t resolution);

  // Seting the angular minima, maxima and mid-point.
  //
  // Depending on whether we're using linear or logarithmic angular binning,
  // we'll need to set the mid-point of the angular bin by hand.
  void SetThetaMin(double theta_min);
  void SetThetaMax(double theta_max);
  void SetTheta(double theta);

  // A set of methods for determining if a given angular scale is within our
  // bounds.  Depending on the application, this scale may be most easily
  // given to us in terms of theta, sin(theta) or cos(theta).
  bool WithinBounds(double theta);
  bool WithinSin2Bounds(double sin2theta);
  bool WithinCosBounds(double costheta);

  // For the pixel-based w(theta), we use two internal variables:
  //
  //  * PixelWtheta, which stores the sum of the products of the overdensities
  //  * PixelWeight, which stores the number of such pixel pairs.
  //
  // w(theta) is then the ratio of these two numbers.
  void AddToPixelWtheta(double dwtheta, double dweight,
			int16_t region_a = -1, int16_t region_b = -1);

  // For the pair-counting, we use the methods in the TreePixel and TreeMap
  // classes.  Those methods are oblivious to the particular data sets they
  // are operating on, so they store the values in Weight (for the sum of the
  // products of the object Weights) and Counter, which stores the raw number
  // of point pairs.
  void AddToWeight(double weight, int16_t region = -1);
  void AddToCounter(uint32_t step=1, int16_t region = -1);

  // For calculating the pair-based w(theta), we use the Landy-Szalay estimator.
  // In the general case of a cross-correlation between two galaxy data sets,
  // there are four terms:
  //
  //  * galaxy-galaxy: the number of pairs between the two data sets.
  //  * random-random: the number of pairs between a randomized version of each
  //                   data set.
  //  * galaxy-random: the number of pairs between the first data set and a
  //                   randomized version of the second data set.
  //  * random-galaxy: the complement of galaxy-random.
  //
  // In the case of an autocorrelation, the last two terms are identical.  Once
  // we have the values of one of these combinations in the Weight and Counter
  // values, we can shift those values to the appropriate internal variable.
  void MoveWeightToGalGal();
  void MoveWeightToGalRand(bool move_to_rand_gal = false);
  void MoveWeightToRandGal(bool move_to_gal_rand = false);
  void MoveWeightToRandRand();

  // If the number of random points is not equal to the number of data points,
  // we will need to rescale the number of pairs accordingly.
  void RescaleGalGal(double pair_ratio);
  void RescaleGalRand(double pair_ratio);
  void RescaleRandGal(double pair_ratio);
  void RescaleRandRand(double pair_ratio);

  // Methods for reseting some or all of our internal data.
  void Reset();
  void ResetPixelWtheta();
  void ResetWeight();
  void ResetCounter();
  void ResetGalGal();
  void ResetGalRand();
  void ResetRandGal();
  void ResetRandRand();

  // Some basic getters for the internal data.
  uint16_t Resolution();
  int16_t NRegion();
  double Theta();
  double ThetaMin();
  double ThetaMax();
  double Sin2ThetaMin();
  double Sin2ThetaMax();
  double CosThetaMin();
  double CosThetaMax();

  // For the angular correlation values, we can access the value for the
  // entire survey (if the method is called without any arguments) or for a
  // given region.
  double Wtheta(int16_t region = -1);
  double WthetaError(int16_t region = -1);
  double WeightedCrossCorrelation(int16_t region = -1);
  double Weight(int16_t region = -1);
  uint32_t Counter(int16_t region = -1);
  double GalGal(int16_t region = -1);
  double GalRand(int16_t region = -1);
  double RandGal(int16_t region = -1);
  double RandRand(int16_t region = -1);

  // For these getters, we calculate the average value over all of the region
  // measurements.
  double MeanWtheta();
  double MeanWthetaError();
  double MeanWeightedCrossCorrelation();
  double MeanWeightedCrossCorrelationError();
  double MeanWeight();
  double MeanCounter();
  double MeanGalGal();
  double MeanGalRand();
  double MeanRandGal();
  double MeanRandRand();

  // Finally, some static methods which the AngularCorrelation method will use
  // to order its vectors of AngularBin objects.
  static bool ThetaOrder(AngularBin theta_a, AngularBin theta_b);
  static bool SinThetaOrder(AngularBin theta_a, AngularBin theta_b);
  static bool ReverseResolutionOrder(AngularBin theta_a,
				     AngularBin theta_b);

 private:
  double theta_min_, theta_max_, theta_;
  double costheta_min_, costheta_max_, sin2theta_min_, sin2theta_max_;
  double weight_, gal_gal_, gal_rand_, rand_gal_, rand_rand_;
  double pixel_wtheta_, pixel_weight_, wtheta_, wtheta_error_;
  uint32_t counter_;
  std::vector<double> weight_region_, gal_gal_region_;
  std::vector<double> gal_rand_region_, rand_gal_region_, rand_rand_region_;
  std::vector<double> pixel_wtheta_region_, pixel_weight_region_;
  std::vector<double> wtheta_region_, wtheta_error_region_;
  std::vector<uint32_t> counter_region_;
  uint16_t resolution_;
  int16_t n_region_;
  bool set_wtheta_error_, set_wtheta_;
};

} // end namespace Stomp

#endif
