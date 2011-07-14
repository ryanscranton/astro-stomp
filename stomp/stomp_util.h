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

#ifndef STOMP_UTIL_H
#define STOMP_UTIL_H

#include <sys/time.h>
#include <string>

namespace Stomp {

class Cosmology;
class StompWatch;
class HistogramBin;
class Histogram;

typedef std::vector<HistogramBin> BinVector;
typedef BinVector::iterator BinIterator;

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
  // A simple class to handle time keeping inside of a program.  The return
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

class HistogramBin {
  // Simple class for storing information related to one of the bins in our
  // Histogram class.

public:
  HistogramBin();
  HistogramBin(double bin_min, double bin_max);
  ~HistogramBin();

  // If the bin is instantiated without bounds, we need to set them manually.
  // This automatically removes any items already added to the bin.
  void SetBounds(double bin_min, double bin_max);

  // Add an object to the bin.  Optionally, a weight can be attached to the
  // object.  The bin stores the total number of objects added as well as the
  // aggregate weight of those objects.  The returned boolean indicates whether
  // or not the bin value given was within the bin's bounds.
  bool AddToBin(double bin_value, double weight = 1.0);

  // Bin range values.
  double BinMinimum();
  double BinMaximum();

  // Our binning scheme may be linear or logrithmic, so we offer both options
  // for indicating the center of the bin.  We can also calculate the center of
  // the bin based on the weighted average of the bin values of the objects that
  // have been added to this bin or the straight average of the bin values.
  double BinCenter();
  double BinLogCenter();
  double BinWeightedCenter();
  double BinAveragedCenter();

  // Simple boolean to indicate whether an input bin value is within the bin's
  // bounds or not.
  bool WithinBin(double bin_value);

  // Some simple accessors for the aggregate state of what's been added to the
  // bin.
  double BinWeight();
  uint32_t BinItems();
  double BinMeanWeight();

private:
  double bin_min_, bin_max_, total_weight_;
  double total_bin_value_, weighted_bin_value_;
  uint32_t total_items_;
};

class Histogram {
  // Simple wrapper class for creating a histogram.  Items can be added to the
  // histogram with an associated weight such that each bin records the number
  // of times that something has been added to the bin as well as the aggregate
  // weight of those items added to the bin.

public:
  Histogram();
  Histogram(double hist_min, double hist_max, uint16_t n_bins,
	    bool log_binning = false);
  ~Histogram();

  // If the histogram is instantiated without bounds, we need to set them
  // manually.  This automatically removes any items already added to the
  // histogram.  Likewise, if the number of bins or type of binning is
  // changed.  If the current bounds are incompatible with logrithmic binning,
  // then the last method will return false.
  void SetBounds(double bin_min, double bin_max);
  void SetNBins(uint16_t n_bins);
  bool SetLogBinning(bool log_binning);
  bool AssignBins();

  // Attempt to add an item to the histogram.  In addition to having a bin
  // value, items can also have an associated weight.  The boolean indicates
  // whether the object was within the histogram's bounds and successfully
  // added to the histogram.
  bool AddToBin(double bin_value, double weight = 1.0);

  // Accessor methods for iterators over the histogram's bins.
  BinIterator Begin();
  BinIterator End();

  // Aggregate statistics for the histogram
  double BoundMin();
  double BoundMax();
  uint16_t NBins();
  bool LogBinning();

  // Aggregate statistics for the items in the histogram.
  uint32_t TotalItems();
  double TotalWeight();
  double MeanItemWeight();
  double MeanBinValue();
  double MeanWeightedBinValue();

private:
  double hist_min_, hist_max_;
  bool log_binning_;
  uint16_t n_bins_;
  BinVector bins_;
};


// A simple utility function for breaking up a string into its component items,
// generally when the string is a list of filenames or such (separated by the
// input delimiter).
void Tokenize(const std::string& str,
              std::vector<std::string>& tokens,
              const std::string& delimiters);

} // end namespace Stomp

#endif
