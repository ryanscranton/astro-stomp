#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include "stomp_core.h"
#include "stomp_util.h"

namespace Stomp {

// Define our default WMAP5 flat, LCDM cosmology.
double Cosmology::omega_m = 0.2736;
double Cosmology::h = 0.705;
const double Cosmology::AA_ = 1.718;
const double Cosmology::BB_ = 0.315;
double Cosmology::a_ = Cosmology::AA_*Cosmology::omega_m;
double Cosmology::b_ = Cosmology::BB_*sqrt(Cosmology::omega_m);

double Cosmology::OmegaM() {
  return omega_m;
}

double Cosmology::HubbleConstant() {
  return h*100.0;
}

double Cosmology::HubbleDistance() {
  return 3000.0/h;
}

double Cosmology::OmegaL() {
  return 1.0 - omega_m;
}

void Cosmology::SetOmegaM(double new_omega_m) {
  omega_m = new_omega_m;
  a_ = AA_*omega_m;
  b_ = BB_*sqrt(omega_m);
}

void Cosmology::SetHubbleConstant(double hubble) {
    h = hubble/100.0;
}

void Cosmology::SetOmegaL(double omega_lambda) {
  omega_m = 1.0 - omega_lambda;
  a_ = AA_*omega_m;
  b_ = BB_*sqrt(omega_m);
}

double Cosmology::ComovingDistance(double z) {
  // In Mpc/h.
  return HubbleDistance()*z/sqrt(1.0 + a_*z + b_*z*z);
}

double Cosmology::AngularDiameterDistance(double z) {
  // In Mpc/h.
  return ComovingDistance(z)/(1.0+z);
}

double Cosmology::LuminosityDistance(double z) {
  // In Mpc/h.
  return ComovingDistance(z)*(1.0+z);
}

double Cosmology::ProjectedDistance(double z, double theta) {
  // where theta is assumed to be in degrees and the return value in Mpc/h.
  return theta*DegToRad*AngularDiameterDistance(z);
}

double Cosmology::ProjectedAngle(double z, double radius) {
  // where radius is in Mpc/h and the return angle is in degrees.
  return RadToDeg*radius/AngularDiameterDistance(z);
}

StompWatch::StompWatch() {
  gettimeofday(&start, NULL);
  stop = start;
}

void StompWatch::StartTimer() {
  gettimeofday(&start, NULL);
}

void StompWatch::StopTimer() {
  gettimeofday(&stop, NULL);
}

double StompWatch::ElapsedTime() {
  timeval interval;
  timersub(&stop, &start, &interval);
  // Need to convert tv_usec from microseconds to seconds.
  return interval.tv_sec + interval.tv_usec/1000000.0;
}

HistogramBin::HistogramBin() {
  total_items_ = 0;
  total_weight_ = 0.0;
  total_bin_value_ = 0.0;
  bin_min_ = -1.0;
  bin_max_ = -1.0;
}

HistogramBin::HistogramBin(double bin_min, double bin_max) {
  total_items_ = 0;
  total_weight_ = 0.0;
  total_bin_value_ = 0.0;
  bin_min_ = bin_min;
  bin_max_ = bin_max;
}

HistogramBin::~HistogramBin() {
  total_items_ = 0;
  total_weight_ = 0.0;
  weighted_bin_value_ = 0.0;
  bin_min_ = -1.0;
  bin_max_ = -1.0;
}

void HistogramBin::SetBounds(double bin_min, double bin_max) {
  total_items_ = 0;
  total_weight_ = 0.0;
  total_bin_value_ = 0.0;
  weighted_bin_value_ = 0.0;
  bin_min_ = bin_min;
  bin_max_ = bin_max;
}

bool HistogramBin::AddToBin(double bin_value, double weight) {
  bool added_to_bin = false;

  if (WithinBin(bin_value)) {
    total_bin_value_ += bin_value;
    weighted_bin_value_ += bin_value*weight;
    total_weight_ += weight;
    total_items_++;
    added_to_bin = true;
  }

  return added_to_bin;
}

double HistogramBin::BinMinimum() {
  return bin_min_;
}

double HistogramBin::BinMaximum() {
  return bin_max_;
}

double HistogramBin::BinCenter() {
  return 0.5*(bin_min_ + bin_max_);
}

double HistogramBin::BinLogCenter() {
  return pow(10.0, 0.5*(log10(bin_min_) + log10(bin_max_)));
}

double HistogramBin::BinWeightedCenter() {
  return weighted_bin_value_/total_weight_;
}

double HistogramBin::BinAveragedCenter() {
  return total_bin_value_/static_cast<double>(total_items_);
}

bool HistogramBin::WithinBin(double bin_value) {
  return (DoubleLE(bin_value, bin_max_) &&
	  DoubleGE(bin_value, bin_min_) ? true : false);
}

uint32_t HistogramBin::BinItems() {
  return total_items_;
}

double HistogramBin::BinWeight() {
  return total_weight_;
}

Histogram::Histogram() {
  hist_min_ = -1.0;
  hist_max_ = -1.0;
  n_bins_ = 0;
  log_binning_ = false;
  bins_.clear();
}

Histogram::Histogram(double hist_min, double hist_max,
		     uint16_t n_bins, bool log_binning) {
  hist_min_ = hist_min;
  hist_max_ = hist_max;
  n_bins_ = n_bins;
  log_binning_ = log_binning;

  AssignBins();
}

Histogram::~Histogram() {
  hist_min_ = -1.0;
  hist_max_ = -1.0;
  n_bins_ = 0;
  log_binning_ = false;
  bins_.clear();
}

bool Histogram::AssignBins() {
  if (!bins_.empty()) bins_.clear();
  bins_.reserve(n_bins_);

  bool assigned_bins = true;

  if (!log_binning_) {
    double dbin = (hist_max_ - hist_min_)/n_bins_;
    for (uint16_t i=0;i<n_bins_;i++) {
      bins_.push_back(HistogramBin(hist_min_ + dbin*i,
				   hist_min_ + dbin*(i+1)));
    }
  } else {
    if (DoubleLT(hist_min_, 0.0)) {
      hist_min_ = 1.0e-10;
      assigned_bins = false;
    }
    if (DoubleLT(hist_max_, 0.0)) {
      hist_min_ = 1.0e-5;
      assigned_bins = false;
    }
    double dlogbin = (log10(hist_max_) - log10(hist_min_))/n_bins_;
    for (uint16_t i=0;i<n_bins_;i++) {
      double log_histmin = log10(hist_min_) + i*dlogbin;
      double log_histmax = log10(hist_min_) + (i+1)*dlogbin;
      bins_.push_back(HistogramBin(pow(10.0, log_histmin),
				   pow(10.0, log_histmax)));
    }
  }

  return assigned_bins;
}

bool Histogram::AddToBin(double bin_value, double weight) {
  bool found_bin = false;

  // Do a simple binary search on the bins to find the appropriate one.
  BinIterator top = bins_.end();
  BinIterator bottom = bins_.begin();
  BinIterator iter;
  --top;

  if ((bin_value < bottom->BinMinimum()) ||
      (bin_value > top->BinMaximum())) {
    iter = bins_.end();
  } else {
    ++top;
    --bottom;
    while (top-bottom > 1) {
      iter = bottom + (top - bottom)/2;
      if (bin_value < iter->BinMinimum()) {
	top = iter;
      } else {
	bottom = iter;
      }
    }
    iter = bottom;
  }

  // Provided that we've found the bin where the input belongs, add the
  // galaxies and area to that bin.
  if (iter!=bins_.end() && iter->WithinBin(bin_value)) {
    found_bin = iter->AddToBin(bin_value);
  }

  return found_bin;
}

BinIterator Histogram::Begin() {
  return bins_.begin();
}

BinIterator Histogram::End() {
  return bins_.end();
}

double Histogram::BoundMin() {
  return hist_min_;
}

double Histogram::BoundMax() {
  return hist_max_;
}

uint16_t Histogram::NBins() {
  return n_bins_;
}

bool Histogram::LogBinning() {
  return log_binning_;
}

uint32_t Histogram::TotalItems() {
  uint32_t total_items = 0;

  for (BinIterator iter=bins_.begin();iter!=bins_.end();++iter) {
    total_items += iter->BinItems();
  }

  return total_items;
}

double Histogram::TotalWeight() {
  double total_weight = 0.0;

  for (BinIterator iter=bins_.begin();iter!=bins_.end();++iter) {
    total_weight += iter->BinWeight();
  }

  return total_weight;
}

double Histogram::MeanItemWeight() {
  uint32_t total_items = 0;
  double total_weight = 0.0;

  for (BinIterator iter=bins_.begin();iter!=bins_.end();++iter) {
    total_items += iter->BinItems();
    total_weight += iter->BinWeight();
  }

  return total_weight/total_items;
}

double Histogram::MeanBinValue() {
  uint32_t total_items = 0;
  double total_bin_value = 0.0;

  for (BinIterator iter=bins_.begin();iter!=bins_.end();++iter) {
    total_items += iter->BinItems();
    total_bin_value += iter->BinAveragedCenter()*iter->BinItems();
  }

  return total_bin_value/total_items;
}

double Histogram::MeanWeightedBinValue() {
  double total_weight = 0.0;
  double weighted_bin_value = 0.0;

  for (BinIterator iter=bins_.begin();iter!=bins_.end();++iter) {
    total_weight += iter->BinWeight();
    weighted_bin_value += iter->BinWeightedCenter()*iter->BinWeight();
  }

  return weighted_bin_value/total_weight;
}

void Tokenize(const std::string& str, std::vector<std::string>& tokens,
	 const std::string& delimiters) {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
  }
}

} // end namespace Stomp
