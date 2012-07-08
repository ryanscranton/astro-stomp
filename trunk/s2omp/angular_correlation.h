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

#ifndef ANGULAR_CORRELATION_H_
#define ANGULAR_CORRELATION_H_

#include <vector>
#include "core.h"
#include "point.h"
#include "angular_bin-inl.h"

namespace s2omp {

class pixel_union; // class declaration in stomp_map.h
class scalar_union; // class declaration in stomp_scalar_map.h
class tree_union; // class declaration in stomp_tree_map.h
class angular_correlation;

class angular_correlation {
	// Class object for calculating auto-correlations and cross-correlations
	// given a set of objects and a Map.  Broadly speaking, this is a
	// container class for a set of AngularBin objects which collectively
	// span some range of angular scales.  Accordingly, the methods are generally
	// intended to package the machinery of the auto-correlation and
	// cross-correlation calculations into simple, one-line calls.

public:
	enum Estimator {
		PAIR, PIXEL, HYBRID
	};

	// The first constructor takes an angular minimum and maximum (in degrees)
	// and constructs a logrithmic binning scheme using the specified number
	// of bins per decade (which can be a non-integer value, obviously).  The
	// bins are such that the minimum angular scale of the first bin will be
	// theta_min and the maximum angular scale of the last bin with be
	// theta_max.  The last boolean argument controls whether or not an
	// pixel resolution will be assigned to the bins.  If it is false, then
	// the resolution values will all be -1.
	angular_correlation* log_binning(double theta_min, double theta_max,
			double bins_per_decade, Estimator e);

	// The alternate constructor is used for a linear binning scheme.  The
	// relationship between theta_min and theta_max remains the same and the
	// spacing of the bins is determined based on the requested number of bins.
	angular_correlation* linear_binning(double theta_min, double theta_max,
			uint16_t n_bins, Estimator e);
	virtual ~angular_correlation() {
		thetabin_.clear();
	}

	// For small angular scales, it's usually faster and more memory
	// efficient to use a pair-based estimator.  To set this scale, we choose
	// a maximum resolution scale we're willing to use our pixel-based estimator
	// on and modify all smaller angular bins to use the pair-based estimator.
	// The boolean indicates to the object whether this break between the two
	// estimators is being set by hand (default) or should be over-ridden if the
	// methods for calculating the correlation functions are called.
	void set_max_level(int level, bool manual_break);

	// Additionally, if we are using regions to calculate correlation functions,
	// we need to set the minimum resolution to match the resolution used to
	// divide the total survey area.
	void set_min_level(int level);

	// If we haven't set the break manually, this method attempts to find a
	// reasonable place for it, based on the number of objects involved in the
	// correlation function calculation and the area involved.
	void auto_max_level(uint32_t n_obj, double area_deg2);

	// If we're going to use regions to find jack-knife errors, then we need
	// to initialize the AngularBins to handle this state of affairs or possibly
	// clear out previous calculations.
	void init_regions(int16_t n_regions);
	void clear_regions();
	int16_t n_region();

	// Some wrapper methods for finding the auto-correlation and cross-correlations
	void find_autocorrelation(const pixel_union& pixels,
			const point_vector& points, uint8_t n_iterations, uint16_t n_regions);
	void find_autocorrelation(const geometric_bound_interface& bound,
			const point_vector& points, uint8_t n_iterations, uint16_t n_regions);
	void find_crosscorrelation(const pixel_union& pixels,
			const point_vector& points_a, const point_vector& points_b,
			uint8_t n_iterations, uint16_t n_regions);
	void find_crosscorrelation(const geometric_bound_interface& bound,
			const point_vector& points_a, const point_vector& points_b,
			uint8_t n_iterations, uint16_t n_regions);

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
	void
	FindPixelCrossCorrelation(ScalarMap& stomp_map_a, ScalarMap& stomp_map_b);
	void FindPairAutoCorrelation(Map& stomp_map, WAngularVector& galaxy,
			uint8_t random_iterations = 1);
	void FindPairCrossCorrelation(Map& stomp_map, WAngularVector& galaxy_a,
			WAngularVector& galaxy_b, uint8_t random_iterations = 1);

	// Once we're done calculating our correlation function, we can write it out
	// to an ASCII file.  The output format will be
	//
	//   THETA   W(THETA)   dW(THETA)
	//
	// where THETA is the angular scale in degrees and dW(THETA) is the jack-knife
	// error based on regionating the data.  If the angular correlation has been
	// calculated without regions, then this column will be omitted.
	bool write(const std::string& output_file_name);

	// In some cases, we want to default to using either the pair-based or
	// pixel-based estimator for all of our bins, regardless of angular scale.
	// These methods allow us to over-ride the default behavior of the
	// correlation code.
	void use_only_pixels();
	void use_only_pairs();

	// Now, some accessor methods for finding the angular range of the bins
	// with a given resolution attached to them (the default value returns the
	// results for all angular bins; for pair-based bins, resolution = -1).
	double theta_min() {
		return theta_min_;
	}
	double theta_min(int level);
	double theta_min(Estimator e);
	double theta_max() {
		return theta_max_;
	}
	double theta_max(int level);
	double theta_max(Estimator e);

	double sin2theta_min() {
		return sin2theta_min_;
	}
	double sin2theta_min(int level);
	double sin2theta_min(Estimator e);
	double sin2theta_max() {
		return sin2theta_max_;
	}
	double sin2theta_max(int level);
	double sin2theta_max(Estimator e);

	theta_iterator begin() {
		return thetabin_.begin();
	}
	theta_iterator begin(int level);
	theta_iterator begin(Estimator e);
	theta_iterator end() {
		return thetabin_.end();
	}
	theta_iterator end(int level);
	theta_iterator end(Estimator e);

	theta_iterator find(theta_iterator begin, theta_iterator end,
			double sin2theta);

	theta_iterator bin_iterator(uint8_t bin_idx);

	inline uint32_t n_bins() {
		return theta_bin_.size();
	}
	inline uint32_t min_level() {
		return min_level_;
	}
	inline uint32_t max_level() {
		return max_level_;
	}

	// The AngularBin accessor methods work fine for calculations that can be
	// done using the data from a single AngularBin, but for the covariance
	// matrix for our measurement we need a method that can access
	// multiple bins at once.  This method returns the (theta_a, theta_b) of the
	// covariance matrix.
	double covariance(uint8_t bin_idx_a, uint8_t bin_idx_b);

	// Alternatively, we can just write the full covariance matrix to a file.
	// The output format will be
	//
	//   THETA_A  THETA_B  Cov(THETA_A, THETA_B)
	//
	bool write_covariance(const std::string& output_file_name);

private:
	// Find the resolution we would use to calculate correlation functions for
	// each of the bins.  If this method is not called, then the resolution
	// for each bin is set to -1, which would indicate that any correlation
	// calculation with that bin should be done using a pair-based estimator.
	void assign_bin_levels(int max_level);

	theta_vector thetabin_;
	Estimator estimator_;
	theta_iterator theta_pixel_begin_, theta_pixel_end_;
	theta_iterator theta_pair_begin_, theta_pair_end_;
	double theta_min_, theta_max_, sin2theta_min_, sin2theta_max_;
	int min_level_, max_level_, regionation_level_;
	int16_t n_region_;
	bool manual_resolution_break_;
};

} // end namespace s2omp


#endif /* ANGULAR_CORRELATION_H_ */
