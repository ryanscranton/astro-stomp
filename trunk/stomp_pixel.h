// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the central class for the library: the hierarchical
// pixelization that makes all of the rest of the spatial classes work.  This
// class defines all of the core Pixel operations: converting angular position
// on the sphere into pixel index, increasing and decreasing pixel resolution,
// finding neighboring pixels and so on.

#ifndef STOMP_PIXEL_H
#define STOMP_PIXEL_H

#include <stdint.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "MersenneTwister.h"

namespace Stomp {

class AngularCoordinate;
class AngularBin;
class Pixel;

typedef std::vector<Pixel> PixelVector;
typedef PixelVector::iterator PixelIterator;
typedef std::pair<PixelIterator, PixelIterator> PixelPair;
typedef std::vector<Pixel *> PixelPtrVector;
typedef PixelPtrVector::iterator PixelPtrIterator;

class Pixel {
  // The core class for this library.  An instance of this class represents
  // a single pixel covering a particular region of the sky, with a particular
  // weight represented by a float.  Pixels can be instantiated with an
  // AngularCoordinate and resolution level or pixel indices or just
  // instantiated with no particular location.

 public:
  Pixel();
  Pixel(const uint32_t resolution, const uint32_t pixnum,
	const double weight = 0.0);
  Pixel(AngularCoordinate& ang, const uint32_t resolution,
	const double weight = 0.0);
  Pixel(const uint32_t x, const uint32_t y,
	const uint32_t resolution, const double weight = 0.0);
  virtual ~Pixel();

  // For the purposes of simple ordering (as would be done in various STL
  // containers and the like), we always use the equivalent of
  // LocalOrdering defined below.
  bool operator<(Pixel& pix);

  // Likewise, when testing equality, we ignore the pixel weight and just
  // concentrate on geometry.
  bool operator==(Pixel& pix);
  bool operator!=(Pixel& pix);

  // For the iterators, we use the same LocalOrdering principle, with the
  // additional caveat that a pixel doesn't go beyond the maximum allowed
  // index+1.
  Pixel operator++();

  // Some basic setters and getters.
  void SetPixnumFromAng(AngularCoordinate& ang);
  void SetResolution(uint32_t resolution);
  void SetLevel(uint8_t level);
  void SetPixnumFromXY(uint32_t x, uint32_t y);
  void SetWeight(double weight);

  uint8_t Level();
  uint32_t Resolution();
  uint32_t PixelX();
  uint32_t PixelY();
  double Weight();

  // And operations on the Weight value associated with the Pixel.
  void ReverseWeight();
  void InvertWeight();

  // These methods all relate the hierarchical nature of the pixelization
  // method.  SetToSuperPix degrades a high resolution pixel into the lower
  // resolution pixel that contains it.  SetToLevel does the same, but takes
  // Level() as an argument.  SubPix returns either a vector of
  // higher resolution pixels contained by this pixel or the X-Y pixel index
  // bounds that will let one iterate through the sub-pixels without
  // instantiating them.
  bool SetToSuperPix(uint32_t lo_resolution);
  bool SetToLevel(uint8_t lo_level);
  void SubPix(uint32_t hi_resolution, PixelVector& pix);
  void SubPix(uint32_t hi_resolution, uint32_t& x_min, uint32_t& x_max,
	      uint32_t& y_min, uint32_t& y_max);

  // CohortPix returns the pixels at the same resolution that would combine
  // with the current pixel to form the pixel at the next coarser resolution
  // level.
  void CohortPix(Pixel& pix_a, Pixel& pix_b, Pixel& pix_c);

  // FirstCohort returns a boolean indicating true if the current pixel would
  // be the first pixel in a sorted list of the cohort pixels and false
  // otherwise.
  bool FirstCohort();

  // Since the pixels are equal-area, we only need to know how many times
  // we've sub-divided to get from the HpixResolution to our current resolution.
  double Area();

  // This returns the index of the pixel that contains the current pixel at
  // a coarser resolution.
  uint32_t SuperPix(uint32_t lo_resolution);

  // The fundamental unit of the spherical pixelization.  There are 7488
  // superpixels covering the sky, which makes isolating any searches to just
  // the pixels within a single superpixel a very quick localization strategy.
  uint32_t Superpixnum();

  // Single index ordering within a superpixel.  The cast finds the x-y
  // position of the superpixel and then we scale that up to the
  // pseudo resolution within the superpixel (where HPixResolution() is
  // effectively 0 and we scale up from there.  This is a cheat that lets
  // us single index pixels up to a higher resolution without going to a
  // 64 bit index.
  uint32_t HPixnum();

  // Single index ordering for the whole sphere.  Unforunately, the limits of
  // a 32 bit integer only take us up to about 17 arcsecond resolution.
  uint32_t Pixnum();

  // Given either the X-Y-resolution, Pixel or AngularCoordinate, return
  // true or false based on whether the implied location is within the current
  // Pixel.
  bool Contains(uint32_t pixel_resolution, uint32_t pixel_x,
		uint32_t pixel_y);
  bool Contains(Pixel& pix);
  bool Contains(AngularCoordinate& ang);

  // Given a set of lon-lat bounds, return true/false if the pixel is within
  // those bounds.
  bool WithinBounds(double lon_min, double lon_max,
		    double lat_min, double lat_max,
		    AngularCoordinate::Sphere sphere);

  // And slightly more permissive version that checks to see if any part of
  // the pixel is within the input bounds.
  bool IntersectsBounds(double lon_min, double lon_max,
			double lat_min, double lat_max,
			AngularCoordinate::Sphere sphere);

  // Given an angle in degrees (or upper and lower angular bounds in degrees),
  // return a list of pixels at the same resolution within those bounds.
  virtual void WithinRadius(double theta_max, PixelVector& pix,
			    bool check_full_pixel=false);
  virtual void WithinAnnulus(double theta_min, double theta_max,
			     PixelVector& pix, bool check_full_pixel=false);
  virtual void WithinAnnulus(AngularBin& theta, PixelVector& pix,
			     bool check_full_pixel=false);

  // While the methods above are useful, sometimes we want a bit less precision.
  // The return vector of pixels include any possible pixels that might be
  // within the input radius.  We also include options for specifiying an
  // AngularCoordinate other than the one at the center of the pixel.
  void BoundingRadius(double theta_max, PixelVector& pix);
  void BoundingRadius(AngularCoordinate& ang, double theta_max,
		      PixelVector& pix);

  // Similar to the previous methods, but the return values here are the
  // X-Y indices rather than the pixels themselves.  The second instance allows
  // the X index bounds to vary with Y index value to take into account the
  // effects of curvature on the sphere, while the first instance just uses
  // the outer limits.  The third and fourth options allow for centering the
  // X-Y bounds at some angular position other than the pixel center.
  void XYBounds(double theta, uint32_t& x_min, uint32_t& x_max,
		uint32_t& y_min, uint32_t& y_max,
		bool add_buffer = false);
  void XYBounds(double theta, std::vector<uint32_t>& x_min,
		std::vector<uint32_t>& x_max,
		uint32_t& y_min, uint32_t& y_max,
		bool add_buffer = false);
  void XYBounds(AngularCoordinate& ang, double theta,
		uint32_t& x_min, uint32_t& x_max,
		uint32_t& y_min, uint32_t& y_max,
		bool add_buffer = false);
  void XYBounds(AngularCoordinate& ang, double theta,
		std::vector<uint32_t>& x_min,
		std::vector<uint32_t>& x_max,
		uint32_t& y_min, uint32_t& y_max,
		bool add_buffer = false);
  uint8_t EtaStep(double theta);

  // Additionally, it can be useful to know the projected angular distance
  // between an angular location and the nearest edge of the pixel.  This
  // routine returns the maximum projected distance to the closest edge (i.e.,
  // if the projected distance to the nearest edge in eta is small, but large
  // in lambda, the returned value is the lambda value, indicating that the
  // effective angular distance to the area covered by pixel is large).  For
  // cases where you want to know whether any part of the pixel is within some
  // fixed angular distance, this is the relavent quantity.
  //
  // The second version does the same, but for the far edge.  The combination
  // of the two should be sufficient for determining if a pixel is within some
  // angular radius as well as determining if a part of the pixel might be
  // within an angular annulus.
  //
  // In both cases, the return value is (sin(theta))^2 rather than just the
  // angle theta.  For small angles (where this method is most likely used),
  // this is a more useful quantity from a computing speed standpoint.
  double NearEdgeDistance(AngularCoordinate& ang);
  double FarEdgeDistance(AngularCoordinate& ang);

  // Likewise, we also want to be able to find the near and far corner
  // distances if necessary.  This is more expensive because there's a trig
  // function involved in generating the distances.
  double NearCornerDistance(AngularCoordinate& ang);
  double FarCornerDistance(AngularCoordinate& ang);

  // Alternatively, we can get both edge distances in a single call.  The
  // returned boolean tells us whether the distances are to pixel edges (true)
  // or pixel corners (false).
  bool EdgeDistances(AngularCoordinate& ang, double& near_edge_distance,
		     double& far_edge_distance);


  // Finishing up the angle-checking methods, we have two more methods that
  // return true or false based on whether the current pixel is within a given
  // angular range of a point on the sphere (specified by either a raw angular
  // coordinate or the center point of another Pixel).  For the slippery
  // case where an annulus might only partially contain the pixel, we have
  // the last two methods.  A return value of 1 indicates that the pixel is
  // fully within the annulus (by calling IsWithinAnnulus with full pixel
  // checking), a return value of 0 indicates that it is fully outside the
  // annulus and a return value of -1 indicates a partial containment.
  bool IsWithinRadius(AngularCoordinate& ang, double theta_max,
		      bool check_full_pixel=false);
  bool IsWithinRadius(Pixel& pix, double theta_max,
		      bool check_full_pixel=false);
  bool IsWithinAnnulus(AngularCoordinate& ang, double theta_min,
		       double theta_max, bool check_full_pixel=false);
  bool IsWithinAnnulus(Pixel& pix, double theta_min, double theta_max,
		       bool check_full_pixel=false);
  bool IsWithinAnnulus(AngularCoordinate& ang, AngularBin& theta,
		       bool check_full_pixel=false);
  bool IsWithinAnnulus(Pixel& pix, AngularBin& theta,
		       bool check_full_pixel=false);
  int8_t IntersectsAnnulus(AngularCoordinate& ang,
				   double theta_min, double theta_max);
  int8_t IntersectsAnnulus(Pixel& pix,
			   double theta_min, double theta_max);
  virtual int8_t IntersectsAnnulus(AngularCoordinate& ang, AngularBin& theta);
  int8_t IntersectsAnnulus(Pixel& pix, AngularBin& theta);

  // A hold-over from the SDSS coordinate system, this converts the current
  // pixel index into an SDSS stripe number.  Although this is generally not
  // useful information in an of itself, stripe number is used as a proxy
  // for constructing roughly square subsections of Maps.
  uint32_t Stripe(uint32_t resolution = HPixResolution);

  // Some methods for extracting the angular position of the pixel center...
  double RA();
  double DEC();
  double GalLon();
  double GalLat();
  void Ang(AngularCoordinate& ang);
  double Lambda();
  double Eta();

  // ... likewise for the Cartesian coordinates on the unit sphere.
  virtual double UnitSphereX();
  virtual double UnitSphereY();
  virtual double UnitSphereZ();

  // Since the pixels are rectangular in survey coordinates, we have meaningful
  // notions of the bounds in lambda-eta space.
  double LambdaMin();
  double LambdaMax();
  double EtaMin();
  double EtaMax();

  // In some cases, it can be handy to have a variant on the previous method
  // where we return a value of EtaMax that might technically violate the
  // proper bounds on Eta (-180 < Eta < 180), but is continous with the
  // value of EtaMin (i.e. EtaMin < EtaMax, regardless of the discontinuity).
  double EtaMaxContinuous();

  // Alternatively, we may just want to know if the pixel crosses the Eta
  // discontinuity.
  bool SurveyContinuous();

  // Corresponding bounds for Equatorial and Galactic coordinates.  These are
  // more expensive since we need to check the bounds at each corner.
  double DECMin();
  double DECMax();
  double RAMin();
  double RAMax();
  double RAMaxContinuous();
  bool EquatorialContinuous();

  double GalLatMin();
  double GalLatMax();
  double GalLonMin();
  double GalLonMax();
  double GalLonMaxContinuous();
  bool GalacticContinuous();

  // And one last function to check continuity for a given angular coordinate
  // system.
  bool ContinuousBounds(AngularCoordinate::Sphere sphere);

  // And it can be useful to be able to quickly extract the x-y-z positions of
  // the pixel corners.
  virtual double UnitSphereX_UL();
  virtual double UnitSphereY_UL();
  virtual double UnitSphereZ_UL();

  virtual double UnitSphereX_UR();
  virtual double UnitSphereY_UR();
  virtual double UnitSphereZ_UR();

  virtual double UnitSphereX_LL();
  virtual double UnitSphereY_LL();
  virtual double UnitSphereZ_LL();

  virtual double UnitSphereX_LR();
  virtual double UnitSphereY_LR();
  virtual double UnitSphereZ_LR();

  void Iterate(bool wrap_pixel = true);

  // Like PixelX, but this returns the first x value for the current pixel's
  // superpixel (useful for knowing the bounds for SuperPixelBasedOrder'd lists.
  uint32_t PixelX0();

  // Same as PixelX0, but for the y index.
  uint32_t PixelY0();

  // This would be the x value just beyond the limit for the current pixel's
  // superpixel.  Hence, all the pixels with PixelX0 <= x < PixelX1 are in
  // the same column of superpixels.  For pixels in
  // superpixel = MaxSuperpixum, this is Nx0*_resolution, so it can be used
  // as an iteration bound for all superpixels.
  uint32_t PixelX1();

  // Same as PixelX1, but for the y index.
  uint32_t PixelY1();

  // Given a requested number of points, return a vector of Poisson random
  // angular positions within the current Pixel's area.
  void GenerateRandomPoints(AngularVector& ang, uint32_t n_point = 1);

  // This next block of code is there to provide backwards compatibility
  // to a straight C interface.  True, we're using references which aren't
  // in C, but these methods still allow users to access the basic interfaces
  // without instantiating a Pixel, which can be handy in some cases.
  static uint8_t Resolution2Level(const uint32_t resolution);
  static uint32_t Level2Resolution(const uint8_t level);
  static void Ang2Pix(const uint32_t resolution, AngularCoordinate& ang,
		      uint32_t& pixnum);
  static void Pix2Ang(uint32_t resolution, uint32_t pixnum,
                      AngularCoordinate& ang);
  static void Pix2HPix(uint32_t input_resolution, uint32_t input_pixnum,
		       uint32_t& output_hpixnum,
		       uint32_t& output_superpixnum);
  static void HPix2Pix(uint32_t input_resolution, uint32_t input_hpixnum,
		       uint32_t input_superpixnum,
		       uint32_t& output_pixnum);
  static void SuperPix(uint32_t hi_resolution, uint32_t hi_pixnum,
                       uint32_t lo_resolution, uint32_t& lo_pixnum);
  static void SubPix(uint32_t lo_resolution, uint32_t hi_pixnum,
		     uint32_t hi_resolution, uint32_t& x_min,
		     uint32_t& x_max, uint32_t& y_min,
		     uint32_t& y_max);
  static void NextSubPix(uint32_t input_resolution, uint32_t input_pixnum,
			 uint32_t& sub_pixnum1,
			 uint32_t& sub_pixnum2,
			 uint32_t& sub_pixnum3,
			 uint32_t& sub_pixnum4);
  static void AreaIndex(uint32_t resolution, double lammin, double lammax,
			double etamin, double etamax, uint32_t& x_min,
			uint32_t& x_max, uint32_t& y_min,
			uint32_t& y_max);
  static void PixelBound(uint32_t resolution, uint32_t pixnum, double& lammin,
			 double& lammax, double& etamin, double& etamax);
  static void CohortPix(uint32_t resolution, uint32_t hpixnum,
			uint32_t& pixnum1, uint32_t& pixnum2,
			uint32_t& pixnum3);
  static double PixelArea(uint32_t resolution);
  static uint8_t Pix2EtaStep(uint32_t resolution, uint32_t pixnum,
			     double theta);
  static void Ang2HPix(uint32_t resolution, AngularCoordinate& ang,
		       uint32_t& hpixnum, uint32_t& superpixnum);
  static void HPix2Ang(uint32_t resolution, uint32_t hpixnum,
                       uint32_t superpixnum, AngularCoordinate& ang);
  static void XY2HPix(uint32_t resolution, uint32_t x, uint32_t y,
                      uint32_t& hpixnum, uint32_t& superpixnum);
  static void HPix2XY(uint32_t resolution, uint32_t hpixnum,
		      uint32_t superpixnum, uint32_t& x,
		      uint32_t& y);
  static void SuperHPix(uint32_t hi_resolution, uint32_t hi_hpixnum,
                        uint32_t lo_resolution, uint32_t& lo_hpixnum);
  static void NextSubHPix(uint32_t resolution, uint32_t hpixnum,
			  uint32_t& hpixnum1, uint32_t& hpixnum2,
			  uint32_t& hpixnum3, uint32_t& hpixnum4);
  static void SubHPix(uint32_t lo_resolution, uint32_t hi_hpixnum,
		      uint32_t hi_superpixnum, uint32_t hi_resolution,
		      uint32_t& x_min, uint32_t& x_max,
		      uint32_t& y_min, uint32_t& y_max);
  static void HPixelBound(uint32_t resolution, uint32_t hpixnum,
			  uint32_t superpixnum, double& lammin,
			  double& lammax, double& etamin, double& etamax);
  static void CohortHPix(uint32_t resolution, uint32_t hpixnum,
			 uint32_t& hpixnum1, uint32_t& hpixnum2,
			 uint32_t& hpixnum3);
  static double HPixelArea(uint32_t resolution);
  static uint8_t HPix2EtaStep(uint32_t resolution, uint32_t hpixnum,
			      uint32_t superpixnum, double theta);
  static void XY2Pix(uint32_t resolution, uint32_t x, uint32_t y,
		     uint32_t& pixnum);
  static void Pix2XY(uint32_t resolution, uint32_t pixnum,
		     uint32_t& x, uint32_t& y);

  // Now we've got the various methods to establish ordering on the pixels.
  // LocalOrder is the the simplest, just arranging all of the pixels in
  // vanilla row-column order.  That's useful for some operations where you
  // want to be able to access nearby pixels simply.  However, if you're
  // doing a search on a large region of the sky, it often helps to be able to
  // limit the search more drastically at the outset.  For that, we have
  // SuperPixelBasedOrder where pixels are grouped by their lowest resolution
  // superpixel and then locally sorted within that bound.  This is the
  // default sorting method for the Map class to make searching on those
  // maps more efficient.  Finally, we have some methods for checking
  // whether or not we're looking at equivalent pixels, one where the weights
  // associated with the pixels matter and one that's purely geometric.
  static bool LocalOrder(const Pixel pix_a, const Pixel pix_b);
  static bool LocalOrderByReference(const Pixel pix_a, const Pixel pix_b);
  static bool SuperPixelBasedOrder(const Pixel pix_a, const Pixel pix_b);
  static bool SuperPixelOrder(const Pixel pix_a, const Pixel pix_b);
  static bool WeightedOrder(Pixel pix_a, Pixel pix_b);
  static bool WeightMatch(Pixel& pix_a, Pixel& pix_b);
  static bool WeightedPixelMatch(Pixel& pix_a, Pixel& pix_b);
  static bool PixelMatch(Pixel& pix_a, Pixel& pix_b);

  // Finally, these methods handle maps consisting of vectors of Pixels.
  // One could make the argument that this should be in the Map class,
  // but one of the primary purposes of these methods is to take a list of
  // pixels where there may be duplication or cases where smaller pixels are
  // within larger pixels and generate a set of pixels that uniquely covers
  // a give region of the sky.  That extra applicability makes it appropriate
  // to put here.  The main method to call is ResolvePixel, which will call
  // ResolveSuperPixel individually for each fo the superpixels covered by
  // the vector of Pixels.  The resulting vector will be sorted by
  // SuperPixelBasedOrder.
  static void ResolveSuperPixel(PixelVector& pix, bool ignore_weight = false);
  static void ResolvePixel(PixelVector& pix, bool ignore_weight = false);
  static void FindUniquePixels(PixelVector& input_pix, PixelVector& unique_pix);

 private:
  double weight_;
  uint32_t x_, y_;
  uint8_t level_;
};

} // end namespace Stomp

#endif
