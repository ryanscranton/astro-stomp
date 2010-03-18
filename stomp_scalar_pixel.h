// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains a variant on the Pixel class.  The goal here is to
// encode some manner of scalar field (galaxy density on the sky, CMB
// temperature, etc.) where we need both the value and an associated noise on
// the value.  Likewise, as will be seen in the ScalarMap class, these pixels
// will generally be used to sample this scalar field uniformly over some
// region.  This is in contrast with the Map object where the goal is to
// accurately describe a region's geometry using pixels of various sizes.

#ifndef STOMP_SCALAR_PIXEL_H
#define STOMP_SCALAR_PIXEL_H

#include "stomp_pixel.h"

namespace Stomp {

class ScalarPixel;

typedef std::vector<ScalarPixel> ScalarVector;
typedef ScalarVector::iterator ScalarIterator;
typedef std::pair<ScalarIterator, ScalarIterator> ScalarPair;
typedef std::vector<ScalarPixel *> ScalarPtrVector;
typedef ScalarPtrVector::iterator ScalarPtrIterator;

class ScalarPixel : public Pixel {
  // In order to do correlation function calculations, we need some
  // functionality beyond the normal Pixel object.  In particular, we want
  // to be able to encode fields, which may take one of three forms:
  //
  // * Pure scalar quantities (e.g. CMB temperature or radio flux).
  // * Point-field densities (e.g. the projected galaxy density over some area).
  // * Point-sampled averages (e.g. the mean galaxy magnitude over some area).
  //
  // In order to accomodate those three cases, we need an extra float and
  // an extra int (the Weight() value in Pixel will encode the fraction of
  // the pixel area that's contained in the survey area).  Pixels are taken
  // to be units of a Map where the geometry of the map is the union of all of
  // the Pixels and the Pixels are not assumed to be at the same resolution
  // level.  ScalarPixels, OTOH, form the basis for ScalarMaps where the
  // map is taken to be a regular sampling of some field over a given area.
  // The total area for the map can be calculated, but operations like
  // determining whether or not a given position is inside or outside the
  // map is not generically available.
 public:
  ScalarPixel();
  ScalarPixel(const uint32_t resolution, const uint32_t pixnum,
	      const double weight = 0.0, const double intensity = 0.0,
	      const uint32_t n_points = 0);
  ScalarPixel(AngularCoordinate& ang, const uint32_t resolution,
	      const double weight = 0.0, const double intensity = 0.0,
	      const uint32_t n_points = 0);
  ScalarPixel(const uint32_t x, const uint32_t y,
	      const uint32_t resolution, const double weight = 0.0,
	      const double intensity = 0.0, const uint32_t n_points = 0);
  virtual ~ScalarPixel();

  // To encode the three usage cases, we use an extra float (intensity) and
  // an int (n_points).  In the pure scalar field case, we ignore the
  // n_points variable and merely modify the intensity.  If we are treating
  // the pixel as a container for something like galaxy density, then we also
  // only care about the intensity, but when it comes to using this for doing
  // correlation calculations, we'll want to use a different form for the
  // over-density than for the scalar field case.  Finally, the point-sampled
  // field case requires us to keep track of the number of points we've added
  // to the pixel so that we can calculate the field average later on.
  void SetIntensity(const double intensity);
  void SetNPoints(const uint32_t n_point);
  double Intensity();
  uint32_t NPoints();

  // If we've gotta scalar field encoded on the pixel, the number of points
  // will be zero, so we just return the intensity value.  Otherwise, return
  // the average intensity for all of the points added to the pixel.
  double MeanIntensity();
  void AddToIntensity(const double intensity, const uint32_t n_point = 1);
  void ScaleIntensity(double scale_factor);

  // Once we're done adding points to the pixel, we may want to normalize to a
  // scalar field value so that we can access this field through the
  // Intensity() method.
  void NormalizeIntensity();

  // For the first and third usage cases, this is the form of the over-density
  // we'll want to use.
  void ConvertToOverDensity(double expected_intensity);

  // For the second use case (something like a galaxy density field), we
  // want to normalize our over-density by the average field (density) value.
  void ConvertToFractionalOverDensity(double expected_intensity);

  // And two complementary methods to take us from over-densities back to raw
  // intensities.
  void ConvertFromOverDensity(double expected_intensity);
  void ConvertFromFractionalOverDensity(double expected_intensity);

  // Finally, a method to tell us whether we're in over-density mode or not.
  bool IsOverDensity();

 private:
  double intensity_;
  uint32_t n_point_;
  bool is_overdensity_;
};

} // end namespace Stomp

#endif
