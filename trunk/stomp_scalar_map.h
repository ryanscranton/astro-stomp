// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the ScalarMap class.  Unlike Maps, the primary
// goal here is to encode a scalar field over some area of the sky.  As such,
// we sacrifice some degree of precision in describing the exact area of the
// field and we use a uniform sampling of the field across the area in
// question.  This makes the class ideal for calculating angular correlation
// functions on the encoded field.

#ifndef STOMP_SCALAR_MAP_H
#define STOMP_SCALAR_MAP_H

#include <stdint.h>
#include <vector>
#include "stomp_core.h"
#include "stomp_angular_bin.h"
#include "stomp_scalar_pixel.h"
#include "stomp_base_map.h"

namespace Stomp {

class AngularCoordinate;           // class def. in stomp_angular_coordinate.h
class WeightedAngularCoordinate;   // class def. in stomp_angular_coordinate.h
class AngularCorrelation;          // class def. in stomp_angular_correlation.h
class Map;                         // class definition in stomp_map.h
class ScalarSubMap;
class ScalarMap;

typedef std::vector<ScalarSubMap> ScalarSubMapVector;
typedef ScalarSubMapVector::iterator ScalarSubMapIterator;
typedef std::pair<ScalarSubMapIterator, ScalarSubMapIterator> ScalarSubMapPair;

typedef std::vector<ScalarMap> ScalarMapVector;
typedef ScalarMapVector::iterator ScalarMapIterator;
typedef std::pair<ScalarMapIterator, ScalarMapIterator> ScalarMapPair;

class ScalarSubMap {
  // For the ScalarMap, the sub-map class is more about book-keeping
  // than operations.  They are still assembled around the superpixels, but
  // since the methods need to talk to pixels from different superpixels so
  // often it doesn't make sense to break up the data to the same degree.

 public:
  ScalarSubMap(uint32_t superpixnum);
  ~ScalarSubMap();
  void AddToArea(uint32_t resolution, double weight);
  void AddToIntensity(const double intensity, const uint32_t n_point = 1);
  void SetIntensity(const double intensity);
  void SetNPoints(const int n_point);
  double Area();
  double Intensity();
  int NPoints();
  double Density();
  double PointDensity();
  void SetBegin(ScalarIterator iter);
  void SetEnd(ScalarIterator iter);
  void SetNull(ScalarIterator iter);
  ScalarIterator Begin();
  ScalarIterator End();
  bool Initialized();
  uint32_t Size();

 private:
  uint32_t superpixnum_;
  ScalarIterator start_;
  ScalarIterator finish_;
  ScalarIterator null_;
  bool initialized_;
  double area_, total_intensity_;
  int total_points_;
  uint32_t size_;
};

class ScalarMap : public BaseMap {
  // Unlike a Map, where the set of Pixels is intended to match the
  // geometry of a particular region, ScalarMaps are intended to be
  // a regular sampling map of a given scalar field over some region.  The
  // area covered by the map will be approximately the same as that covered
  // by the pixels in the map, but each pixel is assumed to have some covering
  // fraction to indicate what percentage of the map is in the underlying
  // region.  To phrase things another way, once you have a Map describing
  // the extent of some data set, a ScalarMap is what you would use to
  // calculate clustering statistics on data contained in that region.

 public:
  // The ScalarPixels in a ScalarMap all have the same resolution.
  // Hence, we can construct a blank map from a Map (essentially
  // re-sampling the Map at a fixed resolution) or another ScalarMap,
  // so long as that map has higher resolution than the one we're trying to
  // construct.
  //
  // As described in the ScalarPixel class, there are three basic use cases:
  //
  // * Pure scalar field (e.g. CMB temperature or radio flux).
  // * Point-based density (e.g. the projected galaxy density over some area).
  // * Point-sampled field (e.g. the mean galaxy magnitude over some area).
  //
  // The way that we'll interact with the map will vary somewhat with each
  // case.  To make sure that we're doing the correct thing, we encode which
  // of these regimes we're operating under with the ScalarMapType enum:
  enum ScalarMapType {
    ScalarField,
    DensityField,
    SampledField
  };

  ScalarMap();

  // Initialize a ScalarMap based on the geometry of an input Map.  If the
  // use_map_weight_as_intensity flag is set to true, the MapType will be
  // set to ScalarField regardless of the value used in calling the
  // constructor.  A warning will be issued if the input value is not
  // "ScalarField".
  ScalarMap(Map& stomp_map,
	    uint32_t resolution,
	    ScalarMapType scalar_map_type = ScalarField,
	    double min_unmasked_fraction = 0.0000001,
	    bool use_map_weight_as_intensity = false);

  // If the map used to initialize the current one is also a ScalarMap, then
  // the ScalarMapType will be based on the input map, as will the geometry.
  ScalarMap(ScalarMap& scalar_map,
	    uint32_t resolution,
	    double min_unmasked_fraction = 0.0000001);

  // Initialize based on a vector of ScalarPixels.  If the input vector contains
  // pixels with heterogeneous resolutions, the code will exit automatically.
  ScalarMap(ScalarVector& pix,
	    ScalarMapType scalar_map_type = ScalarField,
	    double min_unmasked_fraction = 0.0000001);

  // This may seem a bit of an oddity, but needing roughly circular patches
  // from maps comes up more frequently than one might think.  Or not.
  ScalarMap(Map& stomp_map,
	    AngularCoordinate& center,
	    double theta_max,
	    uint32_t resolution,
	    ScalarMapType scalar_map_type = ScalarField,
	    double min_unmasked_fraction = 0.0000001,
	    double theta_min = -1.0);
  virtual ~ScalarMap();

  // Internal method that should probably never be called.
  bool _InitializeSubMap();

  // This is generally set through the constructor.  However, if
  // you want to re-initialize the same object with different parameters or
  // use the constructor without any arguments, this will set the
  // resolution of the map.
  void SetResolution(uint32_t resolution);

  // If you change the map resolution, then you will need to re-initialize
  // the coverage of the map from either a Map or a higher resolution
  // ScalarMap.  These methods allow for that functionality.  In each
  // case, if a resolution is supplied it will over-ride the current map's
  // resolution value.  This will also reset any previously set region
  // information.
  void InitializeFromMap(Map& stomp_map, uint32_t resolution=0,
			 bool use_map_weight_as_intensity = false);
  void InitializeFromScalarMap(ScalarMap& scalar_map, uint32_t resolution=0);

  // If you re-initialize from a vector of ScalarPixels, the resolution of
  // those pixels will automatically over-ride the current map's resolution
  // value.
  void InitializeFromScalarPixels(ScalarVector& pix,
				  ScalarMapType map_type = ScalarField);

  // Once we have our map set up, we'll want to add data points to it.  This
  // method offers two variations on that task.  If the MapType is ScalarField,
  // then the corresponding pixel will take on the value of the weight attached
  // to the input object.  Hence, adding another object which is located in
  // the same pixel will over-ride the old weight value with the new one.  In
  // all cases, the return value is false if the object doesn't localize to any
  // pixel in the map.
  bool AddToMap(AngularCoordinate& ang, double object_weight = 1.0);
  bool AddToMap(WeightedAngularCoordinate& ang);

  // Alternatively, if we are encoding a pure scalar field, then this method
  // will import the weight value from the input pixel into the proper fields.
  // If the input pixel is at a higher resolution than the current resolution
  // of the ScalarMap or the ScalarMap type is not ScalarField, the return
  // value is false.  For wholesale sampling from a Map, use InitializeFromMap.
  bool AddToMap(Pixel& pix);

  // Finally, if we want to set the map values directly, we can add a
  // ScalarPixel to the map.  If the resolution of the pixel doesn't match the
  // map resolution or the pixel can't be found in the map, the return value is
  // false.
  bool AddToMap(ScalarPixel& pix);

  // The instance of the same methods from the BaseMap class.  Note that, if
  // the input pixel is at a higher resolution than the current map,
  // FindUnmaskedFraction will return an invalid value (-1.0).  Likewise,
  // Coverage will over-ride the input resolution value if it is higher than
  // the map resolution.  FindUnmaskedStatus will work for higher resolution
  // pixels, but values indicating that the pixel is within the map don't
  // necessarily mean the same as they do for the Map class.
  virtual void Coverage(PixelVector& superpix,
			uint32_t resolution = HPixResolution);
  bool Covering(Map& stomp_map, uint32_t maximum_pixels);
  virtual double FindUnmaskedFraction(Pixel& pix);
  virtual int8_t FindUnmaskedStatus(Pixel& pix);

  // If we're converting a map from high to low resolution, this method
  // re-calculates the weight and intensity parameters for a given lower
  // resolution pixel.  Exactly how this is done will depend on the
  // ScalarMapType.  For a ScalarField, the resampling will be an
  // area-weighted average.  For DensityField and SampledField, the resampling
  // will be a direct sum of the intensity and points.  Resample will also
  // handle the case where the map has been converted to an over-density (the
  // resampled pixel will be based on the raw values).
  void Resample(ScalarPixel& pix);

  // When we access the field properties of the map, there are three options.
  // We can return the average intensity for the input pixel, the "density"
  // (intensity over the unmasked area), or the "point density" (number of
  // points over the unmasked area).  Which of these three will be most
  // meaningful will depend on the map type.
  double FindIntensity(Pixel& pix);
  double FindDensity(Pixel& pix);
  double FindPointDensity(Pixel& pix);

  // These four methods allow you to sample the area, intensity and density
  // within a given annular area around a point (where the inner radius is
  // set to zero by default.  The angular bounds (theta_min and theta_max) are
  // taken to be in degrees.
  double FindLocalArea(AngularCoordinate& ang, double theta_max,
		       double theta_min = -1.0);
  double FindLocalIntensity(AngularCoordinate& ang, double theta_max,
			    double theta_min = -1.0);
  double FindLocalDensity(AngularCoordinate& ang, double theta_max,
			  double theta_min = -1.0);
  double FindLocalPointDensity(AngularCoordinate& ang, double theta_max,
			       double theta_min = 0.0);

  // In the case where one is adding data points to the map, once this is
  // done, those object counts will need to be translated into a local measure
  // of the fractional over-density.  If the mean density of the map is all
  // that's required, then the first method will do that.  If you want to
  // replace the current data counts with the fractional over-density the
  // second method will do that (and call the first method if you have not
  // already done so).  And if you need to convert back from over-density to
  // raw values, the final method will do that for you.
  void CalculateMeanIntensity();
  void ConvertToOverDensity();
  void ConvertFromOverDensity();

  // Given a Map, we export the scalar field in the current
  // ScalarMap into its weight values.  This will naturally perform an
  // intersection between the areas covered by the two maps.  If there is no
  // overlapping area, then the method returns false and the input Map
  // will not be modified.  A true result means that there was at least some
  // overlap.
  bool ImprintMap(Map& stomp_map);

  // This method calculates the auto-correlation of the scalar field in the
  // current map.  The simpler case is the one where the argument is an
  // iterator for an angular bin.  Given the angular extent of that bin, the
  // code will find the auto-correlation of the field.  If the second option
  // is used, then the code will find the auto-correlation for all of the
  // angular bins whose resolution values match that of the current map.
  void AutoCorrelate(ThetaIterator theta_iter);
  void AutoCorrelate(AngularCorrelation& wtheta);

  // Alternatively, we may want to do the auto-correlation using the regions
  // we've divided the area into as a set of jack-knife samples.  These variants
  // will measure the auto-correlation for each of the jack-knife samples
  // simultaneously.
  void AutoCorrelateWithRegions(AngularCorrelation& wtheta);
  void AutoCorrelateWithRegions(ThetaIterator theta_iter);

  // Same as the auto-correlation methods, except the current map is
  // cross-correlated with another DensityMap.  Only areas of overlap are
  // considered in the cross-correlation.
  void CrossCorrelate(ScalarMap& scalar_map, AngularCorrelation& wtheta);
  void CrossCorrelate(ScalarMap& scalar_map, ThetaIterator theta_iter);
  void CrossCorrelateWithRegions(ScalarMap& scalar_map,
				 AngularCorrelation& wtheta);
  void CrossCorrelateWithRegions(ScalarMap& scalar_map,
				 ThetaIterator theta_iter);

  // Meaningful, since all of the pixels in the map share a common resolution.
  uint32_t Resolution();

  // Some global methods for accessing the aggregate area, intensity, density
  // and so for for the map.  If a superpixnum index is given as an
  // argument, then the values for that pixel are returned.  Otherwise, the
  // values for the whole map are given.  The extent to which these values
  // are meaningful will depend on the type of scalar field encoded in the map.
  double Intensity();
  double Intensity(uint32_t superpixnum);
  int NPoints();
  int NPoints(uint32_t superpixnum);
  double Density();
  double Density(uint32_t superpixnum);
  double PointDensity();
  double PointDensity(uint32_t superpixnum);
  ScalarIterator Begin(uint32_t superpixnum = MaxSuperpixnum);
  ScalarIterator End(uint32_t superpixnum = MaxSuperpixnum);
  double MeanIntensity();
  bool IsOverDensityMap();
  ScalarMapType MapType();

  // We need these methods to comply with the BaseMap signature.
  virtual double Area();
  double Area(uint32_t superpixnum);
  virtual uint32_t Size();
  uint32_t Size(uint32_t superpixnum);
  virtual uint32_t MinResolution();
  virtual uint32_t MaxResolution();
  virtual uint8_t MinLevel();
  virtual uint8_t MaxLevel();
  virtual bool Empty();
  virtual void Clear();


 private:
  ScalarVector pix_;
  ScalarSubMapVector sub_map_;
  ScalarMapType map_type_;
  double area_, mean_intensity_, unmasked_fraction_minimum_, total_intensity_;
  uint32_t resolution_, total_points_;
  bool converted_to_overdensity_, calculated_mean_intensity_;
  bool initialized_sub_map_;
};

} // end namespace Stomp

#endif
