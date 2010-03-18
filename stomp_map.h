// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the Map class.  Maps are intended to describe an
// arbitrary area on the sky as precisely as possible, given the limits of
// pixel resolution, physical memory, etc.  Internally, this information is
// encoded using pixels at a range of resolutions, depending on how the
// boundaries of the area in question and the pixelization scheme interact.
// However, the goal of the class is to abstract away those details, allowing
// the user to treat Maps as a pure representative of spherical geometry.

#ifndef STOMP_MAP_H
#define STOMP_MAP_H

#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"
#include "stomp_base_map.h"

namespace Stomp {

class AngularCoordinate;  // class declaration in stomp_angular_coordinate.h
class Pixel;              // class declaration in stomp_pixel.h
class GeometricBound;     // class declaration in stomp_geometry.h
class SubMap;
class Map;

typedef std::map<const uint32_t, uint32_t> ResolutionDict;
typedef ResolutionDict::iterator ResolutionIterator;
typedef std::pair<ResolutionIterator, ResolutionIterator> ResolutionPair;

typedef std::vector<SubMap> SubMapVector;
typedef SubMapVector::iterator SubMapIterator;
typedef std::pair<SubMapIterator, SubMapIterator> SubMapPair;

typedef std::pair<uint32_t, PixelIterator> MapIterator;
typedef std::pair<MapIterator, MapIterator> MapPair;

class SubMap {
  // While the preferred interface for interacting with a Map is through
  // that class, the actual work is deferred to the SubMap class.  Each
  // instance contains all of the pixels for a corresponding superpixel as well
  // as some summary statistics.  All of the operations done on the Map
  // end up calling corresponding methods for each of the SubMap instances
  // contained in that Map.  See the comments around those methods in the
  // Map class declaration for an explaination of what each method does.

 public:
  SubMap(uint32_t superpixnum);
  ~SubMap();
  void AddPixel(Pixel& pix);
  void Resolve(bool force_resolve = false);
  void SetMinimumWeight(double minimum_weight);
  void SetMaximumWeight(double maximum_weight);
  void SetMaximumResolution(uint32_t maximum_resolution, bool average_weights);
  bool FindLocation(AngularCoordinate& ang, double& weight);
  double FindUnmaskedFraction(Pixel& pix);
  int8_t FindUnmaskedStatus(Pixel& pix);
  double FindAverageWeight(Pixel& pix);
  void FindMatchingPixels(Pixel& pix, PixelVector& match_pix,
			  bool use_local_weights = false);
  double AverageWeight();
  void Soften(PixelVector& softened_pix, uint32_t maximum_resolution,
	      bool average_weights);
  bool Add(Map& stomp_map, bool drop_single);
  bool Multiply(Map& stomp_map, bool drop_single);
  bool Exclude(Map& stomp_map);
  void ScaleWeight(const double weight_scale);
  void AddConstantWeight(const double add_weight);
  void InvertWeight();
  void Pixels(PixelVector& pix);
  void CheckResolution(uint32_t resolution);
  void Clear();
  uint32_t Superpixnum();
  PixelIterator Begin();
  PixelIterator End();
  double Area();
  bool Initialized();
  bool Unsorted();
  uint32_t MinResolution();
  uint32_t MaxResolution();
  uint8_t MinLevel();
  uint8_t MaxLevel();
  double MinWeight();
  double MaxWeight();
  double LambdaMin();
  double LambdaMax();
  double EtaMin();
  double EtaMax();
  double ZMin();
  double ZMax();
  uint32_t Size();
  uint32_t PixelCount(uint32_t resolution);

 private:
  uint32_t superpixnum_, size_;
  PixelVector pix_;
  double area_, lambda_min_, lambda_max_, eta_min_, eta_max_, z_min_, z_max_;
  double min_weight_, max_weight_;
  uint8_t min_level_, max_level_;
  bool initialized_, unsorted_;
  ResolutionDict pixel_count_;
};

class Map : public BaseMap {
  // A Map is intended to function as a region on the sky whose geometry
  // is given by a set of Pixels of various resolutions which combine to
  // cover that area.  Since each Pixel has an associated weight, the map
  // can also encode a scalar field (temperature, observing depth, local
  // seeing, etc.) over that region.  A Map can be combined with other
  // Maps with all of the logical operators you would expect (union,
  // intersection and exclusion as well as addition and multiplication of the
  // weights as a function of position).  Likewise, you can test angular
  // positions and pixels against a Map to see if they are within it or
  // query the angular extent and area of a Map on the Sky.

 public:
  // The preferred constructor for a Map takes a vector of Pixels
  // as its argument.  However, it can be constructed from a properly formatted
  // ASCII text file as well.
  Map();
  Map(PixelVector& pix, bool force_resolve = true);
  Map(std::string& InputFile,
      const bool hpixel_format = true,
      const bool weighted_map = true);

  // Alternatively, we can specify a Map's geometry using a GeometricBound.
  // This object provides an analytic description of a spatial region that
  // the Map then approximates using uniform weight Pixels.  The fidelity of
  // the pixelization is determined by the maximum resolution specified in the
  // constructor.  Likewise, since this is an approximation, a boolean flag is
  // included to provide feedback on whether the pixelization was successful
  // or not and how the area of the resulting Map area compares to the input
  // GeometricBound's area.
  Map(GeometricBound& bound, double weight = 1.0,
      uint32_t maximum_resolution = MaxPixelResolution,
      bool verbose = false);
  virtual ~Map();

  // Initialize is called to organize the Map internally.  Unless the
  // map is being reset with a new set of pixels, as in the second instance of
  // this method, Initialize should probably never be invoked.
  bool Initialize();
  bool Initialize(PixelVector& pix, bool force_resolve = true);

  // Alternatively, we can add pixels one at a time.  If this is done in such
  // a way that the input pixels are not ordered in the same manner as they
  // are stored in the Map object, then you will need to call Initialize()
  // before most of the Map methods will work properly.
  void AddPixel(Pixel& pix);

  // Simple call to determine if a point is within the current Map
  // instance.  Returns true if the point is within the map; false, otherwise.
  bool FindLocation(AngularCoordinate& ang);

  // An variation on FindLocation that also assigns a value to the input
  // "weight" reference variable.  If the location is within the Map, then
  // the weight of value of the map is stored in the "weight" variable.  If
  // not, then the value in "weight" is meaningless.
  bool FindLocation(AngularCoordinate& ang, double& weight);

  // Another variation on FindLocation, this is mostly in place for the Python
  // wrapper.  Instead of returning a boolean, this returns the weight value
  // for the Map at the input location.  If the location is not within
  // the map, then the default value of -1.0e30 is returned.
  double FindLocationWeight(AngularCoordinate& ang);

  // In the same spirit, we can pose similar queries for both Pixels and Maps.
  // Later on, we'll have more sophisticated indicators as to whether these
  // areas are fully, partially or not contained in our Map, but for now we only
  // consider the question of full containment.
  bool Contains(Pixel& pix);
  bool Contains(Map& stomp_map);

  // Given a Pixel, this returns the fraction of that pixel's area that is
  // contained within the current map (0 <= fraction <= 1).  Alternatively, a
  // vector of pixels can be processed in a single call, in which case a
  // vector of coverages is returned or the unmasked fraction is stored in the
  // weight element of the vector of pixels.
  virtual double FindUnmaskedFraction(Pixel& pix);
  void FindUnmaskedFraction(PixelVector& pix,
                            std::vector<double>& unmasked_fraction);
  void FindUnmaskedFraction(PixelVector& pix);
  double FindUnmaskedFraction(Map& stomp_map);

  // Similar to FindUnmaskedFraction, these routines return an integer
  // indicating whether the pixel is fully contained in the map (1),
  // fully outside the map (0) or partially contained in the map (-1).  For
  // cases where we don't have to come up with an accurate number for the
  // map coverage of a given pixel, this should be faster.
  virtual int8_t FindUnmaskedStatus(Pixel& pix);
  void FindUnmaskedStatus(PixelVector& pix,
			  std::vector<int8_t>& unmasked_status);
  int8_t FindUnmaskedStatus(Map& stomp_map);

  // Similar to FindUnmaskedFraction, this returns the area-averaged weight of
  // the map over the area covered by the input pixel (or pixels).
  // AverageWeight does the same task, but over the entire Map.
  double FindAverageWeight(Pixel& pix);
  void FindAverageWeight(PixelVector& pix,
                         std::vector<double>& average_weight);
  void FindAverageWeight(PixelVector& pix);
  double AverageWeight();

  // This is part of the process for finding the intersection between two maps.
  // For a given pixel, we return the pixels in our map that are contained
  // within that test pixel.  If use_local_weights is set to true, then the
  // pixel weights are set to match the weights in the current map.  If not,
  // then the matching pixels are set to the weight from the test Pixel.
  void FindMatchingPixels(Pixel& pix,
			  PixelVector& match_pix,
			  bool use_local_weights = false);
  void FindMatchingPixels(PixelVector& pix,
			  PixelVector& match_pix,
			  bool use_local_weights = false);

  // Return a vector of SuperPixels that cover the Map.  This serves two
  // purposes.  First, it acts as a rough proxy for the area of the current
  // map, which can occasionally be useful.  More importantly, all of the real
  // work in a Map is done on a superpixel-by-superpixel basis, so this
  // becomes an important thing to know when querying the map.
  //
  // If the resolution argument is given, then the resulting set of pixels will
  // be generated at that resolution instead.  In either case, the weight for
  // the pixels will reflect the fraction of the pixel that is within the
  // current map.
  virtual void Coverage(PixelVector& superpix,
			uint32_t resolution = Stomp::HPixResolution);

  // Instead of a set of vectors at the same resolution, we may want to
  // generate a lower resolution version of our current map where the needs of
  // matching the geometry of a given region on the sky can be compromised a
  // bit in the interest of a smaller memory footprint.  In this case, the
  // return object is another Map of greater or equal area which is
  // composed of at most maximum_pixels pixels.  In some cases, the maximum
  // number of pixels is smaller than the number of superpixels in the current
  // map.  In that case, the returned boolean is false and the stomp_map is
  // composed of the superpixels for the current map.  If the method is able
  // to generate a map with at most maximum_pixels the return value is true.
  //
  // Unlike Coverage, the weights in the returned stomp_map will be based on
  // the area-averaged weights in the current Map.
  bool Covering(Map& stomp_map, uint32_t maximum_pixels);

  // As a middle ground between these two cases, we may want a variation of
  // Coverage that outputs a Map where we have reduced the maximum resolution
  // of our current map.  Again, the geometry will be less precise than the
  // current map, but the total number of pixels should be smaller (but not
  // explicitly specified as with Covering).  If the input resolution is larger
  // than the current maximum resolution, then the returned map will be a
  // copy of the current map.  If the average_weights flag is set to true, then
  // the resulting map will retain the current map's weight, averaging the
  // weight for any pixels which are resampled.  If set to false, the resulting
  // map will have unity weight, except for the pixels were resampled, where the
  // value will indicate the included fraction.
  void Soften(Map& stomp_map, uint32_t maximum_resolution,
	      bool average_weights=false);

  // Rather than creating a new Map, we can Soften the current Map
  void Soften(uint32_t maximum_resolution, bool average_weights=false);

  // In addition to modifying the maximum resolution of the Map (which is
  // basically what Soften does), we can also filter the current Map based on
  // the Weight, removing any area that violates the Weight limits.
  void SetMinimumWeight(double min_weight);
  void SetMaximumWeight(double max_weight);

  // After we initialize regions on this Map, we might want to produce a Map
  // corresponding to a given region.  The boolean here is to indicate if the
  // Map was properly constructed or not, depending on whether the specified
  // region index was within the valid range: 0 < region_index < n_region-1.
  bool RegionOnlyMap(int16_t region_index, Map& stomp_map);

  // Conversely, we often want to know what our map looks like when we've
  // excluded a specific region, as you'd do if you were using jack-knife error
  // estimates.  Similar operation as RegionMap with regards to the returned
  // boolean.
  bool RegionExcludedMap(int16_t region_index, Map& stomp_map);

  // Given a requested number of points, return a vector of Poisson random
  // angular positions within the current Map's area.
  //
  // If the use_weighted_sampling flag is set to true, then the local weight is
  // taken into account when generating random points.  In this case, a pixel
  // with the same area but twice the weight as another pixel should, in the
  // limit of infinite realizations, have twice as many points as the
  // lower-weighted one.
  void GenerateRandomPoints(AngularVector& ang,
                            uint32_t n_point = 1,
                            bool use_weighted_sampling = false);

  // Instead of a fixed number of random positions, we may have either a
  // vector of WeightedAngularCoordinate objects where we want to randomize
  // the positions or a vector of weights that need random angular positions
  // to go along with them.  In either case, the number of random positions
  // is taken from the size of the input vector.
  void GenerateRandomPoints(WAngularVector& ang, WAngularVector& input_ang);
  void GenerateRandomPoints(WAngularVector& ang, std::vector<double>& weights);

  // The book-end to the initialization method that takes an ASCII filename
  // as an argument, this method writes the current map to an ASCII file using
  // the same formatting conventions.
  bool Write(std::string& OutputFile, bool hpixel_format = true,
             bool weighted_map = true);

  // Alternatively, if the default constructor is used, this method will
  // re-initialize the map with the contents of the input file.
  bool Read(std::string& InputFile, const bool hpixel_format = true,
	    const bool weighted_map = true);

  // Another option for specifying the Map geometry is to use a GeometricBound
  // object.  This translates from the analytic region described in the
  // GeometricBound to a pixel-based version that we can use as a basis for a
  // Map.  The input weight value will be uniform over the Map and the
  // maximum_resolution value controls the level of Map fidelity to the
  // GeometricBound's area (although higher fidelity comes at the expense of
  // more pixels).  The return value indicates success or failure of the
  // translation.
  bool PixelizeBound(GeometricBound& bound, double weight = 1.0,
		     uint32_t maximum_resolution = MaxPixelResolution);

  // The pixelization method is iteratively adaptive.  First, it tries to find
  // the largest pixels that will likely fit inside the footprint.  Then it
  // checks those pixels against the bound, keeping the ones that are
  // fully inside.  The ones that were at least partially inside are refined
  // to the next resolution level and tested again.  This continues until we
  // reach the maximum resolution level, at which point we keep enough of the
  // partials to match the footprint's area, preferentially keeping those
  // pixels that are most contained by the bound.  These internal methods
  // handle this process.
  uint8_t _FindStartingResolutionLevel(double bound_area);
  bool _FindXYBounds(const uint8_t resolution_level,
		     GeometricBound& bound,
		     uint32_t& x_min, uint32_t& x_max,
		     uint32_t& y_min, uint32_t& y_max);
  double _ScorePixel(GeometricBound& bound, Pixel& pix);

  // Three simple functions for performing the same operation on the weights
  // of all of the pixels in the current map.  These are prelude to the next
  // set of functions for doing logical and arithmetic operations on Maps.
  void ScaleWeight(const double weight_scale);
  void AddConstantWeight(const double add_weight);
  void InvertWeight();

  // Now we begin the core of our class, the ability to treat the Map as
  // an abstract object which we can combine with other maps to form arbitrarily
  // complicated representations on the sphere.
  //
  // Starting simple, IngestMap simply takes the area associated with another
  // map and combines it with the current map.  If pixels overlap between the
  // two maps, then the weights are set to the average of the two maps.
  // Returns true if the procedure succeeded, false otherwise.
  bool IngestMap(PixelVector& pix, bool destroy_copy = true);
  bool IngestMap(Map& stomp_map, bool destroy_copy = true);

  // Now we have intersection.  This method finds the area of intersection
  // between the current map and the argument map and makes that the new area
  // for this map.  Weights are drawn from the current map's values.  If there
  // is no overlapping area, the method returns false and does nothing.  A
  // true response indicates that the area of the current map has changed to
  // the overlapping area between the two and the area is non-zero.
  bool IntersectMap(PixelVector& pix);
  bool IntersectMap(Map& stomp_map);

  // The inverse of IntersectMap, ExcludeMap removes the area associated with
  // the input map from the current map.  If this process would remove all of
  // the area from the current map, then the method returns false and does not
  // change the current map.  A true response indicates that the input maps has
  // been excluded and area remains.
  bool ExcludeMap(PixelVector& pix, bool destroy_copy = true);
  bool ExcludeMap(Map& stomp_map, bool destroy_copy = true);

  // Two sets of methods that operate on the weights between two different
  // maps.  The first set adds the weights of the two maps and the second
  // takes their product.  The drop_single boolean indicates whether the
  // non-overlapping area should be excluded (true) or not (false).  If
  // drop_single is set to false, then the areas where the two maps don't
  // overlap will have their weights set to whatever they are in the map that
  // covers that area.
  bool AddMap(PixelVector& pix, bool drop_single = true);
  bool AddMap(Map& stomp_map, bool drop_single = true);
  bool MultiplyMap(PixelVector& pix, bool drop_single = true);
  bool MultiplyMap(Map& stomp_map, bool drop_single = true);

  // Like IntersectMap, except that the current map takes on the weight
  // values from the map given as the argument.
  bool ImprintMap(PixelVector& pix);
  bool ImprintMap(Map& stomp_map);

  // Simple method for returning the vector representation of the current
  // Map.  If a superpixel index is given as the second argument, then
  // just the pixels for that superpixel are returned.
  void Pixels(PixelVector& pix,
              uint32_t superpixnum = Stomp::MaxSuperpixnum);

  // For a more efficient way to iterate through the pixels that make up the
  // Map, we have these methods.  The standard for loop to iterate through all
  // of the Pixels in a Map would be
  //
  // for (MapIterator iter=stomp_map.Begin();
  //      iter!=stomp_map.End();stomp_map.Iterate(&iter)) {
  //   double weight = iter.second->Weight();
  // }
  MapIterator Begin();
  MapIterator End();
  void Iterate(MapIterator* iter);

  // Resets the Map to a completely clean slate.  No pixels, no area.
  virtual void Clear();
  void Clear(uint32_t superpixnum);

  // Simple method for checking to see if the current map has any area
  // in a given superpixel.
  bool ContainsSuperpixel(uint32_t superpixnum);

  // Some general methods for querying the state of the current map.
  virtual double Area();
  double Area(uint32_t superpixnum);
  virtual uint32_t MinResolution();
  uint32_t MinResolution(uint32_t superpixnum);
  virtual uint32_t MaxResolution();
  uint32_t MaxResolution(uint32_t superpixnum);
  virtual uint8_t MinLevel();
  uint8_t MinLevel(uint32_t superpixnum);
  virtual uint8_t MaxLevel();
  uint8_t MaxLevel(uint32_t superpixnum);
  double MinWeight();
  double MinWeight(uint32_t superpixnum);
  double MaxWeight();
  double MaxWeight(uint32_t superpixnum);
  virtual uint32_t Size();
  uint32_t Size(uint32_t superpixnum);
  virtual bool Empty();
  uint32_t PixelCount(uint32_t resolution);


private:
  SubMapVector sub_map_;
  MapIterator begin_, end_;
  double area_, min_weight_, max_weight_;
  uint8_t min_level_, max_level_;
  uint32_t size_;
  ResolutionDict pixel_count_;
};

} // end namespace Stomp

#endif
