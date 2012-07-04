// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the abstract BaseMap class that serves as the
// basis for all of the *Map objects.  BaseMap sets out the basic functionality
// that all of the *Map classes need to describe a given region on the sky and
// do some basic internal maintenance.  Additionally, BaseMap provides a
// common set of methods for dividing that area up into nearly equal-area,
// similarly-shaped regions.  This functionality is the basis for calculating
// jack-knife errors for our various statistical analyses.

#ifndef ANALYTIC_BOUND_H_
#define ANALYTIC_BOUND_H_


typedef std::map<const uint64, int16_t> region_dict;
typedef region_dict::iterator region_iterator;
typedef std::pair<region_iterator, region_iterator> region_pair;
typedef std::map<const int16_t, double> region_area;

namespace s2omp {
  class region_map {
    // This class provides the functionality for dividing the area subtended by a
    // BaseMap-derived object into roughly equal-area, equal-sized regions.  The
    // class is not intended to be instantiated outside of the BaseMap class.

   public:
    region_map();
    virtual ~region_map();

    uint16_t initialize(analytic_bound* bound, uint16_t n_region);
    uint16_t initialize(analytic_bound* bound, uint16_t n_region, uint32_t level);

    bool initialize(analytic_bound* bound, analytic_bound& reference_bound);

    int16_t find_region(point& p);

    int16_t find_region(pixel& pix);

    void clear();

    void region_covering(int16_t region, pixel_vector& pixels);

    // Given a region index, return the area associated with that region.
    double region_area(int16_t region);

    // Some getter methods to describe the state of the RegionMap.
    uint16_t n_region();
    uint32_t level();
    bool initialized();

    // Return iterators for the set of RegionMap objects.
    region_iterator begin();
    region_iterator end();

   private:
    int find_regionation_level(analytic_bound* bound, uint16_t n_region);

    void regionate(analytic_bound* bound, uint16_t n_region, int level);

    bool verify_regionation(uint16_t n_region);

    region_dict region_map_;
    region_area region_area_;
    int level_;
    uint16_t n_region_;
  };

class pixelized_bound {
public:
  virtual bool is_empty();
  virtual long size();
  virtual void clear();
  virtual void area();

  virtual bool contains(point& p);
  virtual bool contains(pixel& pix);

  virtual double contained_area(pixel& pix);
  virtual bool may_intersect(pixel& pix);

  virtual void covering(pixel_vector& pixels);
  virtual void covering(int max_pixels, pixel_vector& pixels);
  virtual void simple_covering(int level, pixel_vector& pixels);

  virtual circle_bound get_bound();

  virtual point get_random_point();
  virtual void get_random_points(long n_points, pixel_vector& points);

  // These methods all act as wrappers for the RegionMapper object contained
  // in the class.  See that class for documentation.
  uint16_t initialize_regions(uint16_t n_regions);
  uint16_t initialize_regions(uint16_t n_regions, int region_level);
  bool initialize_regions(analytic_bound& bound);
  int16_t find_region(point& p);
  int16_t find_region(pixel& pix);
  void clear_regions();
  void region_covering(int16_t region, pixel_vector& pixels);
  double region_area(int16_t region);
  uint16_t n_region();
  uint32_t region_level();
  bool regions_initialized();
  region_iterator region_begin();
  region_terator region_end();

private:
  pixelized_bound();
  virtual ~pixelized_bound();
  region_map region_map_;
};

} // end namespace s2omp

#endif /* ANALYTIC_BOUND_H_ */
