// Copyright 2012  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the classes used to store point-like data on the
// sphere.  This can be simply a location on the sky (the point
// class) or a location along with additional information about the object at
// that position (Weightedpoint).  CosmoCoordinate extends the
// latter functionality to also allow for treatment of cosmological coordinates.


#ifndef POINT_H_
#define POINT_H_

#include <math.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

// From S2
#include <s2/s2.h>
#include <s2/s2cellid.h>
#include <s2/s2latlng.h>

#include "core.h"

namespace s2omp {

class pixel; // class declaration in pixel.h
class point;

typedef std::vector<point *> point_ptr_vector;
typedef point_ptr_vector::const_iterator point_ptr_iterator;

// Our generic class for handling angular positions.  The idea is that
// locations on the celestial sphere should be abstract objects from which
// you can draw whatever angular coordinate pair is necessary for a given
// use case.  points can be instantiated with a particular
// coordinate system in mind or that can be set later on.
class point {

public:
  // Public enum for defining the coordinate systems that point supports
  enum Coordinate {
    EQUATORIAL, GALACTIC, ECLIPTIC
  };

public:
  point();
  point(double x, double y, double z, double weight);
  point(S2Point s2point, double weight); // Does not normalize the input point.
  point(S2LatLng s2latlng, double weight);
  virtual ~point();

  static point from_latlon_deg(double lat_deg, double lon_deg, Coordinate c);
  static point from_latlon_deg(double lat_deg, double lon_deg, Coordinate c,
      double weight);
  static point from_latlon_rad(double lat_rad, double lon_rad, Coordinate c);
  static point from_latlon_rad(double lat_rad, double lon_rad, Coordinate c,
      double weight);

  static point from_radec_deg(double ra_deg, double dec_deg);
  static point from_radec_deg(double ra_deg, double dec_deg, double weight);
  static point from_radec_rad(double ra_rad, double dec_rad);
  static point from_radec_rad(double ra_rad, double dec_rad, double weight);

  static point* copy_point(const point& p);
  static point copy_point(point* p);

  double lat_deg(Coordinate c) const;
  double lon_deg(Coordinate c) const;
  double lat_rad(Coordinate c) const;
  double lon_rad(Coordinate c) const;

  double ra_deg() const;
  double dec_deg() const;
  double ra_rad() const;
  double dec_rad() const;

  static void equatorial_to_galactic(double ra_rad, double dec_rad,
      double* glon_rad, double* glat_rad);
  static void galactic_to_equatorial(double glon_rad, double glat_rad,
      double* ra_rad, double* dec_rad);
  static void equatorial_to_ecliptic(double ra_rad, double dec_rad,
      double* elon_rad, double* elat_rad);
  static void ecliptic_to_equatorial(double elon_rad, double elat_rad,
      double* ra_rad, double* dec_rad);
  static void galactic_to_ecliptic(double glon_rad, double glat_rad,
      double* elon_rad, double* elat_rad);
  static void ecliptic_to_galactic(double elon_rad, double elat_rad,
      double* glon_rad, double* glat_rad);

  bool is_normalized() const;

  void normalize();

  inline double x() const {
    return point_.x();
  }
  inline double y() const {
    return point_.y();
  }
  inline double z() const {
    return point_.z();
  }

  double unit_sphere_x();
  double unit_sphere_y();
  double unit_sphere_z();

  double unit_sphere_x(Coordinate c);
  double unit_sphere_y(Coordinate c);
  double unit_sphere_z(Coordinate c);

  inline double weight() const {
    return weight_;
  }
  inline void set_weight(double weight) {
    weight_ = weight;
  }

  inline double dot(const point& p) const;
  static double dot(const point& a, const point& b);

  // For cross products, we have two options.  We can either return the point
  // that's in the direction of A x B but which will be a unit normal vector,
  // or we can return the amplitude of A x B, which is |A|B| sin (theta).  For
  // the former, we have the cross() method and for the latter, cross_norm().
  point cross(const point& p) const;
  static point cross(const point& a, const point& b);
  double cross_norm(const point& p) const;
  static double cross_norm(const point& a, const point& b);
  double cross_norm2(const point& p) const;
  static double cross_norm2(const point& a, const point& b);

  // To get an accurate angular distance at small scales, we use both dot()
  // and cross_norm().
  double angular_distance_deg(const point& p) const;
  static double angular_distance_deg(const point& a, const point& b);
  double angular_distance_rad(const point& p) const;
  static double angular_distance_rad(const point& a, const point& b);

  point great_circle(const point& p) const;
  static point great_circle(const point& a, const point& b);

  double position_angle(const point& p, Coordinate c) const;
  static double position_angle(const point& center, const point& target,
      Coordinate c);

  void rotate_about(const point& axis, double rotation_angle_degrees);
  void rotate_about(const point& axis, double rotation_angle_degrees,
      Coordinate c);
  static point rotate_about(const point& p, const point& axis,
      double rotation_angle_degrees, Coordinate c);

  inline uint64 id() const;
  inline uint64 id(int level) const;
  static inline uint64 point_to_id(const point& p);
  static inline uint64 point_to_id(const point& p, int level);
  static inline point s2point_to_point(const S2Point& p);
  static inline S2Point point_to_s2point(const point& p);
  pixel to_pixel() const;
  pixel to_pixel(int level) const;

  inline S2Point s2point() const {
    return point_;
  }

  inline point operator-() const;

protected:
  void set_latlon_degrees(double lat_deg, double lon_deg, Coordinate c);
  void set_latlon_radians(double lat_rad, double lon_rad, Coordinate c);
  void set_xyz(double x, double y, double z) {
    double norm = sqrt(x * x + y * y + z * z);
    point_ = S2Point(x / norm, y / norm, z / norm);
  }

private:
  static double cos_position_angle(const point& center, const point& target,
      Coordinate c);
  static double sin_position_angle(const point& center, const point& target,
      Coordinate c);

  S2Point point_;
  double weight_;
  static double const NORMALIZATION_PRECISION = 1.0e-6;
};

inline bool operator==(point const& a, point const& b) {
  return a.id() == b.id();
}

inline bool operator!=(point const& a, point const& b) {
  return a.id() != b.id();
}

inline bool operator<(point const& a, point const& b) {
  return a.id() < b.id();
}

inline bool operator>(point const& a, point const& b) {
  return a.id() > b.id();
}

inline bool operator<=(point const& a, point const& b) {
  return a.id() <= b.id();
}

inline bool operator>=(point const& a, point const& b) {
  return a.id() >= b.id();
}

inline point point::operator-() const {
  return point(-point_, weight_);
}

inline double point::dot(const point& p) const {
  return point_.x() * p.x() + point_.y() * p.y() + point_.z() * p.z();
}

inline uint64 point::id() const {
  return S2CellId::FromPoint(point_).id();
}

inline uint64 point::id(int level) const {
  return S2CellId::FromPoint(point_).parent(level).id();
}

inline uint64 point::point_to_id(const point& p) {
  return S2CellId::FromPoint(p.s2point()).id();
}

inline uint64 point::point_to_id(const point& p, int level) {
  return S2CellId::FromPoint(p.s2point()).parent(level).id();
}

inline point point::s2point_to_point(const S2Point& p) {
  return point(p, 1.0);
}

inline S2Point point::point_to_s2point(const point& p) {
  return S2Point(p.x(), p.y(), p.z());
}

} // end namespace s2omp
#endif /* POINT_H_ */
