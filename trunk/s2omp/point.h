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
#include "core.h"

// From S2
#include "s2.h"
#include "s2cellid.h"

namespace s2omp {

class pixel; // class declaration in pixel.h
class point;

typedef std::vector<point *> point_ptr_vector;
typedef point_ptr_vector::iterator point_ptr_iterator;

class point {
  // Our generic class for handling angular positions.  The idea is that
  // locations on the celestial sphere should be abstract objects from which
  // you can draw whatever angular coordinate pair is necessary for a given
  // use case.  point's can be instantiated with a particular
  // coordinate system in mind or that can be set later on.
public:
  enum Sphere {
    EQUATORIAL, GALACTIC, ECLIPTIC
  };

public:
  point();
  point(double x, double y, double z, double weight);
  virtual ~point() {};

  static point* from_latlon_deg(double lat_deg, double lon_deg, Sphere s);
  static point* from_latlon_deg(double lat_deg, double lon_deg, Sphere s,
      double weight);
  static point* from_latlon_rad(double lat_rad, double lon_rad, Sphere s);
  static point* from_latlon_rad(double lat_rad, double lon_rad, Sphere s,
      double weight);

  static point* from_radec_deg(double ra_deg, double dec_deg);
  static point* from_radec_deg(double ra_deg, double dec_deg, double weight);
  static point* from_radec_rad(double ra_rad, double dec_rad);
  static point* from_radec_rad(double ra_rad, double dec_rad, double weight);

  double lat_deg(Sphere s) const;
  double lon_deg(Sphere s) const;
  double lat_rad(Sphere s) const;
  double lon_rad(Sphere s) const;

  double ra_deg() const;
  double dec_deg() const;
  double ra_rad() const;
  double dec_rad() const;

  inline double unit_sphere_x() const {
    return point_.x();
  }
  inline double unit_sphere_y() const {
    return point_.y();
  }
  inline double unit_sphere_z() const {
    return point_.z();
  }

  inline double weight() const {
    return weight_;
  }
  inline void set_weight(double weight) {
    weight_ = weight;
  }

  double angular_distance(const point& p) const;
  static double angular_distance(const point& a, const point& b);

  inline double dot(const point& p) const;
  static double dot(const point& a, const point& b);

  point cross(const point& p) const;
  static point cross(const point& a, const point& b);

  point great_circle(const point& p);
  static point great_circle(const point& a, const point& b);

  double position_angle(const point& p, Sphere s);
  static double position_angle(const point& center, const point& target,
      Sphere s);

  void rotate_about(const point& axis, double rotation_angle_degrees);
  void rotate_about(const point& axis, double rotation_angle_degrees, Sphere s);
  static point rotate_about(const point& p, const point& axis,
      double rotation_angle, Sphere s);

  inline uint64 id() const;
  inline uint64 id(int level) const;
  pixel* to_pixel() const;
  pixel* to_pixel(int level) const;

  inline S2Point s2point() {return point_;}

  inline point operator -() const;

protected:
  void set_latlon_degrees(double lat_deg, double lon_deg, Sphere s);
  void set_latlon_radians(double lat_rad, double lon_rad, Sphere s);
  void set_xyz(double x, double y, double z);

private:

  double cos_position_angle(point& p, Sphere s);
  double sin_position_angle(point& p, Sphere s);
  void rotate_about(point& axis, double rotation_angle_degrees, Sphere s,
      double& unit_sphere_x, double& unit_sphere_y, double& unit_sphere_z);

  S2Point point_;
  double weight_;
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

inline point point::operator-() const{
  return point(-unit_sphere_x(), -unit_sphere_y(),
      -unit_sphere_z(), weight());
}

inline double point::dot(const point& p) const {
  return point_.x() * p.unit_sphere_x() + point_.y() * p.unit_sphere_y()
      + point_.z() * p.unit_sphere_z();
}

inline uint64 point::id() const {
  return S2CellId::FromPoint(point_).id();
}

inline uint64 point::id(int level) const {
  return S2CellId::FromPoint(point_).parent(level).id();
}

} // end namespace s2omp
#endif /* POINT_H_ */
