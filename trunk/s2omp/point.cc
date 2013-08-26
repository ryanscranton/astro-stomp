/*
 * point.cc
 *
 *  Created on: Aug 16, 2012
 *      Author: morrison
 */

#include "pixel.h"
#include "point.h"

namespace s2omp {

point::point() {
  point_ = S2Point(0, 0, 0);
  weight_ = 0;
}

point::point(double x, double y, double z, double weight) {
  weight_ = weight;
  double norm = sqrt(x * x + y * y + z * z);
  point_ = S2Point(x / norm, y / norm, z / norm);
}

point::point(S2Point point, double weight) {
  point_ = point;
  weight_ = weight;
}

point::~point() {
}

point point::from_latlon_deg(double lat_deg, double lon_deg, Coordinate c) {
  return from_latlon_rad(lat_deg * DEG_TO_RAD, lon_deg * DEG_TO_RAD, c);
}

point point::from_latlon_deg(double lat_deg, double lon_deg, Coordinate c,
    double weight) {
  return from_latlon_rad(lat_deg * DEG_TO_RAD, lon_deg * DEG_TO_RAD, c, weight);
}

point point::from_latlon_rad(double lat_rad, double lon_rad, Coordinate c) {
  return from_latlon_rad(lat_rad, lon_rad, c, 1.0);
}

point point::from_latlon_rad(double lat_rad, double lon_rad, Coordinate c,
    double weight) {
  double ra_rad = 0.0, dec_rad = 0.0;
  switch (c) {
  case EQUATORIAL:
    ra_rad = lon_rad;
    dec_rad = lat_rad;
    break;
  case GALACTIC:
    galactic_to_equatorial(lon_rad, lat_rad, &ra_rad, &dec_rad);
    break;
  case ECLIPTIC:
    galactic_to_ecliptic(lon_rad, lat_rad, &ra_rad, &dec_rad);
    break;
  }

  return point(sin(lat_rad) * cos(lon_rad), sin(lat_rad) * sin(lon_rad), cos(
      lat_rad), weight);
}

point point::from_radec_deg(double ra_deg, double dec_deg) {
  return from_latlon_rad(dec_deg * DEG_TO_RAD, ra_deg * DEG_TO_RAD, EQUATORIAL,
      1.0);
}

point point::from_radec_deg(double ra_deg, double dec_deg, double weight) {
  return from_latlon_rad(dec_deg * DEG_TO_RAD, ra_deg * DEG_TO_RAD, EQUATORIAL,
      weight);
}

point point::from_radec_rad(double ra_rad, double dec_rad) {
  return from_latlon_rad(dec_rad, ra_rad, EQUATORIAL, 1.0);
}

point point::from_radec_rad(double ra_rad, double dec_rad, double weight) {
  return from_latlon_rad(dec_rad, ra_rad, EQUATORIAL, weight);
}

point* point::copy_point(const point& p) {
  return new point(p.s2point(), p.weight());
}

point point::copy_point(point* p) {
  return point(p->s2point(), p->weight());
}

double point::lat_deg(Coordinate c) const {
  return RAD_TO_DEG * lat_rad(c);
}

double point::lon_deg(Coordinate c) const {
  return RAD_TO_DEG * lon_rad(c);
}

double point::lat_rad(Coordinate c) const {
  switch (c) {
  case EQUATORIAL:
    return acos(point_.z());
  case GALACTIC:
    double glon_rad, glat_rad;
    equatorial_to_galactic(ra_rad(), dec_rad(), &glon_rad, &glat_rad);
    return glat_rad;
  case ECLIPTIC:
    double elon_rad, elat_rad;
    equatorial_to_ecliptic(ra_rad(), dec_rad(), &elon_rad, &elat_rad);
    return elon_rad;
  }

  // Default to equatorial.
  return acos(point_.z());
}

double point::lon_rad(Coordinate c) const {
  switch (c) {
  case EQUATORIAL:
    return atan2(point_.y(), point_.x());
  case GALACTIC:
    double glon_rad, glat_rad;
    equatorial_to_galactic(ra_rad(), dec_rad(), &glon_rad, &glat_rad);
    return glon_rad;
  case ECLIPTIC:
    double elon_rad, elat_rad;
    equatorial_to_ecliptic(ra_rad(), dec_rad(), &elon_rad, &elat_rad);
    return elon_rad;
  }

  // Default to equatorial.
  return atan2(point_.y(), point_.x());
}

double point::ra_deg() const {
  return lon_deg(EQUATORIAL);
}

double point::dec_deg() const {
  return lat_deg(EQUATORIAL);
}

double point::ra_rad() const {
  return lon_rad(EQUATORIAL);
}

double point::dec_rad() const {
  return lat_rad(EQUATORIAL);
}

void point::equatorial_to_galactic(double ra_rad, double dec_rad,
    double* glon_rad, double* glat_rad) {
  static double const G_PSI = 0.57477043300;
  static double const STHETA = 0.88998808748;
  static double const CTHETA = 0.45598377618;
  static double const G_PHI = 4.9368292465;

  double a = ra_rad - G_PHI;
  double b = dec_rad;

  double sb = sin(b);
  double cb = cos(b);
  double cbsa = cb * sin(a);

  b = -1.0 * STHETA * cbsa + CTHETA * sb;
  if (b > 1.0) b = 1.0;

  double bo = asin(b);

  a = atan2(CTHETA * cbsa + STHETA * sb, cb * cos(a));

  double ao = (a + G_PSI + 4.0 * PI);

  while (ao > 2.0 * PI) ao -= 2.0 * PI;

  *glon_rad = ao;
  if (*glon_rad < 0.0) *glon_rad += 2.0 * PI;
  if (*glon_rad > 2.0 * PI) *glon_rad -= 2.0 * PI;

  *glat_rad = bo;
}

void point::galactic_to_equatorial(double glon_rad, double glat_rad,
    double* ra_rad, double* dec_rad) {
  static double const G_PSI = 4.9368292465;
  static double const STHETA = -0.88998808748;
  static double const CTHETA = 0.45598377618;
  static double const G_PHI = 0.57477043300;

  double a = glon_rad - G_PHI;
  double b = glat_rad;

  double sb = sin(b);
  double cb = cos(b);
  double cbsa = cb * sin(a);

  b = -1.0 * STHETA * cbsa + CTHETA * sb;
  if (b > 1.0) b = 1.0;

  double bo = asin(b);

  a = atan2(CTHETA * cbsa + STHETA * sb, cb * cos(a));

  double ao = (a + G_PSI + 4.0 * PI);
  while (ao > 2.0 * PI) ao -= 2.0 * PI;

  *ra_rad = ao;
  if (*ra_rad < 0.0) *ra_rad += 2.0 * PI;
  if (*ra_rad > 2.0 * PI) *ra_rad -= 2.0 * PI;

  *dec_rad = bo;
}

void point::equatorial_to_ecliptic(double ra_rad, double dec_rad,
    double* elon_rad, double* elat_rad) {
  double phi = 0.5 * PI - dec_rad;
  double x_eq = cos(ra_rad) * cos(phi);
  double y_eq = sin(ra_rad) * cos(phi);
  double z_eq = sin(phi);

  double x_ec = x_eq;
  double y_ec = COS_OBLIQUITY * y_eq + SIN_OBLIQUITY * z_eq;
  double z_ec = COS_OBLIQUITY * z_eq - SIN_OBLIQUITY * y_eq;

  *elon_rad = atan2(y_ec, x_ec);
  *elat_rad = acos(z_ec);
}

void point::ecliptic_to_equatorial(double elon_rad, double elat_rad,
    double* ra_rad, double* dec_rad) {
  double phi = 0.5 * PI - elat_rad;
  double x_ec = cos(elon_rad) * cos(phi);
  double y_ec = sin(elon_rad) * cos(phi);
  double z_ec = sin(phi);

  double x_eq = x_ec;
  double y_eq = COS_OBLIQUITY * y_ec - SIN_OBLIQUITY * z_ec;
  double z_eq = COS_OBLIQUITY * z_ec + SIN_OBLIQUITY * y_ec;

  *ra_rad = atan2(y_eq, x_eq);
  *dec_rad = acos(z_eq);
}

double point::dot(const point& a, const point& b) {
  return a.unit_sphere_x() * b.unit_sphere_x() + a.unit_sphere_y()
      * b.unit_sphere_y() + a.unit_sphere_z() * b.unit_sphere_z();
}

double point::angular_distance(const point& p) const {
  return acos(dot(p));
}

double point::angular_distance(const point& a, const point& b) {
  return acos(dot(a, b));
}

// TODO(cbmorrison) Not sure what the rules should be for returning the weight
// on the cross product of the points. Currently I return the weight of the
// first point.
point point::cross(const point& p) const {
  double x = point_.y() * p.unit_sphere_z() - point_.z() * p.unit_sphere_y();
  double y = point_.z() * p.unit_sphere_x() - point_.x() * p.unit_sphere_z();
  double z = point_.x() * p.unit_sphere_y() - point_.y() * p.unit_sphere_x();

  return point(x, y, z, weight_);
}

point point::cross(const point& a, const point& b) {
  double x = a.unit_sphere_y() * b.unit_sphere_z() - a.unit_sphere_z()
      * b.unit_sphere_y();
  double y = -(a.unit_sphere_x() * b.unit_sphere_z() - a.unit_sphere_z()
      * b.unit_sphere_x());
  double z = a.unit_sphere_x() * b.unit_sphere_y() - a.unit_sphere_y()
      * b.unit_sphere_x();

  return point(x, y, z, a.weight());
}

// TODO(cbmorrison) Currently great_circle is just a wrapper for cross that
// sets the weight of the returned point to zero.
point point::great_circle(const point& p) {
  point great = cross(p);
  great.set_weight(0);
  return cross(great);
}

point point::great_circle(const point& a, const point& b) {
  point great = cross(a, b);
  great.set_weight(0);
  return great;
}

pixel point::to_pixel() const {
  return pixel(id());
}

pixel point::to_pixel(int level) const {
  return pixel(id(level));
}

void point::rotate_about(const point& axis, double rotation_angle_degrees) {
  // Time to find that python code I wrote so I can make this easy.
}

} // end namespace s2omp
