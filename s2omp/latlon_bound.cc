#include <s2/s2latlngrect.h>

#include "latlon_bound.h"

#include "circle_bound.h"
#include "coverer.h"
#include "pixel.h"
#include "point.h"

namespace s2omp {

latlon_bound::latlon_bound() {
  latlon_ = S2LatLngRect();
}

latlon_bound::latlon_bound(const point& hi, const point& lo) {
  latlon_ = S2LatLngRect(S2LatLng(hi.s2point()), S2LatLng(lo.s2point()));
}

latlon_bound::latlon_bound(S2LatLngRect bound) {
  latlon_ = bound;
}
latlon_bound::~latlon_bound() {
}

double degrees_to_radians(double angle_degrees) {
  double angle_rad = angle_degrees * DEG_TO_RAD;
  if (angle_rad > PI)
    angle_rad -= 2.0 * PI;
  if (angle_rad < -PI)
    angle_rad += 2.0 * PI;

  return angle_rad;
}

latlon_bound latlon_bound::from_ra_bounds(double ra_min_degrees,
    double ra_max_degrees) {
  S2LatLngRect bound = S2LatLngRect(S2LatLngRect::FullLat(), S1Interval(
      degrees_to_radians(ra_min_degrees), degrees_to_radians(ra_max_degrees)));

  return latlon_bound(bound);
}

latlon_bound latlon_bound::from_dec_bounds(double dec_min_degrees,
    double dec_max_degrees) {
  S2LatLngRect bound = S2LatLngRect(R1Interval(dec_min_degrees * DEG_TO_RAD,
      dec_max_degrees * DEG_TO_RAD), S2LatLngRect::FullLng());

  return latlon_bound(bound);
}

void latlon_bound::add_point(const point& p) {
  latlon_.AddPoint(p.s2point());
}

point latlon_bound::vertex(int k) const {
  return point(latlon_.GetVertex(k), 1.0);
}

bool latlon_bound::is_valid() const {
  return latlon_.is_valid();
}

bool latlon_bound::is_empty() const {
  return latlon_.is_empty();
}

long latlon_bound::size() const {
  return is_empty() ? 0 : 1;
}

void latlon_bound::clear() {
  latlon_ = S2LatLngRect();
}

double latlon_bound::area() const {
  return latlon_.Area() * STRAD_TO_DEG2;
}

bool latlon_bound::contains(const point& p) const {
  return latlon_.Contains(p.s2point());
}

bool latlon_bound::contains(const pixel& pix) const {
  return latlon_.Contains(pix.get_cell());
}

double latlon_bound::contained_area(const pixel& pix) const {
  if (contains(pix)) {
    return pix.exact_area();
  }

  if (!may_intersect(pix)) {
    return 0.0;
  }

  // TODO(scranton): This is something that could probably be done analytically,
  // but for now, we'll do something simpler.
  double total_area = 0.0;
  int sampling_level = min(pix.level() + 5, MAX_LEVEL);
  pixel child_end = pix.child_end(sampling_level);
  for (pixel child_pix = pix.child_begin(sampling_level); child_pix
      != child_end; child_pix.next()) {
    total_area += contains(child_pix) ? child_pix.exact_area() : 0.0;
  }

  return total_area;
}

bool latlon_bound::may_intersect(const pixel& pix) const {
  return latlon_.MayIntersect(pix.get_cell());
}

point latlon_bound::get_center() const {
  return point(latlon_.GetCenter(), 1.0);
}

circle_bound latlon_bound::get_bound() const {
  circle_bound bound = circle_bound(get_center(), 0.0);
  for (int k = 0; k < 4; k++) {
    bound.add_point(vertex(k));
  }
  return bound;
}

} // end namespace s2omp
