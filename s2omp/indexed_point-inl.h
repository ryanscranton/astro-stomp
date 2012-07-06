/*
 * indexed_point.h
 *
 *  Created on: Jul 5, 2012
 *      Author: scranton
 */

#ifndef INDEXED_POINT_H_
#define INDEXED_POINT_H_

// #include <boost/scoped_ptr.hpp>
#include "point.h"

namespace s2omp {

class indexed_point: public point {
public:
	indexed_point(double x, double y, double z, uint32_t index);

	static indexed_point* from_latlon_deg(double lat_deg, double lon_deg,
			Sphere s, uint32_t index);
	static indexed_point* from_latlon_rad(double lat_rad, double lon_rad,
			Sphere s, uint32_t index);

	static indexed_point* from_radec_deg(double ra_deg, double dec_deg,
			uint32_t index);
	static indexed_point* from_radec_rad(double ra_rad, double dec_rad,
			uint32_t index);

	static indexed_point* from_point(const point& p, uint32_t index);

	uint32_t index() const;
	void set_index(uint32_t index);

private:
	uint32_t index_;
};

inline indexed_point::indexed_point(double x, double y, double z,
		uint32_t index) {
	set_xyz(x, y, z);
	set_weight(1.0);
	index_ = index;
}

inline indexed_point* indexed_point::from_latlon_deg(double lat_deg,
		double lon_deg, Sphere s, uint32_t index) {
	return from_latlon_rad(lat_deg * DEG_TO_RAD, lon_deg * DEG_TO_RAD, s, index);
}

inline indexed_point* indexed_point::from_latlon_rad(double lat_rad,
		double lon_rad, Sphere s, uint32_t index) {
	indexed_point* p = new indexed_point(1.0, 0.0, 0.0, index);
	p->set_latlon_rad(lat_rad, lon_deg, s);
	return p;
}

inline indexed_point* indexed_point::from_radec_deg(double ra_deg,
		double dec_deg, Sphere s, uint32_t index) {
	return from_latlon_rad(dec_deg * DEG_TO_RAD, ra_deg * DEG_TO_RAD, EQUATORIAL,
			index);
}

inline indexed_point* indexed_point::from_radec_rad(double ra_rad,
		double dec_rad, Sphere s, uint32_t index) {
	indexed_point* p = new indexed_point(1.0, 0.0, 0.0, index);
	p->set_latlon_rad(dec_rad, ra_deg, EQUATORIAL);
	return p;
}

inline indexed_point* indexed_point::from_point(const point& p, uint32_t index) {
	return new indexed_point(point.unit_sphere_x(), point.unit_sphere_y(),
			point.unit_sphere_z(), index);
}

} // end namespace s2omp

#endif /* INDEXED_POINT_H_ */
