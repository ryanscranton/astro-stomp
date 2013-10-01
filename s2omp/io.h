// Copyright 2013  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains a variety of static methods for reading and writing
// various s2omp objects to and from files.  We intend s2omp to be used with a
// variety of formats and interfaces and therefore want to keep the library
// objects themselves agnostic about how and where the source data come from.
// Hence, we create a separate class to act as a source and sink for all of
// those I/O operations.

#ifndef S2OMP_IO_H_
#define S2OMP_IO_H_

#include <map>
#include <queue>
#include <string>

#include "core.h"

#include "s2omp.pb.h"

namespace s2omp {

class field_pixel;
class pixel_union;
class field_union;
class tree_union;

class io {
public:
  static bool write_ascii(const point_vector& points,
      const string& file_name);
  static bool write_ascii(const pixel_vector& pixels,
      const string& file_name);
  static bool write_ascii(const pixel_union& pix_union,
      const string& file_name);

  static bool write_pb(const point_vector& points,
      const string& output_file);
  static bool write_pb(const pixel_vector& pixels,
      const string& output_file);
  static bool write_pb(const pixel_union& pix_union,
      const string output_file);
  static bool write_pb(const field_union& s, const string& output_file);
  static bool write_pb(const tree_union& t, const string& output_file);

  static bool read_ascii(const string& input_file, point_vector* points);
  static bool read_ascii(const string& input_file, pixel_vector* pixels);
  static bool read_ascii(const string& input_file, pixel_union* pix_union);

  static bool read_pb(const string& input_file, point_vector* points);
  static bool read_pb(const string& input_file, pixel_vector* pixels);
  static bool read_pb(const string& input_file, pixel_union* pix_union);
  static bool read_pb(const string& input_file, field_union* s);
  static bool read_pb(const string& input_file, tree_union* t);

private:
  static point to_point(const PointProto& proto);
  static pixel to_pixel(const PixelProto& proto);
  static field_pixel to_field_pixel(const PixelProto& proto);
};

} // end namespace s2omp

#endif /* IO_H_ */
