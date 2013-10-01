#include <gtest/gtest.h>

#include "io.h"

#include "circle_bound.h"
#include "pixel_union.h"

TEST(io, TestIoPixelAscii) {
  // Test the static methods for reading and writing pixels to ascii files.

  // Start with circle_bound
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Choose a starting resolution that's close to the area of the bound and
  // generate a simple covering.
  int level = s2omp::pixel::get_level_from_area(bound->area() / 10.0);
  s2omp::pixel_vector pixels;
  bound->get_simple_covering(level, &pixels);

  // Write the covering to an ascii file
  string test_file_name = "io_test_pixel_ascii.txt";
  ASSERT_TRUE(s2omp::io::write_ascii(pixels, test_file_name));

  // Re-read the file and verify that the contents are identical.
  s2omp::pixel_vector read_pixels;
  ASSERT_TRUE(s2omp::io::read_ascii(test_file_name, &read_pixels));

  ASSERT_EQ(pixels.size(), read_pixels.size());
  for (int k = 0; k < pixels.size(); k++) {
    ASSERT_EQ(pixels[k].id(), read_pixels[k].id());
  }
}

TEST(io, TestIoPixelUnionAscii) {
  // Test the static methods for reading and writing pixel_unions to
  // ascii files.

  // Start with circle_bound
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Choose a starting resolution that's close to the area of the bound and
  // generate a simple covering to serve as the source of a pixel_union.
  int level = s2omp::pixel::get_level_from_area(bound->area() / 10.0);
  s2omp::pixel_vector pixels;
  bound->get_simple_covering(level, &pixels);
  s2omp::pixel_union pix_union;
  pix_union.init(pixels);

  // Write the pixel_union to an ascii file
  string test_file_name = "io_test_pixel_union_ascii.txt";
  ASSERT_TRUE(s2omp::io::write_ascii(pix_union, test_file_name));

  // Re-read the file and verify that the contents are identical.
  s2omp::pixel_union read_union;
  ASSERT_TRUE(s2omp::io::read_ascii(test_file_name, &read_union));

  ASSERT_EQ(pix_union.size(), read_union.size());
  ASSERT_EQ(pix_union.min_level(), read_union.min_level());
  ASSERT_EQ(pix_union.max_level(), read_union.max_level());
  ASSERT_EQ(pix_union.range_min().id(), read_union.range_min().id());
  ASSERT_EQ(pix_union.range_max().id(), read_union.range_max().id());
  ASSERT_NEAR(pix_union.area(), read_union.area(), 1.0e-10);
  for (s2omp::pixel_iterator iter = pix_union.begin();
      iter != pix_union.end(); ++iter) {
    ASSERT_TRUE(read_union.contains(*iter));
  }
}

TEST(io, TestIoPointAscii) {
  // Test the static methods for reading and writing points to ascii files.

  // Start with circle_bound
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Choose a starting resolution that's close to the area of the bound and
  // generate some random points within the bound.
  s2omp::point_vector points;
  long n_points = 100;
  bound->get_random_points(n_points, &points);

  // Write the points to an ascii file
  string test_file_name = "io_test_point_ascii.txt";
  ASSERT_TRUE(s2omp::io::write_ascii(points, test_file_name));

  // Re-read the file and verify that the contents are identical.
  s2omp::point_vector read_points;
  ASSERT_TRUE(s2omp::io::read_ascii(test_file_name, &read_points));

  ASSERT_EQ(points.size(), read_points.size());
  for (int k = 0; k < points.size(); k++) {
    ASSERT_EQ(points[k].id(), read_points[k].id());
    ASSERT_NEAR(points[k].weight(), read_points[k].weight(), 1.0e-20);
  }
}

TEST(io, TestIoPixelProtoBuffer) {
  // Test the static methods for reading and writing pixels to proto buffers.

  // Start with circle_bound
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Choose a starting resolution that's close to the area of the bound and
  // generate a simple covering.
  int level = s2omp::pixel::get_level_from_area(bound->area() / 10.0);
  s2omp::pixel_vector pixels;
  bound->get_simple_covering(level, &pixels);

  // Write the covering to a proto buffer
  string test_file_name = "io_test_pixel_vector.pb";
  ASSERT_TRUE(s2omp::io::write_pb(pixels, test_file_name));

  // Re-read the file and verify that the contents are identical.
  s2omp::pixel_vector read_pixels;
  ASSERT_TRUE(s2omp::io::read_pb(test_file_name, &read_pixels));

  ASSERT_EQ(pixels.size(), read_pixels.size());
  for (int k = 0; k < pixels.size(); k++) {
    ASSERT_EQ(pixels[k].id(), read_pixels[k].id());
  }
}
