#include <gtest/gtest.h>

#include "circle_bound.h"

#include "angular_bin-inl.h"
#include "pixel.h"
#include "point.h"

TEST(circle_bound, TestCircleBoundDefaultConstructor) {
  // The default constructor without any arguments should produce an invalid
  // circle_bound
  s2omp::circle_bound bound;
  ASSERT_FALSE(bound.is_valid());
  ASSERT_TRUE(bound.is_empty());
  ASSERT_EQ(bound.size(), 0);

  // Create a circle_bound centered on the north Equatorial pole with 0 radius.
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double height = 0.0;
  bound = s2omp::circle_bound(axis, height);

  ASSERT_DOUBLE_EQ(bound.area(), 0.0);
  ASSERT_DOUBLE_EQ(bound.axis().dot(axis), 1.0);
  ASSERT_DOUBLE_EQ(bound.get_center().dot(axis), 1.0);
  ASSERT_DOUBLE_EQ(bound.height(), height);
  ASSERT_DOUBLE_EQ(bound.radius(), 0.0);

  // Create a circle_bound with the same center but encompassing the northern
  // hemisphere.
  height = 1.0;
  bound = s2omp::circle_bound(axis, height);

  ASSERT_FALSE(bound.is_empty());
  ASSERT_EQ(bound.size(), 1);
  ASSERT_DOUBLE_EQ(bound.area(), 2.0 * s2omp::PI * s2omp::STRAD_TO_DEG2);
  ASSERT_DOUBLE_EQ(bound.radius(), 90.0);
}

TEST(circle_bound, TestCircleBoundAltConstructors) {
  // Test the alternate constructors.
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double theta_min_deg = 0.01;
  double theta_max_deg = 10.0;
  s2omp::angular_bin bin(theta_min_deg, theta_max_deg);

  s2omp::circle_bound* bound = s2omp::circle_bound::from_angular_bin(axis, bin);

  // Default checks to make sure that we've got a valid bound.
  ASSERT_FALSE(bound->is_empty());
  ASSERT_EQ(bound->size(), 1);
  ASSERT_DOUBLE_EQ(bound->radius(), theta_max_deg);

  // Since the angular bin has an inner radius, we expect that the area of
  // our circle_bound is larger than the angular_bin.
  ASSERT_GT(bound->area(), bin.area());

  // Likewise, the circle_bound should contain the axis, while the angular_bin
  // generally wouldn't.
  ASSERT_TRUE(bound->contains(axis));

  delete bound;
}

TEST(circle_bound, TestCircleBoundContainment) {
  // Test the containment routines for points and pixels.

  // Start with an empty bound.
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double height = 0.0;
  s2omp::point anti_axis = -axis;
  s2omp::pixel pix(axis.id());
  s2omp::circle_bound bound(axis, height);

  ASSERT_TRUE(bound.contains(axis));
  ASSERT_FALSE(bound.contains(anti_axis));
  ASSERT_FALSE(bound.contains(pix));

  // Now replace that bound with one that encloses the entire northern
  // hemisphere.
  height = 1.0;
  bound = s2omp::circle_bound(axis, height);

  ASSERT_TRUE(bound.contains(axis));
  ASSERT_FALSE(bound.contains(anti_axis));
  ASSERT_TRUE(bound.contains(pix));
  // Since our axis is at the north pole, this should also include the
  // corresponding face pixel.
  s2omp::pixel face = pix.parent(0);
  ASSERT_TRUE(bound.contains(face));
  // Since we contain these pixels, we also expect that the contained
  // area should be equal to the pixel area
  ASSERT_DOUBLE_EQ(bound.contained_area(pix), pix.exact_area());
  ASSERT_DOUBLE_EQ(bound.contained_area(face), face.exact_area());

  // Now shrink the circle_bound to a case that we'd expect for small
  // angle pair finding.
  double theta_deg = 0.001;
  height = s2omp::circle_bound::get_height_for_angle(theta_deg);
  bound = s2omp::circle_bound(axis, height);

  ASSERT_TRUE(bound.contains(axis));
  ASSERT_FALSE(bound.contains(anti_axis));

  // Leaf pixels should still be contained, but a pixel on the scale of the
  // circle bound won't be entirely contained.
  s2omp::pixel small_pix = pix.parent(15);
  ASSERT_TRUE(bound.contains(pix));
  ASSERT_FALSE(bound.contains(small_pix));
  ASSERT_TRUE(bound.may_intersect(small_pix));

  // Larger ones won't be at all
  s2omp::pixel mid_pix = pix.parent(10);
  ASSERT_FALSE(bound.contains(mid_pix));
  ASSERT_FALSE(bound.contains(face));
  ASSERT_TRUE(bound.may_intersect(mid_pix));
  ASSERT_TRUE(bound.may_intersect(face));

  // The contained area for the larger pixels should be the area of the
  // circle_bound.  For the pixel on the scale of the circle_bound, we expect
  // to have partial overlap, so the contained area should be less than the
  // circle_bound's area.
  ASSERT_LT(bound.contained_area(small_pix), bound.area());
  ASSERT_DOUBLE_EQ(bound.contained_area(mid_pix), bound.area());
  ASSERT_DOUBLE_EQ(bound.contained_area(face), bound.area());
}

TEST(circle_bound, TestCircleBoundRandomPoints) {
  // Test the routines for generating random points with the circle_bound

  // Start with an axis at the north pole.
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double theta_deg = 10.0;
  double height = s2omp::circle_bound::get_height_for_angle(theta_deg);
  s2omp::circle_bound bound(axis, height);

  s2omp::point p = bound.create_random_point();
  ASSERT_TRUE(bound.contains(p));

  s2omp::point_vector points;
  int n_points = 1000;
  bound.get_random_points(n_points, &points);
  for (int k = 0; k < n_points; k++) {
    ASSERT_TRUE(bound.contains(points[k]));
  }

  // Now switch to a circle bound on the equator
  axis = s2omp::point::from_radec_deg(0.0, 0.0);
  s2omp::circle_bound eq_bound(axis, height);

  p = eq_bound.create_random_point();
  ASSERT_TRUE(eq_bound.contains(p));
  ASSERT_FALSE(bound.contains(p));

  eq_bound.get_random_points(n_points, &points);
  for (int k = 0; k < n_points; k++) {
    ASSERT_TRUE(eq_bound.contains(points[k]));
    ASSERT_FALSE(bound.contains(points[k]));
  }
}

TEST(circle_bound, TestCircleBoundAddPoint) {
  // Test the routines that potentially alter the bound by making sure that
  // it includes input points.

  // Start with an axis at the north pole.
  s2omp::point axis = s2omp::point::from_radec_deg(0.0, 90.0);
  double theta_deg = 10.0;
  double height = s2omp::circle_bound::get_height_for_angle(theta_deg);
  s2omp::circle_bound bound(axis, height);
  ASSERT_NEAR(bound.axis().dot(axis), 1.0, 1.0e-10);
  ASSERT_NEAR(bound.radius(), theta_deg, 1.0e-10);

  // Now add a point that's within the current bound and verify that the bound
  // remains unchanged.
  s2omp::point p = s2omp::point::from_radec_deg(axis.ra_deg(),
          axis.dec_deg() - 0.5 * theta_deg);
  bound.add_point(p);
  ASSERT_TRUE(bound.contains(p));
  ASSERT_NEAR(bound.axis().dot(axis), 1.0, 1.0e-10);
  ASSERT_NEAR(bound.radius(), theta_deg, 1.0e-10);

  // Add a point that's outside the bound.  The bound should expand to include
  // the point.
  p = s2omp::point::from_radec_deg(axis.ra_deg(),
        axis.dec_deg() - 2.0 * theta_deg);
  bound.add_point(p);
  ASSERT_TRUE(bound.contains(p));
  ASSERT_NEAR(bound.axis().dot(axis), 1.0, 1.0e-10);
  ASSERT_NEAR(bound.radius(), 2.0 * theta_deg, 1.0e-10);
}

TEST(circle_bound, TestCircleBoundAddBound) {
  // Test the method for adding a circle_bound to a current circle_bound.

  // Start with an axis at the north pole.
  s2omp::point axis = s2omp::point::from_radec_deg(0.0, 90.0);
  double theta_deg = 10.0;
  double height = s2omp::circle_bound::get_height_for_angle(theta_deg);
  s2omp::circle_bound bound(axis, height);
  ASSERT_NEAR(bound.radius(), theta_deg, 1.0e-10);

  // Start by adding a circle_bound that's inside our current bound.  Verify
  // that this leaves our current bound unchanged.
  bound.add_circle_bound(s2omp::circle_bound(axis, 0.5 * height));
  ASSERT_NEAR(bound.axis().dot(axis), 1.0, 1.0e-10);
  ASSERT_NEAR(bound.radius(), theta_deg, 1.0e-10);

  // Now add a circle_bound that's outside our original bound.
  s2omp::point p = s2omp::point::from_radec_deg(axis.ra_deg(),
          axis.dec_deg() - 3.0 * theta_deg);
  ASSERT_FALSE(bound.contains(p));

  s2omp::circle_bound outside_bound(p, height);
  // Random points generated inside this bound should not be contained by the
  // original circle_bound.
  int n_random = 1000;
  s2omp::point_vector points;
  outside_bound.get_random_points(n_random, &points);
  for (int k = 0; k < n_random; k++) {
    ASSERT_FALSE(bound.contains(points[k]));
  }

  bound.add_circle_bound(outside_bound);
  // The resulting circle_bound should have the same axis, but have a radius
  // 4 times as large.
  ASSERT_NEAR(bound.axis().dot(axis), 1.0, 1.0e-10);
  ASSERT_NEAR(bound.radius(), 4.0 * theta_deg, 1.0e-10);
  ASSERT_TRUE(bound.contains(p));

  // The random points we generated previously should now all be contained
  // within the modified circle_bound.
  for (int k = 0; k < n_random; k++) {
    ASSERT_TRUE(bound.contains(points[k]));
  }
}
