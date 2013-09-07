#include <gtest/gtest.h>

#include "annulus_bound.h"

#include "angular_bin-inl.h"
#include "circle_bound.h"
#include "pixel.h"
#include "point.h"

TEST(annulus_bound, TestAnnulusBoundDefaultConstructor) {
  // The default constructor without any arguments should produce an invalid
  // annulus_bound
  s2omp::annulus_bound bound;
  ASSERT_FALSE(bound.is_valid());
  ASSERT_TRUE(bound.is_empty());
  ASSERT_EQ(bound.size(), 0);

  // Create an annulus bound centered around the north pole.
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double inner_height = 0.0;
  double outer_height = 0.0;
  bound = s2omp::annulus_bound(axis, inner_height, outer_height);

  ASSERT_DOUBLE_EQ(bound.area(), 0.0);
  ASSERT_DOUBLE_EQ(bound.axis().dot(axis), 1.0);
  ASSERT_DOUBLE_EQ(bound.get_center().dot(axis), 1.0);
  ASSERT_DOUBLE_EQ(bound.inner_height(), inner_height);
  ASSERT_DOUBLE_EQ(bound.outer_height(), outer_height);

  // Create an annulus_bound with the same center but encompassing the northern
  // hemisphere.
  outer_height = 1.0;
  bound = s2omp::annulus_bound(axis, inner_height, outer_height);

  ASSERT_FALSE(bound.is_empty());
  ASSERT_EQ(bound.size(), 1);
  ASSERT_DOUBLE_EQ(bound.area(), 2.0 * s2omp::PI * s2omp::STRAD_TO_DEG2);
  ASSERT_DOUBLE_EQ(bound.inner_radius(), 0.0);
  ASSERT_DOUBLE_EQ(bound.outer_radius(), 90.0);

  // Now increase the inner radius to the same value and verify that we have
  // a non-empty, but zero area bound.
  inner_height = outer_height;
  bound = s2omp::annulus_bound(axis, inner_height, outer_height);

  ASSERT_FALSE(bound.is_empty());
  ASSERT_EQ(bound.size(), 1);
  ASSERT_NEAR(bound.area(), 0.0, 1.0e-6);
  ASSERT_DOUBLE_EQ(bound.outer_radius(), 90.0);
  ASSERT_DOUBLE_EQ(bound.inner_radius(), 90.0);

  // Verify that an annulus_bound with inner_height > outer_height is invalid.
  bound = s2omp::annulus_bound(axis, 1.0, 0.0);
  ASSERT_FALSE(bound.is_valid());
}

TEST(annulus_bound, TestAnnulusBoundAltConstructors) {
  // Test the alternate constructors.

  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double theta_min_deg = 0.01;
  double theta_max_deg = 10.0;
  double inner_height =
      s2omp::circle_bound::get_height_for_angle(theta_min_deg);
  double outer_height =
      s2omp::circle_bound::get_height_for_angle(theta_max_deg);
  double expected_precision = 1.0e-10;
  s2omp::angular_bin bin(theta_min_deg, theta_max_deg);

  s2omp::annulus_bound* bound =
      s2omp::annulus_bound::from_angular_bin(axis, bin);

  // Default checks to make sure that we've got a valid bound.
  ASSERT_FALSE(bound->is_empty());
  ASSERT_EQ(bound->size(), 1);
  ASSERT_NEAR(bound->inner_radius(), theta_min_deg, expected_precision);
  ASSERT_NEAR(bound->outer_radius(), theta_max_deg, expected_precision);

  // The areas for the annulus_bound and angular_bin should be equal
  ASSERT_NEAR(bound->area(), bin.area(), expected_precision);

  // Likewise, the circle_bound should contain the axis, while the angular_bin
  // generally wouldn't.
  ASSERT_FALSE(bound->contains(axis));

  bound->clear();
  delete bound;

  // Construct what should be an equivalent bound from innner and outer heights.
  bound = s2omp::annulus_bound::from_heights(axis, inner_height, outer_height);
  ASSERT_FALSE(bound->is_empty());
  ASSERT_EQ(bound->size(), 1);
  ASSERT_NEAR(bound->inner_radius(), theta_min_deg, expected_precision);
  ASSERT_NEAR(bound->outer_radius(), theta_max_deg, expected_precision);
  ASSERT_NEAR(bound->area(), bin.area(), expected_precision);
  ASSERT_FALSE(bound->contains(axis));
  bound->clear();
  delete bound;

  // And finally, the same from the radii directly.
  bound = s2omp::annulus_bound::from_radii(axis, theta_min_deg, theta_max_deg);
  ASSERT_FALSE(bound->is_empty());
  ASSERT_EQ(bound->size(), 1);
  ASSERT_NEAR(bound->inner_radius(), theta_min_deg, expected_precision);
  ASSERT_NEAR(bound->outer_radius(), theta_max_deg, expected_precision);
  ASSERT_NEAR(bound->area(), bin.area(), expected_precision);
  ASSERT_FALSE(bound->contains(axis));
  bound->clear();
  delete bound;
}

TEST(annulus_bound, TestAnnulusBoundContainment) {
  // Test the containment routines for points and pixels.

  // Start with an empty bound.
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double inner_height = 0.0;
  double outer_height = 0.0;
  s2omp::point anti_axis = -axis;
  s2omp::pixel pix = axis.to_pixel();
  s2omp::annulus_bound bound(axis, inner_height, outer_height);

  ASSERT_FALSE(bound.contains(axis));
  ASSERT_FALSE(bound.contains(anti_axis));
  ASSERT_FALSE(bound.contains(pix));

  // Now replace that bound with one that encloses the entire northern
  // hemisphere.
  outer_height = 1.0;
  bound = s2omp::annulus_bound(axis, inner_height, outer_height);

  ASSERT_FALSE(bound.contains(axis));
  ASSERT_FALSE(bound.contains(anti_axis));
  ASSERT_FALSE(bound.contains(pix));

  // If we instead use a point that's at RA=45,DEC=45, that should be
  // contained, likewise the pixel from that point.
  s2omp::point off_axis = s2omp::point::from_radec_deg(45.0, 45.0);
  s2omp::pixel off_axis_pix = off_axis.to_pixel();
  ASSERT_TRUE(bound.contains(off_axis));
  ASSERT_TRUE(bound.contains(off_axis_pix));
  ASSERT_TRUE(bound.contains(off_axis_pix.parent(20)));

  // Since we contain these pixels, we also expect that the contained
  // area should be equal to the pixel area
  ASSERT_DOUBLE_EQ(bound.contained_area(off_axis_pix),
      off_axis_pix.exact_area());

  // Now shrink the annulus_bound to a case that we'd expect for small
  // angle pair finding.
  double theta_min_deg = 0.001;
  double theta_max_deg = 0.002;
  double theta = 0.5 * (theta_min_deg + theta_max_deg);
  inner_height = s2omp::circle_bound::get_height_for_angle(theta_min_deg);
  outer_height = s2omp::circle_bound::get_height_for_angle(theta_max_deg);
  bound = s2omp::annulus_bound(axis, inner_height, outer_height);
  off_axis = s2omp::point::from_radec_deg(45.0, 90.0 - theta);
  off_axis_pix = off_axis.to_pixel();

  ASSERT_FALSE(bound.contains(axis));
  ASSERT_FALSE(bound.contains(anti_axis));
  ASSERT_TRUE(bound.contains(off_axis));

  // Leaf pixels should still be contained, but a pixel on the scale of the
  // circle bound won't be entirely contained.
  s2omp::pixel small_pix = off_axis_pix.parent(15);
  ASSERT_TRUE(bound.contains(off_axis_pix));
  ASSERT_FALSE(bound.contains(small_pix));
  ASSERT_TRUE(bound.may_intersect(small_pix));

  // Larger ones won't be at all
  s2omp::pixel mid_pix = off_axis_pix.parent(10);
  ASSERT_FALSE(bound.contains(mid_pix));
  ASSERT_TRUE(bound.may_intersect(mid_pix));

  // The contained area for the larger pixels should be the area of the
  // annulus_bound.  For the pixel on the scale of the annulus_bound, we expect
  // to have partial overlap, so the contained area should be less than the
  // annulus_bound's area.
  ASSERT_LT(bound.contained_area(small_pix), bound.area());
  ASSERT_DOUBLE_EQ(bound.contained_area(mid_pix), bound.area());
}

TEST(annulus_bound, TestAnnulusBoundRandomPoints) {
  // Test the routines for generating random points with the annulus_bound

  // Start with an axis at the north pole.
  s2omp::point axis(0.0, 0.0, 1.0, 1.0);
  double theta_min_deg = 0.01;
  double theta_max_deg = 10.0;

  s2omp::annulus_bound* bound =
      s2omp::annulus_bound::from_radii(axis, theta_min_deg, theta_max_deg);

  s2omp::point p = bound->get_random_point();
  ASSERT_TRUE(bound->contains(p));

  s2omp::point_vector points;
  int n_points = 1000;
  bound->get_random_points(n_points, &points);
  for (int k = 0; k < n_points; k++) {
    ASSERT_TRUE(bound->contains(points[k]));
  }

  // Now switch to a circle bound on the equator
  axis = s2omp::point::from_radec_deg(0.0, 0.0);
  s2omp::annulus_bound* eq_bound =
      s2omp::annulus_bound::from_radii(axis, theta_min_deg, theta_max_deg);

  p = eq_bound->get_random_point();
  ASSERT_TRUE(eq_bound->contains(p));
  ASSERT_FALSE(bound->contains(p));

  eq_bound->get_random_points(n_points, &points);
  for (int k = 0; k < n_points; k++) {
    ASSERT_TRUE(eq_bound->contains(points[k]));
    ASSERT_FALSE(bound->contains(points[k]));
  }

  bound->clear();
  eq_bound->clear();
  delete bound;
  delete eq_bound;
}
