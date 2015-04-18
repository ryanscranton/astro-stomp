#include <gtest/gtest.h>

#include "polygon_bound.h"

#include "angular_bin-inl.h"
#include "circle_bound.h"
#include "pixel.h"
#include "point.h"

void get_default_loop(s2omp::point_vector& points) {
  points.clear();
  points.reserve(4);
  points.push_back(s2omp::point::from_radec_deg(-1.0, 1.0));
  points.push_back(s2omp::point::from_radec_deg(-1.0, -1.0));
  points.push_back(s2omp::point::from_radec_deg(1.0, -1.0));
  points.push_back(s2omp::point::from_radec_deg(1.0, 1.0));
}

void get_intersecting_loop(s2omp::point_vector& points) {
  points.clear();
  points.reserve(4);
  points.push_back(s2omp::point::from_radec_deg(0.0, 0.0));
  points.push_back(s2omp::point::from_radec_deg(0.0, -2.0));
  points.push_back(s2omp::point::from_radec_deg(2.0, -2.0));
  points.push_back(s2omp::point::from_radec_deg(2.0, 0.0));
}

void get_disjoint_loop(s2omp::point_vector& points) {
  points.clear();
  points.reserve(4);
  points.push_back(s2omp::point::from_radec_deg(9.0, 1.0));
  points.push_back(s2omp::point::from_radec_deg(9.0, -1.0));
  points.push_back(s2omp::point::from_radec_deg(11.0, -1.0));
  points.push_back(s2omp::point::from_radec_deg(11.0, 1.0));
}

s2omp::polygon_bound* get_default_bound() {
  s2omp::point_vector points;
  get_default_loop(points);
  return s2omp::polygon_bound::from_points(points);
}

s2omp::polygon_bound* get_intersecting_bound() {
  s2omp::point_vector points;
  get_intersecting_loop(points);
  return s2omp::polygon_bound::from_points(points);
}

s2omp::polygon_bound* get_disjoint_bound() {
  s2omp::point_vector points;
  get_disjoint_loop(points);
  return s2omp::polygon_bound::from_points(points);
}

TEST(polygon_bound, TestPolygonBoundDefaultConstructor) {
  // The default constructor without any points should produce an invalid bound.
  s2omp::point_vector points;
  s2omp::polygon_bound* bound = s2omp::polygon_bound::from_points(points);
  ASSERT_FALSE(bound->is_valid());
  ASSERT_TRUE(bound->is_empty());
  ASSERT_EQ(bound->size(), 0);
  delete bound;

  // Create a polygon_bound from the default loop.
  bound = get_default_bound();
  ASSERT_EQ(bound->num_loops(), 1);
  ASSERT_EQ(bound->size(), 4);
  ASSERT_EQ(bound->size(), bound->num_vertices());
  ASSERT_FALSE(bound->is_empty());
  ASSERT_TRUE(bound->is_valid());
  delete bound;

  // Verify that our other two default polygons are valid.
  bound = get_intersecting_bound();
  ASSERT_TRUE(bound->is_valid());
  delete bound;

  bound = get_disjoint_bound();
  ASSERT_TRUE(bound->is_valid());
}

TEST(polygon_bound, TestAddDisjointLoop) {
  // Start with the default polygon and add a disjoint loop to it.
  s2omp::polygon_bound* bound = get_default_bound();
  double area = bound->area();

  s2omp::point_vector points;
  get_disjoint_loop(points);
  ASSERT_TRUE(bound->add_loop(points));
  ASSERT_GT(bound->area(), area);
  ASSERT_EQ(bound->num_loops(), 2);
  ASSERT_EQ(bound->num_vertices(), 8);
}

TEST(polygon_bound, TestAddInteriorLoop) {
  // Start with the default polygon and add an interior loop to it.
  s2omp::polygon_bound* bound = get_default_bound();
  double starting_area = bound->area();

  s2omp::point_vector points;
  points.push_back(s2omp::point::from_radec_deg(-0.1, 0.1));
  points.push_back(s2omp::point::from_radec_deg(-0.1, -0.1));
  points.push_back(s2omp::point::from_radec_deg(0.1, -0.1));
  points.push_back(s2omp::point::from_radec_deg(0.1, 0.1));

  ASSERT_TRUE(bound->add_loop(points));
  ASSERT_LT(bound->area(), starting_area);
  ASSERT_EQ(bound->num_loops(), 2);
  ASSERT_EQ(bound->num_vertices(), 8);
}

TEST(polygon_bound, TestIntersects) {
  // Start with the default polygon.
  s2omp::polygon_bound* bound = get_default_bound();

  // This should intersect with the default intersecting polygon.
  s2omp::polygon_bound* bound2 = get_intersecting_bound();
  ASSERT_TRUE(bound->intersects(*bound2));

  // The original bound should not intersect the default disjoint polygon.
  s2omp::polygon_bound* bound3 = get_disjoint_bound();
  ASSERT_FALSE(bound->intersects(*bound3));
}

TEST(polygon_bound, TestContainsPoint) {
  s2omp::polygon_bound* bound = get_default_bound();

  // Should contain the origin.
  ASSERT_TRUE(bound->contains(s2omp::point::from_radec_deg(0.0, 0.0)));

  // Shouldn't contain a point just outside.
  ASSERT_FALSE(bound->contains(s2omp::point::from_radec_deg(1.0, 1.00001)));

  // This particular bound should contain its center.
  ASSERT_TRUE(bound->contains(bound->get_center()));
}

TEST(polygon_bound, TestContainsPixel) {
  // Our default bound should contain a max resolution pixel at the origin.
  s2omp::polygon_bound* bound = get_default_bound();
    s2omp::pixel pix =
        s2omp::pixel::from_point(s2omp::point::from_radec_deg(0.0, 0.0));
  ASSERT_TRUE(bound->contains(pix));
  ASSERT_TRUE(bound->may_intersect(pix));

  // If we consider the much larger level 0 parent, this should not be the case.
  ASSERT_FALSE(bound->contains(pix.parent(0)));
  ASSERT_TRUE(bound->may_intersect(pix.parent(0)));

  // Likewise, if we consider a near leaf cell on one of the vertices, this
  // should not be contained but may intersect.
  pix = s2omp::pixel::from_point(s2omp::point::from_radec_deg(-1.0, -1.0), 25);
  ASSERT_FALSE(bound->contains(pix));
  ASSERT_TRUE(bound->may_intersect(pix));
}

TEST(polygon_bound, TestContainedAreaInPixel) {
  // A fully contained pixel should have a contained area equal to its exact
  // area.
  s2omp::polygon_bound* bound = get_default_bound();
  s2omp::pixel pix =
      s2omp::pixel::from_point(s2omp::point::from_radec_deg(0.0, 0.0));
  ASSERT_DOUBLE_EQ(pix.exact_area(), bound->contained_area(pix));

  // Likewise, if the pixel is disjoint, then we expect the contained area to
  // be zero.
  pix = s2omp::pixel::from_point(s2omp::point::from_radec_deg(10.0, 0.0), 25);
  ASSERT_DOUBLE_EQ(bound->contained_area(pix), 0.0);

  // If we have partial overlap, then the precise intersection area is
  // calculated.
  s2omp::point_vector vertices;
  for (int k = 0; k < 4; k++) {
    vertices.push_back(pix.vertex(k));
  }
  delete bound;
  bound = s2omp::polygon_bound::from_points(vertices);
  ASSERT_DOUBLE_EQ(bound->area(), bound->contained_area(pix.parent(20)));
}

TEST(polygon_bound, TestGetBound) {
  s2omp::polygon_bound* bound = get_default_bound();
  s2omp::circle_bound cap = bound->get_bound();
  ASSERT_TRUE(bound->contains(cap.axis()));
}
