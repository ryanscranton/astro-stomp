#include <gtest/gtest.h>
#include <set>

#include "tree_pixel.h"

#include "angular_bin-inl.h"
#include "annulus_bound.h"
#include "util.h"

class TreePixelTest : public testing::Test {
 protected:
  virtual void SetUp() {
    s2omp::point p(0.0, 0.0, 1.0, 1.0);
    level_ = 20;
    max_points_ = s2omp::DEFAULT_NODE_CAPACITY;
    pix_ = s2omp::tree_pixel::from_point(p, level_, max_points_);
  }

  virtual void TearDown() {
    delete pix_;
  }

  s2omp::tree_pixel* pix_;
  int level_;
  int max_points_;
};

TEST(tree_pixel, TestTreePixelConstructors) {
  // The default constructor is private, so we need a point or pixel to
  // instantiate a tree_pixel.
  s2omp::point p(0.0, 0.0, 1.0, 1.0);
  int level = 20;
  s2omp::tree_pixel pix(p.id(level));
  // Basic pixel construction.
  ASSERT_TRUE(pix.is_valid());
  ASSERT_EQ(pix.level(), level);
  ASSERT_EQ(pix.id(), p.id(level));
  ASSERT_TRUE(pix.contains(p));

  // Default point and node counts for an empty pixel.
  ASSERT_EQ(pix.node_capacity(), s2omp::DEFAULT_NODE_CAPACITY);
  ASSERT_EQ(pix.n_nodes(), 1);
  ASSERT_EQ(pix.n_points(), 0);
  ASSERT_DOUBLE_EQ(pix.weight(), 0.0);
  ASSERT_TRUE(pix.is_empty());
  ASSERT_EQ(pix.size(), 0);

  // Verify that our other constructors also yield valid tree_pixels.
  int max_points = 2 * s2omp::DEFAULT_NODE_CAPACITY;
  s2omp::tree_pixel* pix_ptr =
      s2omp::tree_pixel::from_pixel(p.to_pixel(level), max_points);
  ASSERT_TRUE(pix_ptr->is_valid());
  ASSERT_EQ(pix_ptr->node_capacity(), max_points);
  delete pix_ptr;

  pix_ptr = s2omp::tree_pixel::from_point(p, level, max_points);
  ASSERT_TRUE(pix_ptr->is_valid());
  ASSERT_EQ(pix_ptr->node_capacity(), max_points);
  delete pix_ptr;
}

TEST_F(TreePixelTest, TestTreePixelAddPoint) {
  // Test the method for adding points to a tree_pixel.

  // Now try adding a single point to the reference tree_pixel and verify that
  // the tree_pixel has the expected number of points/nodes/weight, etc.
  double expected_precision = 1.0e-10;
  s2omp::point p = pix_->get_random_point();
  ASSERT_TRUE(pix_->add_point(p));
  ASSERT_TRUE(pix_->has_points());
  ASSERT_EQ(pix_->n_points(), 1);
  ASSERT_FALSE(pix_->has_nodes());
  ASSERT_EQ(pix_->n_nodes(), 1);
  ASSERT_NEAR(pix_->weight(), p.weight(), expected_precision);
  ASSERT_FALSE(pix_->is_empty());
  ASSERT_EQ(pix_->size(), 1);

  // Now add max_points points to the pixel.  This should trigger the creation
  // of the 4 immediate subnodes.
  s2omp::point_vector points;
  pix_->get_random_points(pix_->node_capacity(), &points);
  for (s2omp::point_iterator iter = points.begin();
      iter != points.end(); ++iter) {
    ASSERT_TRUE(pix_->add_point(*iter));
  }
  ASSERT_EQ(pix_->n_points(), max_points_ + 1);
  ASSERT_EQ(pix_->n_nodes(), 5);

  // Clear the pixel and verify that we've reverted to the original pixel state.
  pix_->clear();
  ASSERT_TRUE(pix_->is_valid());
  ASSERT_EQ(pix_->level(), level_);
  ASSERT_TRUE(pix_->contains(p));
  ASSERT_EQ(pix_->node_capacity(), max_points_);
  ASSERT_EQ(pix_->n_nodes(), 1);
  ASSERT_EQ(pix_->n_points(), 0);
  ASSERT_DOUBLE_EQ(pix_->weight(), 0.0);
  ASSERT_TRUE(pix_->is_empty());
  ASSERT_EQ(pix_->size(), 0);

  // Now, repeat the same task using a child pixel of the original tree_pixel
  // to generate the random points.  This should spawn a new set of 4 subnodes
  // for each of the parent pixels leading from the child pixel to the
  // reference tree_pixel.
  int child_level = pix_->level() + 3;
  s2omp::pixel child = pix_->child_begin(child_level);
  child.get_random_points(max_points_ + 1, &points);
  for (s2omp::point_iterator iter = points.begin();
      iter != points.end(); ++iter) {
    ASSERT_TRUE(pix_->add_point(*iter));
  }
  ASSERT_EQ(pix_->n_points(), points.size());
  ASSERT_EQ(pix_->n_nodes(), 4 * (child_level - pix_->level() + 1) + 1);
}

TEST_F(TreePixelTest, TestTreePixelCoveringFraction) {
  // Test the methods for calculating covering fraction

  // Verify that an empty pixel has covering_fraction of 0.  Likewise, for
  // parents and children of the reference pixel and pixels that do not overlap.
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(), 0.0);
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(*pix_), 0.0);
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(pix_->parent()), 0.0);
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(pix_->child_begin()), 0.0);
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(pix_->next()), 0.0);

  // Now add a single point.  The covering fraction for the pixel and any
  // children that contain the point should be unity.  For the parent, it should
  // be the ratio of the pixel areas and the non-overlapping pixel should still
  // give 0.
  double expected_precision = 1.0e-10;
  ASSERT_TRUE(pix_->add_point(pix_->get_random_point()));
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(), 1.0);
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(*pix_), 1.0);
  ASSERT_NEAR(pix_->covering_fraction(pix_->parent()),
      pix_->exact_area() / pix_->parent().exact_area(), expected_precision);
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(pix_->next()), 0.0);

  // If we now add enough random points over the whole pixel to trigger
  // creation of subnodes, then results should be unchanged.
  s2omp::point_vector points;
  pix_->get_random_points(pix_->node_capacity(), &points);
  for (int k = 0; k < points.size(); k++) {
    ASSERT_TRUE(pix_->add_point(points[k]));
  }
  ASSERT_NEAR(pix_->covering_fraction(), 1.0, expected_precision);
  ASSERT_NEAR(pix_->covering_fraction(*pix_), 1.0, expected_precision);
  ASSERT_NEAR(pix_->covering_fraction(pix_->parent()),
      pix_->exact_area() / pix_->parent().exact_area(), expected_precision);
  ASSERT_NEAR(pix_->covering_fraction(pix_->child_begin()), 1.0,
      expected_precision);
  ASSERT_DOUBLE_EQ(pix_->covering_fraction(pix_->next()), 0.0);

  // Finally, if we clear the pixel and then add enough points drawn from one
  // of the children to trigger sub-nodes, then we should get a covering
  // fraction that's proportional to the ratio of the child's area to that of
  // the reference pixel.
  pix_->clear();
  s2omp::pixel child = pix_->child_begin();
  child.get_random_points(pix_->node_capacity() + 1, &points);
  for (int k = 0; k < points.size(); k++) {
    ASSERT_TRUE(pix_->add_point(points[k]));
  }
  // Note that the ratio is determined by the level of subnodes created, not
  // the pixel used to generate the random points
  ASSERT_NEAR(pix_->covering_fraction(),
      pix_->child_begin().exact_area() / pix_->exact_area(),
      expected_precision);
  ASSERT_NEAR(pix_->covering_fraction(child), 1.0, expected_precision);
  ASSERT_NEAR(pix_->covering_fraction(pix_->child_begin()), 1.0,
      expected_precision);
  ASSERT_NEAR(pix_->covering_fraction(pix_->child_begin().next()), 0.0,
      expected_precision);
}

TEST_F(TreePixelTest, TestTreePixelCopyPoints) {
  // Test the methods for extracting a copy of the points added to a tree_pixel.

  s2omp::point_vector points;

  // First, verify that extracting the points from an empty tree_pixel yields
  // no points, even if the input point_vector contained points.
  pix_->copy_points(&points);
  ASSERT_EQ(points.size(), 0);
  pix_->copy_points(pix_->parent(), &points);
  ASSERT_EQ(points.size(), 0);
  pix_->copy_points(pix_->child_begin(), &points);
  ASSERT_EQ(points.size(), 0);

  points.push_back(pix_->get_random_point());
  pix_->copy_points(&points);
  ASSERT_EQ(points.size(), 0);

  // Now add a number of points to the tree_pixel and verify that the right
  // number of points come out and that they all are contained in the
  // tree_pixel.
  int n_points = 10 * pix_->node_capacity();
  s2omp::point_vector random_points;
  pix_->get_random_points(n_points, &random_points);
  for (int k = 0; k < random_points.size(); k++) {
    ASSERT_TRUE(pix_->add_point(random_points[k]));
  }
  pix_->copy_points(&points);
  ASSERT_EQ(points.size(), random_points.size());
  for (int k = 0; k < points.size(); k++) {
    ASSERT_TRUE(pix_->contains(points[k]));
  }

  // Now verify that the points extracted for parent and child pixels of the
  // reference pixel behave as expected.
  pix_->copy_points(*pix_, &points);
  ASSERT_EQ(points.size(), pix_->n_points());
  pix_->copy_points(pix_->parent(), &points);
  ASSERT_EQ(points.size(), pix_->n_points());

  int child_points = 0;
  s2omp::pixel child = pix_->child_begin();
  for (int k = 0; k < random_points.size(); k++) {
    if (child.contains(random_points[k])) {
      child_points++;
    }
  }
  pix_->copy_points(child, &points);
  ASSERT_EQ(points.size(), child_points);
}

TEST_F(TreePixelTest, TestTreePixelPairCounts) {
  // Test the method for counting the number of points within a tree_pixel that
  // are within an angular bin from a reference point.

  // Start with an angular bin that contains the entire reference tree_pixel.
  double scale = sqrt(s2omp::pixel::average_area(pix_->level()));
  s2omp::angular_bin bin(0.0, 10.0 * scale);
  s2omp::point p = pix_->next().center_point();
  s2omp::annulus_bound* bound = s2omp::annulus_bound::from_angular_bin(p, bin);
  ASSERT_TRUE(bound->contains(*pix_));

  // Now add points to the reference pixel and verify that all points are within
  // the bound.
  s2omp::timer timer;
  int n_points = 100 * pix_->node_capacity();
  s2omp::point_vector points;
  pix_->get_random_points(n_points, &points);
  for (int k = 0; k < points.size(); k++) {
    ASSERT_TRUE(pix_->add_point(points[k]));
  }
  timer.start_timer();
  ASSERT_EQ(points.size(), pix_->find_pairs(*bound));
  timer.stop_timer();
  double full_pixel_time = timer.elapsed_time();
  ASSERT_NEAR(1.0 * points.size(), pix_->find_weighted_pairs(*bound), 1.0e-10);

  // Now change the bound so that it doesn't enclose the entire pixel and
  // verify that the pair counts match the counts from a direct check.
  delete bound;
  bin.set_theta_min(scale);
  bin.set_theta_max(1.5 * scale);
  bound = s2omp::annulus_bound::from_angular_bin(p, bin);
  ASSERT_FALSE(bound->contains(*pix_));

  timer.start_timer();
  int n_contained = 0;
  for (int k = 0; k < points.size(); k++) {
    if (bound->contains(points[k])) {
      n_contained++;
    }
  }
  timer.stop_timer();
  double direct_time = timer.elapsed_time();
  timer.start_timer();
  ASSERT_EQ(n_contained, pix_->find_pairs(*bound));
  timer.stop_timer();
  std::cout <<
      "full pixel: " << full_pixel_time * 1000.0 << " msec; " <<
      "raw counts: " << direct_time * 1000.0 << " msec; " <<
      "find_pairs: " << timer.elapsed_time() * 1000.0 << " msec\n";
}

TEST_F(TreePixelTest, TestTreePixelNearestNeighbor) {
  // Test the methods for finding the nearest neighbor to an input point.

  // Start by adding some random points to a tree_pixel.
  int n_points = 1000;
  s2omp::point_vector points;
  pix_->get_random_points(n_points, &points);

  for (int k = 0; k < points.size(); k++) {
    ASSERT_TRUE(pix_->add_point(points[k]));
  }

  // Now try finding the nearest neighbor to one of the points
  s2omp::point test_point = points[n_points/2];
  s2omp::point p = pix_->find_nearest_neighbor(test_point);
  ASSERT_EQ(p.id(), test_point.id());

  // Now find the nearest neighbor to a point outside the tree_pixel.  Verify
  // that the distance between the reference point and all other points in the
  // pixel is larger than that for the returned point.
  s2omp::point outside = pix_->next().center_point();
  p = pix_->find_nearest_neighbor(outside);
  double min_ang_dist = p.angular_distance_deg(outside);
  for (int k = 0; k < points.size(); k++) {
    if (points[k] != p) {
      ASSERT_GT(points[k].angular_distance_deg(outside), min_ang_dist);
    }
  }
}

TEST_F(TreePixelTest, TestTreePixelClosestMatch) {
  // Test the methods for finding the closest match to an input point.

  // Start by adding some random points to a tree_pixel.
  int n_points = 1000;
  s2omp::point_vector points;
  pix_->get_random_points(n_points, &points);

  for (int k = 0; k < points.size(); k++) {
    ASSERT_TRUE(pix_->add_point(points[k]));
  }

  // Now try finding the closest match to one of the points
  double tol_deg = 1.0e-10;
  s2omp::point test_point = points[n_points/2];
  s2omp::point match;
  ASSERT_TRUE(pix_->closest_match(test_point, tol_deg, match));
  ASSERT_LT(test_point.angular_distance_deg(match), tol_deg);

  // Now find the closest match to a point outside the pixel.  With the same
  // tolerance, this should fail.  Likewise, the input point should remain
  // unchanged.
  s2omp::point outside = pix_->next().center_point();
  ASSERT_FALSE(pix_->closest_match(outside, tol_deg, match));
  ASSERT_LT(test_point.angular_distance_deg(match), tol_deg);

  // If we increase the allowable distance to the distance between the two
  // pixels, we should find an acceptable match.
  tol_deg = pix_->center_point().angular_distance_deg(outside);
  ASSERT_TRUE(pix_->closest_match(outside, tol_deg, match));
  ASSERT_LT(outside.angular_distance_deg(match), tol_deg);

  // Verify that the returned point is, in fact, closer to the input point than
  // any other point in the pixel.
  double min_ang_dist = match.angular_distance_deg(outside);
  for (int k = 0; k < points.size(); k++) {
    if (points[k] != match) {
      ASSERT_GT(points[k].angular_distance_deg(outside), min_ang_dist);
    }
  }
}
