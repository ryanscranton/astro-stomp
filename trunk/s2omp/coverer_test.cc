#include <gtest/gtest.h>

#include "coverer.h"

#include "circle_bound.h"
#include "pixel_union.h"

TEST(coverer, TestCovererBasics) {
  // The default constructor without any arguments should produce a coverer
  // with min and max levels set to the full range of possible values.
  s2omp::coverer cover;
  ASSERT_EQ(cover.min_level(), 0);
  ASSERT_EQ(cover.max_level(), s2omp::MAX_LEVEL);

  // From here, we can set the minimum and maximum levels.
  int min_level = 10;
  int max_level = 20;
  ASSERT_TRUE(cover.set_min_level(min_level));
  ASSERT_TRUE(cover.set_max_level(max_level));
  ASSERT_EQ(cover.min_level(), min_level);
  ASSERT_EQ(cover.max_level(), max_level);

  ASSERT_TRUE(cover.set_min_max_level(min_level, min_level));
  ASSERT_EQ(cover.min_level(), min_level);
  ASSERT_EQ(cover.max_level(), min_level);

  // If we try to set the levels to invalid values (min > max), then the
  // return value should be false and the level should remain unchanged.
  ASSERT_FALSE(cover.set_min_max_level(max_level, min_level));
  ASSERT_EQ(cover.min_level(), min_level);
  ASSERT_EQ(cover.max_level(), min_level);

  // Alternatively, we can initialize with min and max level values.
  cover = s2omp::coverer(min_level, max_level);
  ASSERT_EQ(cover.min_level(), min_level);
  ASSERT_EQ(cover.max_level(), max_level);
}

void profile_covering(std::string tag, s2omp::pixel_vector& covering,
    s2omp::bound_interface& bound) {
  std::cout << tag << " profile:\n";
  double area = 0.0;
  int min_level = s2omp::MAX_LEVEL;
  int max_level = 0;
  int n_contained = 0;
  int n_center_contained = 0;
  for (s2omp::pixel_iterator iter = covering.begin();
      iter != covering.end(); ++iter) {
    area += iter->exact_area();
    int level = iter->level();
    if (level < min_level) min_level = level;
    if (level > max_level) max_level = level;
    if (bound.contains(*iter)) n_contained++;
    if (bound.contains(iter->center_point())) n_center_contained++;
  }
  std::cout << "\tn_pixel = " << covering.size() << ", n_contained = " <<
      n_contained << ", n_center_contained = " << n_center_contained << "\n";
  std::cout << "\tArea: " << area << " sq. deg.\n";
  std::cout << "\t" << min_level << " <= level <= " << max_level << "\n";
  std::cout << "\t\tBound area: " << bound.area() << " sq. deg.\n";
}

TEST(coverer, TestCovererSimpleCovering) {
  // Test the method for generating a simple covering of a
  // bound_interface-derived object.

  // Start by constructing a simple circle_bound.
  s2omp::point axis = s2omp::point::from_radec_deg(0.0, 0.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Choose a starting resolution that's close to the area of the bound:
  int level = s2omp::pixel::get_level_from_area(bound->area());

  // Create a simple covering on that scale.
  s2omp::pixel_vector covering;
  s2omp::coverer::get_simple_covering(*bound, level, &covering);
  profile_covering("SimpleCovering", covering, *bound);
  ASSERT_FALSE(covering.empty());
  for (int k = 0; k < covering.size(); k++) {
    ASSERT_TRUE(bound->may_intersect(covering[k]));
    ASSERT_TRUE(s2omp::double_le(bound->contained_area(covering[k]),
        covering[k].exact_area()));
    ASSERT_TRUE(s2omp::double_ge(bound->contained_area(covering[k]), 0.0));
    // Covering pixels should be sorted.
    if (k > 0) {
      ASSERT_TRUE(covering[k - 1] < covering[k]);
    }
  }

  // Verify that none of the neighbors of our covering pixels intersect the
  // bound, unless they're in the covering as well.
  for (int k = 0; k < covering.size(); k++) {
    s2omp::pixel_vector neighbors;
    covering[k].neighbors(&neighbors);
    for (int j = 0; j < neighbors.size(); j++) {
      if (!std::binary_search(covering.begin(), covering.end(), neighbors[j])) {
        ASSERT_FALSE(bound->may_intersect(neighbors[j]));
      }
    }
  }

  // If we make a pixel_union from the covering, then we can verify that the
  // simple covering does contain the entire bound.
  s2omp::pixel_union covering_union;
  int n_covering = covering.size();
  covering_union.init(covering);
  ASSERT_GT(covering_union.area(), bound->area());
  ASSERT_LE(covering_union.size(), n_covering);
  int n_random = 1000;
  for (int k = 0; k < n_random; k++) {
    ASSERT_TRUE(covering_union.contains(bound->create_random_point()));
  }

  // Now increase the resolution used for the covering.
  level += 5;
  s2omp::pixel_vector fine_covering;
  s2omp::coverer::get_simple_covering(*bound, level, &fine_covering);
  profile_covering("FineSimpleCovering", fine_covering, *bound);
  ASSERT_FALSE(fine_covering.empty());
  ASSERT_GT(fine_covering.size(), n_covering);
  for (int k = 0; k < fine_covering.size(); k++) {
    ASSERT_TRUE(bound->may_intersect(fine_covering[k]));
    ASSERT_TRUE(s2omp::double_le(bound->contained_area(fine_covering[k]),
        fine_covering[k].exact_area()));
    ASSERT_TRUE(s2omp::double_ge(bound->contained_area(fine_covering[k]), 0.0));
  }

  // Re-initialize the pixel_union and run the same verification tests.
  int n_fine_covering = fine_covering.size();
  covering_union.init(fine_covering);
  ASSERT_GT(covering_union.area(), bound->area());
  ASSERT_LE(covering_union.size(), n_fine_covering);
  for (int k = 0; k < n_random; k++) {
    ASSERT_TRUE(covering_union.contains(bound->create_random_point()));
  }
}

TEST(coverer, TestCovererCenterCovering) {
  // Test the method for generating a center covering of a
  // bound_interface-derived object.

  // Start by constructing a simple circle_bound.
  s2omp::point axis = s2omp::point::from_radec_deg(0.0, 0.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Choose a starting resolution that's close to the area of the bound:
  int level = s2omp::pixel::get_level_from_area(bound->area());

  // Create a center covering on that scale.
  s2omp::pixel_vector covering;
  s2omp::coverer::get_center_covering(*bound, level, &covering);
  profile_covering("CenterCovering", covering, *bound);
  ASSERT_FALSE(covering.empty());
  for (int k = 0; k < covering.size(); k++) {
    ASSERT_TRUE(bound->contains(covering[k].center_point()));
  }

  // Verify that all of the pixels in the center covering would also be
  // included in a simple covering and that the number of center covering
  // pixels is <= the number of simple covering cells.
  s2omp::pixel_vector simple_covering;
  s2omp::coverer::get_simple_covering(*bound, level, &simple_covering);
  ASSERT_GE(simple_covering.size(), covering.size());
  for (int k = 0; k < covering.size(); k++) {
    ASSERT_TRUE(std::binary_search(simple_covering.begin(),
        simple_covering.end(), covering[k]));
  }

  // Now increase the resolution used for the covering.
  level += 5;
  s2omp::pixel_vector fine_covering;
  s2omp::coverer::get_center_covering(*bound, level, &fine_covering);
  profile_covering("FineCenterCovering", fine_covering, *bound);
  ASSERT_FALSE(fine_covering.empty());
  ASSERT_GT(fine_covering.size(), covering.size());
  for (int k = 0; k < fine_covering.size(); k++) {
    ASSERT_TRUE(bound->contains(fine_covering[k].center_point()));
  }
}

TEST(coverer, TestCovererStandardCovering) {
  // Test the method for generating a multi-resolution covering with a fixed
  // number of pixels.

  // Start by constructing a simple circle_bound.
  s2omp::point axis = s2omp::point::from_radec_deg(0.0, 0.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Create a coverer with level bounds that are reasonable, given the bound's
  // area, and generate a covering with the default number of pixels.
  s2omp::coverer cover;
  cover.set_levels_from_area(bound->area());
  s2omp::pixel_vector covering;
  ASSERT_TRUE(cover.get_covering(*bound, &covering));
  profile_covering("StandardCovering", covering, *bound);

  ASSERT_FALSE(covering.empty());
  ASSERT_LE(covering.size(), s2omp::DEFAULT_COVERING_PIXELS);

  // If we make a pixel_union from the covering, then we can verify that the
  // covering does contain the entire bound.
  s2omp::pixel_union covering_union;
  int n_covering = covering.size();
  covering_union.init(covering);
  ASSERT_GT(covering_union.area(), bound->area());
  ASSERT_LE(covering_union.size(), n_covering);
  int n_random = 1000;
  for (int k = 0; k < n_random; k++) {
    ASSERT_TRUE(covering_union.contains(bound->create_random_point()));
  }

  // Now verify that, if we increase the allowed number of pixels, we can get
  // a closer approximation to the starting bound.
  long max_pixels = 10 * s2omp::DEFAULT_COVERING_PIXELS;
  ASSERT_TRUE(cover.get_size_covering(max_pixels, *bound, &covering));
  profile_covering("FineStandardCovering", covering, *bound);
  ASSERT_GT(covering.size(), s2omp::DEFAULT_COVERING_PIXELS);
  ASSERT_LE(covering.size(), max_pixels);
  double area = 0.0;
  for (int k = 0; k < covering.size(); k++) {
    area += covering[k].exact_area();
  }
  ASSERT_LT(area - bound->area(), covering_union.area() - bound->area());
}

TEST(coverer, TestCovererAreaCovering) {
  // Test the method for generating a multi-resolution covering based on
  // matching the covering area to the bound_interface's area.

  // Start by constructing a simple circle_bound.
  s2omp::point axis = s2omp::point::from_radec_deg(0.0, 0.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Create a coverer with level bounds that are reasonable, given the bound's
  // area, and generate an area covering.
  s2omp::coverer cover;
  cover.set_levels_from_area(bound->area());
  double precision = 0.2;
  s2omp::pixel_vector covering;
  EXPECT_TRUE(cover.get_area_covering(precision, *bound, &covering));
  profile_covering("AreaCovering", covering, *bound);

  ASSERT_FALSE(covering.empty());

  // To verify our covering, create a pixel_union.
  s2omp::pixel_union covering_union;
  int n_covering = covering.size();
  covering_union.init(covering);
  ASSERT_LT(covering_union.area(), bound->area() * (1.0 + precision));
  ASSERT_GT(covering_union.area(), bound->area() * (1.0 - precision));
  int n_random = 1000;
  for (int k = 0; k < n_random; k++) {
    ASSERT_TRUE(covering_union.contains(bound->create_random_point()));
  }

  // Now increase the precision and verify that we've also increased the number
  // of pixels in the covering.
  s2omp::pixel_vector fine_covering;
  precision *= 0.1;
  EXPECT_TRUE(cover.get_area_covering(precision, *bound, &fine_covering));
  profile_covering("FineAreaCovering", fine_covering, *bound);
  ASSERT_GE(fine_covering.size(), n_covering);
  double area = 0.0;
  for (int k = 0; k < fine_covering.size(); k++) {
    area += fine_covering[k].exact_area();
  }
  ASSERT_LT(abs(area - bound->area()),
      abs(covering_union.area() - bound->area()));
}

TEST(coverer, TestCovererInteriorCovering) {
  // Test the method for generating a multi-resolution interior covering.

  // Start by constructing a simple circle_bound.
  s2omp::point axis = s2omp::point::from_radec_deg(0.0, 0.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Create a coverer with level bounds that are reasonable, given the bound's
  // area, and generate an interior covering.
  s2omp::coverer cover;
  cover.set_levels_from_area(bound->area());
  s2omp::pixel_vector covering;
  ASSERT_TRUE(cover.get_interior_covering(*bound, &covering));
  profile_covering("InteriorCovering", covering, *bound);

  ASSERT_FALSE(covering.empty());

  // Verify that all of the covering pixels are contained by the bound and
  // that the total area is less than the bound's area.
  double area = 0.0;
  for (int k = 0; k < covering.size(); k++) {
    ASSERT_TRUE(bound->contains(covering[k]));
    area += covering[k].exact_area();
  }
  ASSERT_LT(area, bound->area());
}
