#include <gtest/gtest.h>

#include "pixel.h"

#include "circle_bound.h"
#include "point.h"

TEST(pixel, TestPixelDefaultConstructor) {
  // The default constructor without any arguments should produce an invalid
  // pixel.
  s2omp::pixel pix;
  ASSERT_FALSE(pix.is_valid());

  // A pixel with id = 1 should be a valid leaf cell.
  pix = s2omp::pixel(1ULL);
  ASSERT_TRUE(pix.is_valid());
  ASSERT_EQ(pix.id(), 1L);
  ASSERT_EQ(pix.level(), s2omp::MAX_LEVEL);
  ASSERT_TRUE(pix.is_leaf());
  ASSERT_FALSE(pix.is_face());
}

TEST(pixel, TestPixelAltConstructors) {
  // Test constructing a pixel from an input point.
  s2omp::point p(0.0, 0.0, 1.0, 0.0);
  s2omp::pixel pix = s2omp::pixel::from_point(p);

  // This should be a leaf cell
  ASSERT_TRUE(pix.is_leaf());

  // Now use the constructor to generate a coarser resolution pixel.
  int level = 15;
  s2omp::pixel pix_mid = s2omp::pixel::from_point(p, level);
  ASSERT_EQ(pix_mid.level(), level);
  ASSERT_FALSE(pix_mid.is_leaf());
  ASSERT_EQ(pix_mid.id(), pix.parent(level).id());
}

TEST(pixel, TestPixelChildren) {
  // Test generation and containment ofchild pixels.

  // Create a face pixel.
  s2omp::point p(0, 0, 1, 0);
  s2omp::pixel pix = s2omp::pixel::from_point(p, 0);
  s2omp::pixel_vector pixels;

  // Verify that the pixel contains itself.
  ASSERT_TRUE(pix.contains(pix));

  // Generate the 4 immediate children of this pixel.
  pix.children(&pixels);
  ASSERT_EQ(pixels.size(), 4);

  // Verify that all children are contained within the starting pixel but not
  // vice versa
  for (s2omp::pixel_iterator iter = pixels.begin();
      iter != pixels.end(); ++iter) {
    ASSERT_TRUE(iter->is_valid());

    ASSERT_TRUE(pix.contains(*iter));
    ASSERT_FALSE(iter->contains(pix));

    // Intersection should be true both ways.
    ASSERT_TRUE(pix.intersects(*iter));
    ASSERT_TRUE(iter->intersects(pix));
  }

  // Now go a few levels deeper and verify that we've picked up a factor of 4
  // pixels for each level.
  int level = 5;
  pix.children(level, &pixels);
  ASSERT_EQ(pixels.size(), 1024); // 4^5 = 1024
}

TEST(pixel, TestPixelParents) {
  // Test parent creation and containment.

  // Create a leaf pixel
  s2omp::point p(0, 0, 1, 0);
  s2omp::pixel pix = s2omp::pixel::from_point(p, s2omp::MAX_LEVEL);

  // Verify that the immediate parent satisfies the expect level and containment
  // relations
  s2omp::pixel parent = pix.parent();
  ASSERT_TRUE(parent.is_valid());
  ASSERT_FALSE(parent.is_leaf());
  ASSERT_EQ(parent.level(), s2omp::MAX_LEVEL - 1);
  ASSERT_TRUE(parent.contains(pix));
  ASSERT_FALSE(pix.contains(parent));

  // Now go up a few levels and check again.
  int level = 10;
  s2omp::pixel uber_parent = pix.parent(level);
  ASSERT_TRUE(uber_parent.is_valid());
  ASSERT_EQ(uber_parent.level(), level);
  ASSERT_TRUE(uber_parent.contains(pix));
  ASSERT_TRUE(uber_parent.contains(parent));
  ASSERT_EQ(uber_parent.id(), parent.parent(level).id());
}

TEST(pixel, TestPixelIteration) {
  // Verify that pixel id ordering works as expected

  // Create a leaf cell in the middle of a face.
  s2omp::point p(0, 0, 1, 0);
  s2omp::pixel pix = s2omp::pixel::from_point(p, s2omp::MAX_LEVEL);
  ASSERT_LT(pix.id(), pix.next().id());
  ASSERT_GT(pix.id(), pix.prev().id());

  // Reset to a corner pixel and this no longer holds.
  pix = s2omp::pixel(1ULL);
  ASSERT_FALSE(pix.prev().id() == pix.prev_wrap().id());
  ASSERT_GT(pix.prev_wrap().id(), pix.id());
}

TEST(pixel, TestPixelNeighbors) {
  // Test creation of pixel neighbors

  // Create a leaf cell in the middle of a face.
  s2omp::point p(0, 0, 1, 0);
  s2omp::pixel pix = s2omp::pixel::from_point(p, s2omp::MAX_LEVEL);
  s2omp::pixel_vector neighbors;

  // For most pixels, the expected number of neighbors is 8
  pix.neighbors(&neighbors);
  EXPECT_EQ(neighbors.size(), 8);
  for (s2omp::pixel_iterator iter = neighbors.begin();
      iter != neighbors.end(); ++iter) {
    ASSERT_TRUE(iter->is_valid());
    ASSERT_EQ(pix.level(), iter->level());

    ASSERT_FALSE(pix.contains(*iter));
    ASSERT_FALSE(iter->contains(pix));

    ASSERT_FALSE(pix.intersects(*iter));
  }

  // For face pixels, there should only be 4 neighbors.
  s2omp::pixel face = pix.parent(0);
  face.neighbors(&neighbors);
  ASSERT_TRUE(face.is_face());
  EXPECT_EQ(neighbors.size(), 8);

  // If we ask for the neighbors at a finer resolution than our current pixel,
  // then we should expect more of them.  At one level finer, we should get
  // twice as many edge neighbors and an equal number of vertex neighbors,
  // raising the number from 8 to 12.
  s2omp::pixel parent = pix.parent();
  parent.neighbors(s2omp::MAX_LEVEL, &neighbors);
  EXPECT_EQ(neighbors.size(), 12);
}

TEST(pixel, TestPixelBoundBasics) {
  // Testing the basic inherited interface from bound_interface.

  // An invalid pixel should be empty with size == 0.
  s2omp::pixel pix;
  ASSERT_TRUE(pix.is_empty());
  ASSERT_EQ(pix.size(), 0);

  // A valid pixel should not be empty and should have size == 1.
  pix = s2omp::pixel(1);
  ASSERT_FALSE(pix.is_empty());
  ASSERT_EQ(pix.size(), 1);

  // A pixel should contain its center and center_point.
  ASSERT_TRUE(pix.contains(pix.center_point()));
  ASSERT_TRUE(pix.contains(pix.get_center()));

  // Verify some simple contained_area calculations.
  s2omp::pixel parent = pix.parent();
  ASSERT_DOUBLE_EQ(pix.exact_area(), pix.contained_area(pix));
  ASSERT_DOUBLE_EQ(pix.exact_area(), pix.contained_area(parent));
  ASSERT_DOUBLE_EQ(pix.exact_area(), parent.contained_area(pix));

  // After clearing a pixel, it should be empty.
  pix.clear();
  ASSERT_TRUE(pix.is_empty());
  ASSERT_FALSE(pix.is_valid());
}

TEST(pixel, TestPixelCenter) {
  // Test the routines for generating the pixel center.

  // Start with a pixel at the north pole.
  s2omp::point p(0.0, 0.0, 1.0, 1.0);
  s2omp::pixel pix(p.id());

  s2omp::point center = pix.center_point();
  ASSERT_TRUE(center.is_normalized());
  ASSERT_NEAR(p.dot(center), 1.0, 1.0e-5);
  ASSERT_TRUE(pix.contains(center));

  // Now use a point that's well away from the pole.
  p = s2omp::point(1.0, 1.0, 1.0, 1.0);
  pix = s2omp::pixel::from_point(p);
  center = pix.center_point();
  ASSERT_TRUE(center.is_normalized());
  ASSERT_NEAR(p.dot(center), 1.0, 1.0e-5);
  ASSERT_TRUE(pix.contains(center));
}

TEST(pixel, TestPixelBound) {
  // Test the virtual methods for finding the circle_bound enclosing the
  // pixel.

  // Start with a pixel at the north pole.
  int level = 20;
  s2omp::point p(0.0, 0.0, 1.0, 1.0);
  s2omp::pixel pix(p.id(level));

  s2omp::circle_bound bound = pix.get_bound();
  ASSERT_FALSE(bound.is_empty());
  ASSERT_TRUE(bound.is_valid());
  ASSERT_TRUE(bound.contains(p));
  ASSERT_TRUE(bound.contains(pix));
  ASSERT_TRUE(bound.contains(pix.center_point()));
}

TEST(pixel, TestPixelScales) {
  // Not really a test, but rather a delineation of pixel scales.
  for (int level = 0; level <= s2omp::MAX_LEVEL; level++) {
    std::cout << "Level " << level << ": Average area = " <<
        s2omp::pixel::average_area(level) << " sq. deg.; Average scale = " <<
        sqrt(s2omp::pixel::average_area(level)) << " deg.\n";
  }
}
