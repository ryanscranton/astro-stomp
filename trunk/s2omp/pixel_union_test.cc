#include <gtest/gtest.h>

#include "pixel_union.h"

#include "point.h"

TEST(pixel_union, TestPixelUnionDefaultConstructor) {
  // The default constructor without any arguments should produce an empty
  // pixel_union.
  s2omp::pixel_union pix_union;

  ASSERT_TRUE(pix_union.is_empty());
  ASSERT_EQ(pix_union.size(), 0);
  ASSERT_DOUBLE_EQ(pix_union.area(), 0.0);

  // This shouldn't contain anything.
  s2omp::point p(1.0, 0.0, 0.0, 1.0);
  s2omp::pixel pix = p.to_pixel();
  ASSERT_FALSE(pix_union.contains(p));
  ASSERT_FALSE(pix_union.contains(pix));
}

void clean_up_union(s2omp::pixel_union* pix_union) {
  pix_union->clear();
  delete pix_union;
}

TEST(pixel_union, TestPixelUnionInitialization) {
  // Test the basic initialization routines.

  // If we make a pixel_vector out of the children of a pixel, then the
  // resulting pixel_union should be the original pixel.
  int level = 20;
  int child_level = level + 5;
  s2omp::point p(1.0, 0.0, 0.0, 1.0);
  s2omp::pixel pix = p.to_pixel(level);
  s2omp::pixel_vector children;
  pix.children(child_level, &children);
  ASSERT_EQ(children.size(), 1024);

  s2omp::pixel_union* pix_union = s2omp::pixel_union::from_covering(children);
  ASSERT_FALSE(pix_union->is_empty());
  ASSERT_TRUE(children.empty());
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(p));
  ASSERT_TRUE(pix_union->contains(pix));
  ASSERT_NEAR(pix_union->contained_area(pix), pix.exact_area(), 1.0e-10);
  pix_union->clear();

  // Now, create a new union where we exclude one of the child pixels.  This
  // should result in a pixel_union with more than one internal pixel.
  pix.children(&children);
  s2omp::pixel removed_pix = children.back();
  children.pop_back();
  // Initializing the pixel_union will empty the input pixel_vector, so we need
  // to grab any pixels we'll use later before we initialize.
  s2omp::pixel first_child = children[0];
  pix_union->init(children);

  ASSERT_FALSE(pix_union->is_empty());
  ASSERT_EQ(pix_union->size(), 3);
  ASSERT_FALSE(pix_union->contains(pix));
  ASSERT_FALSE(pix_union->contains(removed_pix));
  ASSERT_TRUE(pix_union->contains(first_child));
  ASSERT_NEAR(pix_union->contained_area(first_child), first_child.exact_area(),
      1.0e-10);
  ASSERT_NEAR(pix_union->area(), pix.exact_area() - removed_pix.exact_area(),
      1.0e-10);

  // Now verify with a case spanning multiple resolution levels.
  pix.children(child_level, &children);
  children.pop_back();
  pix_union->init(children);
  // Expect 3 pixels at each level that the union spans.
  ASSERT_EQ(pix_union->size(), 3 * (child_level - level));
  ASSERT_EQ(pix_union->max_level(), child_level);
  ASSERT_EQ(pix_union->min_level(), level + 1);
  ASSERT_TRUE(pix_union->contains(p));
  ASSERT_FALSE(pix_union->contains(pix));

  clean_up_union(pix_union);
}

TEST(pixel_union, TestPixelUnionSoften) {
  // Test the routine for softening the pixel_union to a specified level.

  // If we make a pixel_union out of the children of a given pixel excluding
  // one child, then soften to the original pixel's level, we should end up
  // with a pixel union equivalent to our original pixel.
  int level = 20;
  s2omp::point p(1.0, 0.0, 0.0, 1.0);
  s2omp::pixel pix = p.to_pixel(level);

  s2omp::pixel_vector children;
  pix.children(&children);
  s2omp::pixel removed_pix = children.back();
  children.pop_back();

  s2omp::pixel_union* pix_union = s2omp::pixel_union::from_covering(children);
  ASSERT_EQ(pix_union->size(), 3);
  ASSERT_FALSE(pix_union->contains(pix));
  ASSERT_FALSE(pix_union->contains(removed_pix));

  pix_union->soften(level);
  ASSERT_EQ(pix_union->max_level(), level);
  ASSERT_EQ(pix_union->min_level(), level);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));
  ASSERT_TRUE(pix_union->contains(removed_pix));

  // Now verify with a case spanning multiple resolution levels.
  int child_level = level + 5;
  pix.children(child_level, &children);
  children.pop_back();
  pix_union->init(children);

  int soften_level = child_level - 1;
  // Even though we're softening a single level, adding the parent of the
  // missing child pixel is enough to collapse us back to the original pixel.
  pix_union->soften(soften_level);
  ASSERT_EQ(pix_union->max_level(), level);
  ASSERT_EQ(pix_union->min_level(), level);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));
  ASSERT_TRUE(pix_union->contains(removed_pix));

  clean_up_union(pix_union);
}

s2omp::pixel_union* create_union_from_pixel(s2omp::pixel pix,
    bool include_neighbors) {
  s2omp::pixel_vector pixels;
  if (include_neighbors) {
    pix.neighbors(&pixels);
  }
  pixels.push_back(pix);
  return s2omp::pixel_union::from_covering(pixels);
}

TEST(pixel_union, TestPixelUnionCombine) {
  // Test the methods for combining pixel_unions.

  // Start by verifying that combining a pixel_union with itself leaves the
  // union unchanged.
  int level = 20;
  s2omp::point p(1.0, 0.0, 0.0, 1.0);
  s2omp::pixel pix = p.to_pixel(level);
  s2omp::pixel_union* pix_union = create_union_from_pixel(pix, false);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));

  pix_union->combine_with(*pix_union);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));

  // Create a pixel_union from our starting pixel's nearest neighbor.
  // Verify that the new union has size == 2 and contains both pixels.
  s2omp::pixel neighbor = pix.next();
  s2omp::pixel_union* neighbor_union = create_union_from_pixel(neighbor, false);

  pix_union->combine_with(*neighbor_union);
  ASSERT_EQ(pix_union->size(), 2);
  ASSERT_TRUE(pix_union->contains(pix));
  ASSERT_TRUE(pix_union->contains(neighbor));
  ASSERT_NEAR(pix_union->area(), pix.exact_area() + neighbor.exact_area(),
      1.0e-10);

  // Verify that we get the same behavior from the init_from_combination
  // method.
  s2omp::pixel_union combined_union;
  combined_union.init_from_combination(*pix_union, *neighbor_union);
  ASSERT_EQ(combined_union.size(), 2);
  ASSERT_TRUE(combined_union.contains(pix));
  ASSERT_TRUE(combined_union.contains(neighbor));
  ASSERT_NEAR(combined_union.area(), pix_union->area(), 1.0e-10);

  clean_up_union(pix_union);
  clean_up_union(neighbor_union);
}

TEST(pixel_union, TestPixelUnionIntersection) {
  // Test the various methods for finding intersections with pixel_unions.

  // Start by verifying that the intersection of a union with itself leaves
  // the union unchanged;
  int level = 20;
  s2omp::point p(1.0, 0.0, 0.0, 1.0);
  s2omp::pixel pix = p.to_pixel(level);
  s2omp::pixel_union* pix_union = create_union_from_pixel(pix, false);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));

  pix_union->intersect_with(*pix_union);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));

  // Now take the intersection of the union and a child of the original pixel.
  // The result should be equivalent to the child pixel.
  s2omp::pixel_union* child_union =
      create_union_from_pixel(pix.child_begin(), false);
  ASSERT_EQ(child_union->size(), 1);
  ASSERT_TRUE(child_union->contains(pix.child_begin()));
  ASSERT_FALSE(child_union->contains(pix));
  ASSERT_TRUE(pix_union->may_intersect(*child_union));
  ASSERT_TRUE(pix_union->intersects(*child_union));

  s2omp::pixel_union intersect_union;
  intersect_union.init_from_intersection(*pix_union, *child_union);
  ASSERT_EQ(intersect_union.size(), 1);
  ASSERT_TRUE(intersect_union.contains(pix.child_begin()));
  ASSERT_FALSE(intersect_union.contains(pix));

  // Finally, verify that the intersection between two unions that don't
  // intersect is empty.
  s2omp::pixel_union* outside_union =
      create_union_from_pixel(pix.next(), false);
  intersect_union.init_from_intersection(*pix_union, *outside_union);
  ASSERT_TRUE(intersect_union.is_empty());
  ASSERT_EQ(intersect_union.size(), 0);

  // Now consider the methods for intersecting pixels and unions.  First,
  // verify that the intersection between a child of the reference pixel and
  // the reference union is the child pixel.
  s2omp::pixel_vector pixels;
  pix_union->pixel_intersection(pix.child_begin(), &pixels);
  ASSERT_EQ(pixels.size(), 1);
  ASSERT_TRUE(pixels[0] == pix.child_begin());

  // Alternatively, if we construct a union from a pixel missing a child, then
  // the intersection between that union and the pixel should be equivalent to
  // the original input pixel_vector, but compacted by the pixel_union.
  int child_level = level + 5;
  s2omp::pixel_vector children;
  pix.children(child_level, &children);
  children.pop_back();
  pix_union->init(children);
  pix_union->pixel_intersection(pix, &pixels);
  ASSERT_EQ(pix_union->size(), pixels.size());

  clean_up_union(pix_union);
  clean_up_union(child_union);
  clean_up_union(outside_union);
}

TEST(pixel_union, TestPixelUnionExclusion) {
  // Test the various methods for finding exclusions with pixel_unions.

  // Start by verifying that excluding a non-intersecting union leaves our
  // reference union unchanged.
  int level = 20;
  s2omp::point p(1.0, 0.0, 0.0, 1.0);
  s2omp::pixel pix = p.to_pixel(level);
  s2omp::pixel_union* pix_union = create_union_from_pixel(pix, false);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));

  s2omp::pixel_union* outside_union =
      create_union_from_pixel(pix.next(), false);
  pix_union->exclude_from(*outside_union);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));

  // Now exclude a union constructed from one of the child pixels.  The result
  // should be a union consisting of the other three children.
  s2omp::pixel_union* child_union =
      create_union_from_pixel(pix.child_begin(), false);
  s2omp::pixel_union exclude_union;
  exclude_union.init_from_exclusion(*pix_union, *child_union);
  ASSERT_EQ(exclude_union.size(), 3);
  ASSERT_FALSE(exclude_union.contains(pix));
  ASSERT_FALSE(exclude_union.contains(pix.child_begin()));
  ASSERT_TRUE(exclude_union.intersects(pix));
  ASSERT_FALSE(exclude_union.intersects(pix.child_begin()));
  for (s2omp::pixel c = pix.child_begin().next();
      c != pix.child_end(); c = c.next()) {
    ASSERT_TRUE(exclude_union.contains(c));
  }

  // Finally, exclude a union that contains the current union and verify that
  // the resulting union is empty.
  s2omp::pixel_union* parent_union =
      create_union_from_pixel(pix, true);
  ASSERT_TRUE(parent_union->contains(pix));
  // The initial pixel_vector would have contained 9 pixels.  4 of those,
  // including the center pixel, should be combined into the parent of the
  // original pixel, leaving 6 pixels in total.
  ASSERT_TRUE(parent_union->contains(pix.parent()));
  ASSERT_EQ(parent_union->size(), 6);

  exclude_union.init_from_exclusion(*pix_union, *parent_union);
  ASSERT_TRUE(exclude_union.is_empty());

  // Now consider the methods for finding the area in an input pixel that is
  // outside an input pixel.  Start with the case where the input pixel is
  // equivalent to the union itself.
  s2omp::pixel_vector pixels;
  pix_union->pixel_exclusion(pix, &pixels);
  ASSERT_TRUE(pixels.empty());

  // Now try a parent of the input pixel.  The number of output pixels should
  // be a function of the difference in level.
  int parent_level = level - 5;
  pix_union->pixel_exclusion(pix.parent(parent_level), &pixels);
  ASSERT_EQ(pixels.size(), 3 * (level - parent_level));
  for (s2omp::pixel_iterator iter = pixels.begin();
      iter != pixels.end(); ++iter) {
    ASSERT_FALSE(pix_union->contains(*iter));
    ASSERT_FALSE(pix_union->intersects(*iter));
  }

  clean_up_union(pix_union);
  clean_up_union(child_union);
  clean_up_union(outside_union);
  clean_up_union(parent_union);
}

TEST(pixel_union, TestPixelUnionBound) {
  // Test the virtual methods for finding the circle_bound enclosing the
  // pixel_union and its center.

  // Start with a pixel_bound consisting of a single pixel.
  int level = 18;
  s2omp::point p(0.0, 0.0, 1.0, 1.0);
  s2omp::pixel pix = p.to_pixel(level);
  s2omp::pixel_union* pix_union = create_union_from_pixel(pix, false);
  ASSERT_EQ(pix_union->size(), 1);
  ASSERT_TRUE(pix_union->contains(pix));
  ASSERT_TRUE(pix_union->contains(p));

  s2omp::circle_bound bound = pix_union->get_bound();
  ASSERT_FALSE(bound.is_empty());
  ASSERT_TRUE(bound.is_valid());
  ASSERT_TRUE(s2omp::double_ge(bound.radius(), pix.get_bound().radius()));
  ASSERT_TRUE(pix.get_bound().contains(p));
  // TODO(scranton): This test passes, but there are some numerical issues here
  // for finer pixels wherein the pixel_union bound may not contain the point
  // that the original pixel was based on even if the bound for the
  // corresponding pixel does.  I *think* this is is just flakiness based on
  // the way that we're instantiating the point and that this won't be a real
  // problem, but it's worth keeping an eye on.
  ASSERT_TRUE(bound.contains(p));
  ASSERT_TRUE(bound.contains(pix));

  s2omp::point center = pix_union->get_center();
  ASSERT_TRUE(bound.contains(center));
  // This isn't generically true, but should be true for our case.
  ASSERT_TRUE(pix_union->contains(center));

  // Now verify that this also works for a pixel_union consisting of multiple
  // pixels at different resolutions.
  s2omp::pixel_union* neighbor_union = create_union_from_pixel(pix, true);
  bound = neighbor_union->get_bound();
  ASSERT_FALSE(bound.is_empty());
  ASSERT_TRUE(bound.is_valid());
  ASSERT_TRUE(bound.contains(p));
  ASSERT_TRUE(bound.contains(pix));
  for (s2omp::pixel_iterator iter = neighbor_union->begin();
      iter != neighbor_union->end(); ++iter) {
    ASSERT_TRUE(bound.contains(*iter));
  }

  center = pix_union->get_center();
  ASSERT_TRUE(bound.contains(center));
  // This isn't generically true, but should be true for our case.
  ASSERT_TRUE(neighbor_union->contains(center));
}
