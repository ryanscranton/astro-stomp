/*
 * pixel_test.cc
 *
 *  Created on: Aug 19, 2012
 *      Author: cbmorrison
 */


#include "circle_bound.h"
#include "core.h"
#include "pixel.h"
#include "point.h"
#include "s2.h"
#include "s2cell.h"

void test_pixel_construction() {
  std::cout << "Testing Pixel Construction...\n";
  std::cout << "\tTesting Default Constructor...\n";
  // If my knowledge of pixel IDs is correct this should give us
  // a leaf pixel on the first face, the "first" quadrant aka the corner.
  pixel pix = pixel(1);

  // For this pixel these tests should return the id(1), level(30),
  // is_leaf(true), is_face(false), and area(really tiny);
  std::cout << "\t\tPixel id=" << pix.id() << ": Should be 1\n"
            << "\t\tPixel level=" << pix.level() << ": Should be 30\n"
            << "\t\tis_leaf() : " << pix.is_leaf() << ": Should be true.\n"
            << "\t\tis_face() : " << pix.is_face() << ": Should be false.\n"
            << "\t\tarea=" << pix.exact_area() << " should be tiny\n";

  // Now that we've tested the default constructor we should test to see if we
  // can properly make a pixel from a point at various levels
  std::cout << "\tTesting from_point Constructor...\n";
  pixel* pix_leaf = pixel::from_point(point(0, 0, 1, 0));
  pixel* pix_mid = pixel::from_point(point(0, 0, 1, 0), 15);
  pixel* pix_face = pixel::from_point(point(0, 0, 1, 0), 0);

  // Even though we are saving the bulk of the geometry tests for later we would
  // still like to see how the pixel centers compare to the input point. If they
  // are in the correct hemisphere and close to the z axis we declare success.
  // To do this we test that the dot product between the centers and the z axis
  // is greater than zero and also close to 1.
  std::cout << "\t\tPixel Leaf: id="
      << pix_leaf->id() << "; level=" << pix_leaf->level() << "; is_leaf="
      << pix_leaf->is_leaf() << "; is_face=" << pix_leaf->is_face() << "\n";
  std::cout << "\t\tPixel Mid: id="
        << pix_mid->id() << "; level=" << pix_mid->level() << "\n";
  std::cout << "\t\tPixel Face: id="
        << pix_face->id() << "; level=" << pix_face->level() << "\n";
  if (pix_leaf->center_point().dot(point(0, 0, 1, 0)) > 0 &&
      pix_mid->center_point().dot(point(0, 0, 1, 0)) > 0 &&
      pix_face->center_point().dot(point(0, 0, 1, 0)) > 0) {
    std::cout << "Succefuly created pixels in the correct quadrant.\n";
  } else {
    std::cout << "Failed to created pixels in the correct quadrant.\n"
        << "\tCheck pixel from point creation, point, or dot product.";
  }
  delete pix_leaf;
  delete pix_mid;
  delete pix_face;
}

void test_pixel_refinement() {
  std::cout << "Testing Pixel Refinement...\n";
  std::cout << "\tTesting Children...\n";
  // Now we test our ablity return both children and parents of the current
  // pixel as well as iterating through the pixels.

  //First we create a face pixel.
  pixel* pix = pixel::from_point(point(1, 0, 0, 0), 0);
  pixel_vector* pixels;

  // Now we fill the pixels up with a resonable level so we can test that it is
  // the right length.
  pix->children(pixels);
  if (pixels->size() == 4) {
    std::cout << "\t\tSuccessfully created child pixels.\n";
  } else {
    std::cout << "\t\tFailed to create expected number of child pixels.\n"
        << "\t\tCheck children method in pixel.";
    std::cout << "\t\tExpected 4, Got " << pixel->size() << "\n";
  }

  // Now we fill the pixels up with a resonable level so we can test that it is
  // the right length.
  pix->children(5, pixels);
  if (pixels->size() == 1024) {
    std::cout << "Successfully created child pixels.\n";
  } else {
    std::cout << "\t\tFailed to create expected number of child pixels.\n"
              << "\t\tCheck children method in pixel.\n";
    std::cout << "\t\tExpected 1024, Got " << pixel->size() << "\n";
  }

  std::cout << "\tTesting Parents...\n";
  // Now we test so see if you can successfully create parent pixels.
  delete pix;
  pixel* pix = pixel::from_point(point(1, 0, 0, 0), 30);
  pixel parent = pix->parent();
  std::cout << "\t\tParent level is " << parent.level();
  std::cout << "\t\tparent id | child id:"
      << parent.id() " | " << pix.id() << "\n";

  pixel parent = pix->parent(10);
  std::cout << "\t\tParent level is " << parent.level();
  std::cout << "\t\tparent id | child id:"
      << parent.id() " | " << pix.id() << "\n";
  delete pix;

  std::cout << "\tTesting Child Iteration...\n";
  pixel* pix = pixel::from_point(point(1, 0, 0, 0), 15);
  int n_child = 0;
  pixel end = pix->child_end();
  pxiel child = pix->child_begin();
  while (child != end) {
    n_child++;
    child.next();
  }
  if (n_child == 4) {
    std::cout << "\t\tSuccessfully created 4 children.\n";
  } else {
    std::cout << "\t\tFailed to create 4 children."
        << "Check next and child_end/begin.\n";
    std::cout << "\t\tExpected 4, Got " << n_child << "\n";
  }

  std::cout << "\tTesting Neighbors...\n";
  pix->neighbors(pixels);
  if (pixels->size() == 4) {
    std::cout << "\t\tSuccessfully created 4 neighbors.";
  } else {
    std::cout << "\t\tFailed to create 4 children. "
        << "Check next and child_end/begin.\n"
        << "Expected 4, Got " << pixels->size() << "\n";
  }
}

void test_pixel_geometry() {
  std::cout << "Testing Pixel Geometry...\n";
  std::cout << "\tTesting Pixel Area...\n";
  pixel* pix = pixel::from_point(point(0, 1, 0, 0));
  std::cout << "\t\tAverage are at level 15: " << pix.average_area() <<
      " 14: " << pix.average_area(14) <<
      " 16: " << pix.average_area(16) << "\n";

  std::cout << "\tTesting Pixel Contains/Intersect...\n";
  pixel child_end = pix->child_end(17);
  if (pix.contains(point(0, 1, 0, 0))) {
    std::cout << "\t\tSuccessfully Contains point.";

  } else {
    std::cout << "\t\Faield on contains point.";
  }
  if (pix->contains(child_end)) {
    std::cout << "\t\tSuccessfully Contains point.";
  } else {
    std::cout << "\t\Faield on contains pixel.";
  }
  pixel parent = pix->parent(4)
  if (pix->intersect(parent)) {
    std::cout << "\t\tSuccessfully Contains point.";
  } else {
    std::cout << "\t\Faield on intersect pixel.";
  }

  std::cout << "\tTesting Pixel Vertices/Edges...\n";
  point_vector vertices;
  ppint_vector edges;
  for (int k = 0; k < 4; ++k) {
    vertices.push_back(pix->vertex(k));
    edges.push_back(pix->edge(k));
  }

  if (vertices.size() == 4) {
    std::cout << "\t\tSuccessfully Created 4 Vertices.";
  } else {
    std::cout << "\t\Failed on create vertices.";
  }

  if (edges.size() == 4) {
    std::cout << "\t\tSuccessfully Created 4 Edges.";
  } else {
    std::cout << "\t\Failed on create edges.";
  }

  std::cout << "\tTesting Pixel Edge Distances...\n";
  std::cout << "\t\tNearest Edge distance, should be about 1 :"
      << pix->nearest_edge_distance(point(1, 0, 0, 0)) << "\n";
  std::cout << "\t\tFarthest Edge distance, should be about 1 :"
        << pix->farthest_edge_distance(point(1, 0, 0, 0)) << "\n";

  if (pix->edge_distances(point(1, 0, 0, 0))) {
    std::cout << "\t\tSuccessful edge distance type\n"
  } else {
    std::cout << "\t\Failed edge distance type\n"
  }
}

void pixel_unit_tests() {
  test_pixel_construction(); // test all of pixels constructor methods and
                             // get/setters.
  test_pixel_refinement(); // test get children and pixel area form the face
                           // pixel to the leaf(and back) as well as accessing
                           // neighbor pixels.
  test_pixel_geometry(); // Test that all pixel methods for accessing
                         // vertices/edges as well as contains and may_intersect

}
