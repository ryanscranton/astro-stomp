/*
 * pixel_test.cc
 *
 *  Created on: Aug 19, 2012
 *      Author: cbmorrison
 */

#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include "pixel.h"
#include "point.h"

using namespace s2omp;

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
  std::cout << "\t\tPixel Leaf: id=" << hex << pix_leaf->id() <<
      "; level=" << dec << pix_leaf->level() <<
      "; is_leaf=" << pix_leaf->is_leaf() <<
      "; is_face=" << pix_leaf->is_face() << "\n";
  std::cout << "\t\tPixel Mid: id=" << hex << pix_mid->id() <<
      "; level=" << dec << pix_mid->level() << "\n";
  std::cout << "\t\tPixel Face: id=" << hex << pix_face->id() <<
      "; level=" << dec << pix_face->level() << "\n";
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
  pixel_vector pixels;

  // Now we fill the pixels up with a resonable level so we can test that it is
  // the right length.
  pix->children(&pixels);
  if (pixels.size() == 4) {
    std::cout << "\t\tSuccessfully created child pixels.\n";
  } else {
    std::cout << "\t\tFailed to create expected number of child pixels.\n"
        << "\t\tCheck children method in pixel.";
    std::cout << "\t\tExpected 4, Got " << pixels.size() << "\n";
  }

  // Now we fill the pixels up with a resonable level so we can test that it is
  // the right length.
  pix->children(5, &pixels);
  if (pixels.size() == 1024) {
    std::cout << "Successfully created child pixels.\n";
  } else {
    std::cout << "\t\tFailed to create expected number of child pixels.\n"
              << "\t\tCheck children method in pixel.\n";
    std::cout << "\t\tExpected 1024, Got " << pixels.size() << "\n";
  }

  std::cout << "\tTesting Parents...\n";
  // Now we test so see if you can successfully create parent pixels.
	delete pix;
  pixel* pix2 = pixel::from_point(point(1, 0, 0, 0), 30);
  pixel parent = pix2->parent();
  std::cout << "\t\tParent level is " << parent.level();
  std::cout << "\t\tParent id | child id:"
      << parent.id() << " | " << pix2->id() << "\n";

  pixel parent2 = pix2->parent(10);
  std::cout << "\t\tParent level is " << parent2.level();
  std::cout << "\t\tParent id | child id:"
      << parent2.id() << " | " << pix2->id() << "\n";
  delete pix2;

  std::cout << "\tTesting Child Iteration...\n";
  pixel* pix3 = pixel::from_point(point(1, 0, 0, 0), 15);
  int n_child = 0;
  pixel end = pix3->child_end();
  pixel child = pix3->child_begin();
  while (child != end) {
    n_child++;
    child = child.next();
  }
  if (n_child == 4) {
    std::cout << "\t\tSuccessfully created 4 children.\n";
  } else {
    std::cout << "\t\tFailed to create 4 children."
        << "Check next and child_end/begin.\n";
    std::cout << "\t\tExpected 4, Got " << n_child << "\n";
  }

  std::cout << "\tTesting Neighbors...\n";
  pix3->neighbors(&pixels);
  if (pixels.size() == 8) {
    std::cout << "\t\tSuccessfully created 8 neighbors.";
  } else {
    std::cout << "\t\tFailed to create 8 children. "
        << "Check next and child_end/begin.\n"
        << "Expected 8, Got " << pixels.size() << "\n";
  }
}

void test_pixel_geometry() {
  std::cout << "Testing Pixel Geometry...\n";
  std::cout << "\tTesting Pixel Area...\n";
  pixel* pix = pixel::from_point(point(0, 1, 0, 0), 15);
  std::cout << "\t\tAverage are at level 15: " << pix->average_area() <<
      " 14: " << pix->average_area(14) <<
      " 16: " << pix->average_area(16) << "\n";

  std::cout << "\tTesting Pixel Contains/Intersect...\n";
  pixel child_begin = pix->child_begin(17);
  if (pix->contains(point(0, 1, 0, 0))) {
    std::cout << "\t\tSuccessfully Contains point.\n";
  } else {
    std::cout << "\t\tFaield on contains point.\n";
  }
  if (pix->contains(child_begin)) {
    std::cout << "\t\tSuccessfully Contains point.\n";
  } else {
    std::cout << "\t\tFaield on contains pixel.\n";
  }
  pixel parent = pix->parent(4);
  if (pix->intersects(parent)) {
    std::cout << "\t\tSuccessfully Contains point.\n";
  } else {
    std::cout << "\t\tFaield on intersect pixel.\n";
  }

  std::cout << "\tTesting Pixel Vertices/Edges...\n";
  point_vector vertices;
  point_vector edges;
  for (int k = 0; k < 4; ++k) {
    vertices.push_back(pix->vertex(k));
    edges.push_back(pix->edge(k));
  }

  if (vertices.size() == 4) {
    std::cout << "\t\tSuccessfully Created 4 Vertices.\n";
  } else {
    std::cout << "\t\tFailed on create vertices.\n";
  }

  if (edges.size() == 4) {
    std::cout << "\t\tSuccessfully Created 4 Edges.\n";
  } else {
    std::cout << "\t\tFailed on create edges.\n";
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

int main(int argc, char **argv) {
	pixel_unit_tests();

	return 0;
}
