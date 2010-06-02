#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_angular_bin.h"
#include "stomp_angular_correlation.h"
#include "stomp_itree_pixel.h"
#include "stomp_util.h"

void IndexedTreePixelBasicTests() {
  // Before moving on to the map tests, we need to verify that the
  // Stomp::TreePixel, a derived class from Stomp::Pixel, is working properly.
  //  Let's start with some initialization.
  std::cout << "\n";
  std::cout << "************************************\n";
  std::cout << "*** IndexedTreePixel Basic Tests ***\n";
  std::cout << "************************************\n";
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points = 50;
  Stomp::IndexedTreePixel tree_pix(ang, Stomp::HPixResolution, n_points);

  if (tree_pix.Contains(ang)) {
    std::cout << "\tSuccessfully verified the position used to set pixnum.\n";
  } else {
    std::cout << "\tFailed to verify the position used to set pixnum.\n";
  }

  if (tree_pix.AddPoint(ang, 0)) {
    std::cout << "\tSuccessfully added the position used to set the pixnum.\n";
  } else {
    std::cout << "\tFailed to add the position used to set the pixnum.\n";
    exit(1);
  }

  // Now we generate a bunch of random points within the pixel and attempt to
  // add them to the pixel.  This is much greater than the pixel capacity (set
  // above), so we will be testing the generation of sub-pixels as well.
  Stomp::AngularVector angVec;
  tree_pix.GenerateRandomPoints(angVec, 200*n_points);

  std::cout << "\tAttempting to add more points than pixel capacity.\n";
  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if (!tree_pix.AddPoint(*iter, idx)) {
      std::cout << "\t\tFailed to add point: " <<
	iter->Lambda() << ", " << iter->Eta() << "\n";
      exit(1);
    }
    idx++;
  }
  std::cout << "\tAdded " << tree_pix.NPoints() << "/" << 200*n_points <<
    " points to this pixel\n";
}

void IndexedTreePixelPairTests() {
  std::cout << "\n";
  std::cout << "***********************************\n";
  std::cout << "*** IndexedTreePixel Pair Tests ***\n";
  std::cout << "***********************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  uint32_t n_points = 10000;
  Stomp::IndexedTreePixel tree_pix(ang, Stomp::HPixResolution,
				   n_points_per_node);

  Stomp::AngularVector angVec;
  tree_pix.GenerateRandomPoints(angVec, n_points);
  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if (!tree_pix.AddPoint(*iter, idx))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
  }

  std::cout << "Added " << n_points << " points to pixel...\n";

  tree_pix.Ang(ang);
  // First, we choose a radius that is certain to encompass the entire pixel.
  double theta = 15.0;
  Stomp::IndexVector pair_indices;
  tree_pix.FindPairs(ang, theta, pair_indices);
  std::cout << "Found " << pair_indices.size() <<
    " pairs (" << n_points << ")\n";

  // Next, we choose a smaller radius where we should begin to probe some of
  // the tree structure in the pixel.
  theta = 0.1;
  uint32_t direct_pairs = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if (angVec[0].AngularDistance(*iter) <= theta) direct_pairs++;
  }
  tree_pix.FindPairs(angVec[0], theta, pair_indices);
  std::cout << "Found " << pair_indices.size() << " pairs; " <<
    direct_pairs << " pairs from brute force.\n";

  // Now, we check the annulus finding.
  double theta_max = 0.15;
  double theta_min = 0.1;
  double annulus_area =
    (1.0 - cos(theta_max*Stomp::DegToRad))*
    2.0*Stomp::Pi*Stomp::StradToDeg;
  annulus_area -=
    (1.0 - cos(theta_min*Stomp::DegToRad))*
    2.0*Stomp::Pi*Stomp::StradToDeg;
  direct_pairs = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if ((angVec[0].AngularDistance(*iter) <= theta_max) &&
	(angVec[0].AngularDistance(*iter) >= theta_min)) direct_pairs++;
  }
  tree_pix.FindPairs(angVec[0], theta_min, theta_max, pair_indices);
  std::cout << "Found " << pair_indices.size() << " pairs; " <<
    direct_pairs << " pairs from brute force.\n";

  Stomp::AngularBin theta_bin(theta_min, theta_max);

  uint32_t mismatched_pairs = 0;
  uint32_t overcounted_pairs = 0;
  uint32_t undercounted_pairs = 0;
  for (Stomp::AngularIterator ang_iter=angVec.begin();
       ang_iter!=angVec.end();++ang_iter) {
    tree_pix.FindPairs(*ang_iter, theta_bin, pair_indices);
    uint32_t n_pairs = pair_indices.size();
    direct_pairs = 0;
    for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
      if (theta_bin.WithinCosBounds(ang_iter->DotProduct(*iter)))
	direct_pairs++;
    }
    if (n_pairs != direct_pairs) {
      mismatched_pairs++;
      if (n_pairs > direct_pairs) overcounted_pairs++;
      if (n_pairs < direct_pairs) undercounted_pairs++;
      if (mismatched_pairs == 1) ang = *ang_iter;
    }
  }
  std::cout << "Found discrepancy between counts for " << mismatched_pairs <<
    " of " << angVec.size() << " points (" << overcounted_pairs <<
    " over, " << undercounted_pairs << " under).\n";
}

void IndexedTreePixelCoverageTests() {
  // Since we don't know a priori what the underlying distribution of the
  // points added to our TreePixel is, we estimate it using the Coverage
  // method.  This module verifies that Coverage is working properly.
  std::cout << "\n";
  std::cout << "***************************************\n";
  std::cout << "*** IndexedTreePixel Coverage Tests ***\n";
  std::cout << "***************************************\n";
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points = 20;
  Stomp::IndexedTreePixel tree_pix(ang, Stomp::HPixResolution, n_points);

  // Now we find the sub-pixels for the current tree_pix a few levels down
  // and start generating random points from them.
  Stomp::PixelVector pixVec;
  uint32_t hi_resolution = Stomp::HPixResolution*32;
  tree_pix.SubPix(hi_resolution, pixVec);

  Stomp::AngularVector angVec;
  pixVec[0].GenerateRandomPoints(angVec, n_points);
  bool added_point = false;
  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    added_point = tree_pix.AddPoint(*iter, idx);
    if (!added_point)
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
  }

  // At this point our pixel should be full at the top level, so the coverage
  // should still be unity even though we've only been adding points in a
  // small part of it.
  std::cout << "First sub-pixel test:\n";
  std::cout << "\t" << tree_pix.NPoints() << " points: Coverage = " <<
    tree_pix.Coverage() << ", Expected Coverage = 1.0\n";

  // Now, we begin adding more points, triggering the formation of more
  // sub-nodes.  Eventually, we should find that the Coverage equals the
  // ratio of the TreePixel and sub-pixel areas.
  for (int k=0;k<5;k++) {
    pixVec[0].GenerateRandomPoints(angVec, n_points/2);
    bool added_point = false;
    for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
      added_point = tree_pix.AddPoint(*iter, idx);
      if (!added_point)
	std::cout << "\t\tFailed to add point: " <<
	  iter->RA() << ", " << iter->DEC() << "\n";
      idx++;
    }

    std::cout << "\t" << tree_pix.NPoints() << " points: Coverage = " <<
      tree_pix.Coverage() << ", Expected Coverage = " <<
      pixVec[0].Area()/tree_pix.Area() << "\n";
  }

  // Just to be sure, now let's do the same exercise with a sub-pixel on the
  // other side of our TreePixel.
  std::cout << "\nSecond sub-pixel test:\n";
  for (int k=0;k<5;k++) {
    pixVec[pixVec.size()-1].GenerateRandomPoints(angVec, n_points/2);
    bool added_point = false;
    for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
      added_point = tree_pix.AddPoint(*iter, idx);
      if (!added_point)
	std::cout << "\t\tFailed to add point: " <<
	  iter->RA() << ", " << iter->DEC() << "\n";
      idx++;
    }

    std::cout << "\t" << tree_pix.NPoints() << " points: Coverage = " <<
      tree_pix.Coverage() << ", Expected Coverage = " <<
      2.0*pixVec[0].Area()/tree_pix.Area() << "\n";
  }

  // And finally, we test the version of Coverage that takes a Pixel as input.
  // First, we verify that the pixels used to add the random points
  // have unity Coverage values and another pixel not used to generate any
  // random points has Coverage equal to zero.
  std::cout << "\nFirst pixel Coverage = " << tree_pix.Coverage(pixVec[0]) <<
    " (1.0), Second pixel Coverage = " <<
    tree_pix.Coverage(pixVec[pixVec.size()-1]) << " (1.0)\n";
  std::cout << "Random other pixel Coverage = " <<
    tree_pix.Coverage(pixVec[1]) << " (0.0)\n";

  // Now we check to make sure that the super-pixels leading from the original
  // TreePixel to our current sub-pixels have the expected Coverages.
  std::cout << "\nResolution tests:\n";
  for (uint32_t resolution=tree_pix.Resolution()*2;resolution<=hi_resolution;
       resolution*=2) {
    tree_pix.SubPix(resolution, pixVec);
    std::cout << "\tResolution: " << resolution <<
      "\n\t\tFirst pixel Coverage = " << tree_pix.Coverage(pixVec[0]) <<
      " (" << 1.0*resolution*resolution/(hi_resolution*hi_resolution) << ")" <<
      "\n\t\tSecond pixel Coverage = " <<
      tree_pix.Coverage(pixVec[pixVec.size()-1]) << " (" <<
      1.0*resolution*resolution/(hi_resolution*hi_resolution) << ")\n";
    std::cout << "\t\tRandom other pixel Coverage = " <<
      tree_pix.Coverage(pixVec[1]) << " (0.0)\n";
  }
}

void IndexedTreePixelNeighborTests() {
  // Checking nearest neighbor finding routines.
  std::cout << "\n";
  std::cout << "***********************************************\n";
  std::cout << "*** IndexedTreePixel Nearest Neighbor Tests ***\n";
  std::cout << "***********************************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 100;
  uint32_t n_points = 50000;
  Stomp::IndexedTreePixel tree_pix(ang, Stomp::HPixResolution,
				   n_points_per_node);

  Stomp::AngularVector angVec;
  Stomp::IAngularVector i_angVec;
  std::cout << "Adding " << n_points << "...\n";
  tree_pix.GenerateRandomPoints(angVec, n_points);
  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    Stomp::IndexedAngularCoordinate tmp_ang(iter->Lambda(), iter->Eta(), idx,
					    Stomp::AngularCoordinate::Survey);
    i_angVec.push_back(tmp_ang);
    if (!tree_pix.AddPoint(tmp_ang))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
  }

  tree_pix.Ang(ang);

  // Quick global accounting check.
  std::cout << "\t" << tree_pix.NPoints() << " points added; " <<
    1.0*tree_pix.NPoints()/tree_pix.Area() << " points/sq. degree.\n";

  uint16_t n_test_points = 1000;
  Stomp::AngularVector test_angVec;
  tree_pix.GenerateRandomPoints(test_angVec, n_test_points);

  // Now, we check the annulus finding.
  Stomp::StompWatch stomp_watch;
  Stomp::WeightedAngularCoordinate test_point;

  std::cout <<
    "\nFinding nearest neighbor distances using " << n_test_points <<
    " points in the tree...\n";
  stomp_watch.StartTimer();
  uint16_t total_nodes = tree_pix.Nodes();
  uint16_t nodes_visited = 0;
  uint16_t failed_matches = 0;
  double mean_neighbor_distance = 0.0;
  double mean_nodes_visited = 0.0;
  Stomp::IndexedAngularCoordinate nearest_neighbor;
  for (uint16_t i=0;i<n_test_points;i++) {
    nodes_visited = tree_pix.FindNearestNeighbor(angVec[i], nearest_neighbor);
    mean_nodes_visited += static_cast<double>(nodes_visited);
    mean_neighbor_distance += angVec[i].AngularDistance(nearest_neighbor);
    if (i != nearest_neighbor.Index()) failed_matches++;
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << " (0)\n" <<
    "\t\tFailed matches: " << failed_matches << "\n" <<
    "\t\tMean nodes visited = " << mean_nodes_visited/n_test_points <<
    "/" << total_nodes << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime()/n_test_points << "s\n";


  int k = 10;
  std::cout << "\nFinding " << k << "th nearest neighbor distance using " <<
    n_test_points << " points in the tree...\n";
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  nodes_visited = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_pix.KNearestNeighborDistance(angVec[i], k, nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to " << k << "th nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";


  std::cout << "\nFinding nearest neighbor distances using " <<
    n_test_points << " points in the pixel...\n";
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  nodes_visited = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_pix.NearestNeighborDistance(test_angVec[i], nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";


  std::cout << "\nFinding " << k << "th nearest neighbor distance using " <<
    n_test_points << " points in the pixel...\n";
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  nodes_visited = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_pix.KNearestNeighborDistance(test_angVec[i], k, nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to " << k << "th nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";


  std::cout << "\nFinding nearest neighbor distances using " <<
    n_test_points << " points outside the pixel...\n";
  for (uint16_t i=0;i<n_test_points;i++)
    test_angVec[i].SetSurveyCoordinates(test_angVec[i].Lambda()-10.0,
					test_angVec[i].Eta());

  double min_edge_distance, max_edge_distance;
  double mean_edge_distance = 0.0;
  for (uint16_t i=0;i<n_test_points;i++) {
    tree_pix.EdgeDistances(test_angVec[i], min_edge_distance,
			    max_edge_distance);
    mean_edge_distance += Stomp::RadToDeg*asin(sqrt(min_edge_distance));
  }

  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  nodes_visited = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_pix.NearestNeighborDistance(test_angVec[i], nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << " (" <<
    mean_edge_distance/n_test_points << ")\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";


  std::cout << "\nFinding " << k << "th nearest neighbor distance using " <<
    n_test_points << " points outside the pixel...\n";
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  nodes_visited = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_pix.KNearestNeighborDistance(test_angVec[i], k, nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to " << k << "th nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << " (" <<
    mean_edge_distance/n_test_points << ")\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";
}

void IndexedTreePixelMatchTests() {
  // Checking closest match finding routines.
  std::cout << "\n";
  std::cout << "********************************************\n";
  std::cout << "*** IndexedTreePixel Closest Match Tests ***\n";
  std::cout << "********************************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 50;
  uint32_t n_points = 50000;
  Stomp::IndexedTreePixel tree_pix(ang, Stomp::HPixResolution,
				   n_points_per_node);

  Stomp::AngularVector angVec;
  Stomp::IAngularVector i_angVec;
  std::cout << "Adding " << n_points << "...\n";
  tree_pix.GenerateRandomPoints(angVec, n_points);
  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    Stomp::IndexedAngularCoordinate tmp_ang(iter->Lambda(), iter->Eta(), idx,
					    Stomp::AngularCoordinate::Survey);
    i_angVec.push_back(tmp_ang);
    if (!tree_pix.AddPoint(tmp_ang))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
  }

  tree_pix.Ang(ang);

  // Quick global accounting check.
  std::cout << "\t" << tree_pix.NPoints() << " points added; " <<
    1.0*tree_pix.NPoints()/tree_pix.Area() << " points/sq. degree.\n";

  uint16_t n_test_points = 1000;
  Stomp::AngularVector test_angVec;
  tree_pix.GenerateRandomPoints(test_angVec, n_test_points);

  // Now, we check the annulus finding.
  Stomp::StompWatch stomp_watch;
  Stomp::IndexedAngularCoordinate test_point;

  double match_radius = 3.0;
  std::cout <<
    "\nFinding closest matches using " << n_test_points <<
    " points in the tree and " << match_radius << " arcsecond radius...\n";
  stomp_watch.StartTimer();
  uint16_t n_match = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_pix.ClosestMatch(angVec[i], match_radius/3600.0, test_point) &&
	i == test_point.Index()) n_match++;
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points << " points.\n" <<
    "\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";

  match_radius = 0.5;
  std::cout <<
    "\nFinding closest matches using " << n_test_points <<
    " points in the tree and " << match_radius << " arcsecond radius...\n";
  stomp_watch.StartTimer();
  n_match = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_pix.ClosestMatch(angVec[i], match_radius/3600.0, test_point))
      n_match++;
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points << " points.\n" <<
    "\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";

  match_radius = 3.0;
  std::cout << "\nFinding matches using " << n_test_points <<
    " points in the pixel with " << match_radius << " arcsecond radius...\n";
  stomp_watch.StartTimer();
  n_match = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_pix.ClosestMatch(test_angVec[i], match_radius/3600.0, test_point))
      n_match++;
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points << " points.\n" <<
    "\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";

  match_radius = 0.5;
  std::cout << "\nFinding matches using " << n_test_points <<
    " points in the pixel with " << match_radius << " arcsecond radius...\n";
  stomp_watch.StartTimer();
  n_match = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_pix.ClosestMatch(test_angVec[i], match_radius/3600.0, test_point))
      n_match++;
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points << " points.\n" <<
    "\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";

  match_radius = 3.0;
  std::cout << "\nFinding matches using " << n_test_points <<
    " points outside the pixel and " << match_radius <<
    " arcsecond radius...\n";
  for (uint16_t i=0;i<n_test_points;i++)
    test_angVec[i].SetSurveyCoordinates(test_angVec[i].Lambda()-40.0,
					test_angVec[i].Eta());

  stomp_watch.StartTimer();
  n_match = 0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_pix.ClosestMatch(test_angVec[i], match_radius/3600.0, test_point))
      n_match++;
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points << " points.\n" <<
    "\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";
}

// Define our command line flags
DEFINE_bool(all_itree_pixel_tests, false, "Run all class unit tests.");
DEFINE_bool(itree_pixel_basic_tests, false, "Run IndexedTreePixel basic tests");
DEFINE_bool(itree_pixel_pair_tests, false, "Run IndexedTreePixel pair tests");
DEFINE_bool(itree_pixel_coverage_tests, false,
	    "Run IndexedTreePixel coverage tests");
DEFINE_bool(itree_pixel_neighbor_tests, false,
            "Run IndexedTreePixel nearest neighbor tests");
DEFINE_bool(itree_pixel_match_tests, false,
            "Run IndexedTreePixel closest match tests");

void IndexedTreePixelUnitTests(bool run_all_tests) {
  void IndexedTreePixelBasicTests();
  void IndexedTreePixelPairTests();
  void IndexedTreePixelCoverageTests();
  void IndexedTreePixelNeighborTests();
  void IndexedTreePixelMatchTests();

  if (run_all_tests) FLAGS_all_itree_pixel_tests = true;

  // Check that the Stomp::IndexedTreePixel class is able to add points and
  // automatically generate sub-pixels.
  if (FLAGS_all_itree_pixel_tests || FLAGS_itree_pixel_basic_tests)
    IndexedTreePixelBasicTests();

  // Check the Stomp::IndexedTreePixel pair-finding routines.
  if (FLAGS_all_itree_pixel_tests || FLAGS_itree_pixel_pair_tests)
    IndexedTreePixelPairTests();

  // Check the Stomp::IndexedTreePixel Coverage method works as advertised.
  if (FLAGS_all_itree_pixel_tests || FLAGS_itree_pixel_coverage_tests)
    IndexedTreePixelCoverageTests();

  // Checking nearest neighbor finding routines.
  if (FLAGS_all_itree_pixel_tests || FLAGS_itree_pixel_neighbor_tests)
    IndexedTreePixelNeighborTests();

  // Checking closest match finding routines.
  if (FLAGS_all_itree_pixel_tests || FLAGS_itree_pixel_match_tests)
    IndexedTreePixelMatchTests();
}
