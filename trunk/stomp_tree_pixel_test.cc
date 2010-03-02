#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_angular_bin.h"
#include "stomp_angular_correlation.h"
#include "stomp_tree_pixel.h"
#include "stomp_util.h"

void StompTreePixelTests() {
  // Before moving on to the map tests, we need to verify that the
  // Stomp::TreePixel, a derived class from Stomp::Pixel, is working properly.
  //  Let's start with some initialization.
  std::cout << "\n";
  std::cout << "******************************\n";
  std::cout << "*** Stomp::TreePixel Tests ***\n";
  std::cout << "******************************\n";
  Stomp::WeightedAngularCoordinate* ang =
    new Stomp::WeightedAngularCoordinate(60.0, 0.0, 1.0,
					 Stomp::AngularCoordinate::Survey);
  uint16_t n_points = 50;
  Stomp::TreePixel tree_pix(*ang, Stomp::HPixResolution, n_points);

  if (tree_pix.Contains(*ang)) {
    std::cout << "\tSuccessfully verified the position used to set pixnum.\n";
  } else {
    std::cout << "\tFailed to verify the position used to set pixnum.\n";
  }

  if (tree_pix.AddPoint(ang)) {
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
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if (!tree_pix.AddPoint(*iter)) {
      std::cout << "\t\tFailed to add point: " <<
	iter->Lambda() << ", " << iter->Eta() << "\n";
      exit(1);
    }
  }
  std::cout << "\tAdded " << tree_pix.NPoints() << "/" << 200*n_points <<
    " points to this pixel\n";
}

void StompTreePixelPairTests() {
  std::cout << "\n";
  std::cout << "*********************************\n";
  std::cout << "*** StompTreePixel Pair Tests ***\n";
  std::cout << "*********************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 500;
  uint32_t n_points = 50000;
  Stomp::TreePixel tree_pix(ang, Stomp::HPixResolution,
			    n_points_per_node);

  Stomp::AngularVector angVec;
  tree_pix.GenerateRandomPoints(angVec, n_points);
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if (!tree_pix.AddPoint(*iter))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }

  tree_pix.Ang(ang);
  // First, we choose a radius that is certain to encompass the entire pixel.
  double theta = 15.0;
  uint32_t n_pairs = tree_pix.FindPairs(ang, theta);
  std::cout << "Found " << n_pairs << " pairs (" << n_points << ")\n";

  // Next, we choose a smaller radius where we should begin to probe some of
  // the tree structure in the pixel.
  theta = 0.1;
  uint32_t direct_pairs = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if (angVec[0].AngularDistance(*iter) <= theta) direct_pairs++;
  }
  std::cout << "Found " << tree_pix.FindPairs(angVec[0], theta) <<
    " pairs; " << direct_pairs << " pairs from brute force.\n";

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
  std::cout << "Found " <<
    tree_pix.FindPairs(angVec[0], theta_min, theta_max) <<
    " pairs; "<< direct_pairs << " pairs from brute force.\n";

  Stomp::AngularBin theta_bin(theta_min, theta_max);

  uint32_t mismatched_pairs = 0;
  uint32_t overcounted_pairs = 0;
  uint32_t undercounted_pairs = 0;
  uint32_t counter_offset_pairs = 0;
  for (Stomp::AngularIterator ang_iter=angVec.begin();
       ang_iter!=angVec.end();++ang_iter) {
    n_pairs = tree_pix.FindPairs(*ang_iter, theta_bin);
    direct_pairs = 0;
    for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
      if (theta_bin.WithinCosBounds(ang_iter->DotProduct(*iter)))
	direct_pairs++;
    }
    if (n_pairs != direct_pairs || n_pairs != theta_bin.Counter()) {
      mismatched_pairs++;
      if (n_pairs > direct_pairs) overcounted_pairs++;
      if (n_pairs < direct_pairs) undercounted_pairs++;
      if (n_pairs != theta_bin.Counter()) counter_offset_pairs++;
      if (mismatched_pairs == 1) ang = *ang_iter;
    }
    theta_bin.ResetCounter();
  }
  std::cout << "Found discrepancy between counts for " << mismatched_pairs <<
    " of " << angVec.size() << " points (" << overcounted_pairs <<
    " over, " << undercounted_pairs << " under, " << counter_offset_pairs <<
    " counter offset).\n";

  std::cout << "Angular Bin full pairs test:\n";
  Stomp::AngularCorrelation wtheta(0.01, 5.0, 5.0, false);
  for (Stomp::ThetaIterator iter=wtheta.Begin();iter!=wtheta.End();++iter)
    iter->ResetCounter();

  std::cout << "\tCalculating with FindPairs...\n";
  tree_pix.FindWeightedPairs(angVec, wtheta);

  std::cout << "\tDone.  Starting brute force calculation...\n";
  Stomp::AngularCorrelation wtheta_brute(0.01, 5.0, 5.0, false);
  Stomp::ThetaIterator theta_begin = wtheta_brute.Begin();
  Stomp::ThetaIterator theta_end = wtheta_brute.End();
  Stomp::ThetaIterator theta_iter;
  double costheta = 0.0;
  for (Stomp::AngularIterator ang_iter=angVec.begin();
       ang_iter!=angVec.end();++ang_iter) {
    for (Stomp::AngularIterator inner_iter=angVec.begin();
	 inner_iter!=angVec.end();++inner_iter) {
      costheta = ang_iter->DotProduct(*inner_iter);
      theta_iter = wtheta_brute.Find(theta_begin, theta_end,
				     1.0-costheta*costheta);
      if (theta_iter != theta_end) theta_iter->AddToCounter();
    }
  }

  for (Stomp::ThetaIterator iter=wtheta.Begin();iter!=wtheta.End();++iter) {
    double sin2theta = 0.5*(iter->Sin2ThetaMin()+iter->Sin2ThetaMax());
    theta_iter = wtheta_brute.Find(theta_begin, theta_end, sin2theta);
    std::cout << "\t" << iter->ThetaMin() << " - " << iter->ThetaMax() <<
      ": " << iter->Counter() << " pairs; " << theta_iter->Counter() <<
      " brute force pairs.\n";
  }
}

void StompTreePixelCoverageTests() {
  // Since we don't know a priori what the underlying distribution of the
  // points added to our TreePixel is, we estimate it using the Coverage
  // method.  This module verifies that Coverage is working properly.
  std::cout << "\n";
  std::cout << "***************************************\n";
  std::cout << "*** Stomp::TreePixel Coverage Tests ***\n";
  std::cout << "***************************************\n";
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points = 20;
  Stomp::TreePixel tree_pix(ang, Stomp::HPixResolution, n_points);

  // Now we find the sub-pixels for the current tree_pix a few levels down
  // and start generating random points from them.
  Stomp::PixelVector pixVec;
  uint16_t hi_resolution = Stomp::HPixResolution*32;
  tree_pix.SubPix(hi_resolution, pixVec);

  Stomp::AngularVector angVec;
  pixVec[0].GenerateRandomPoints(angVec, n_points);
  bool added_point = false;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    added_point = tree_pix.AddPoint(*iter);
    if (!added_point)
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
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
      added_point = tree_pix.AddPoint(*iter);
      if (!added_point)
	std::cout << "\t\tFailed to add point: " <<
	  iter->RA() << ", " << iter->DEC() << "\n";
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
      added_point = tree_pix.AddPoint(*iter);
      if (!added_point)
	std::cout << "\t\tFailed to add point: " <<
	  iter->RA() << ", " << iter->DEC() << "\n";
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
  for (uint16_t resolution=tree_pix.Resolution()*2;resolution<=hi_resolution;
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

void StompTreePixelFieldPairTests() {
  // Checking pair finding routines with Field values.
  std::cout << "\n";
  std::cout << "***************************************\n";
  std::cout << "*** StompTreePixel Field Pair Tests ***\n";
  std::cout << "***************************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 500;
  uint32_t n_points = 50000;
  Stomp::TreePixel tree_pix(ang, Stomp::HPixResolution,
			    n_points_per_node);

  Stomp::AngularVector angVec;
  Stomp::WAngularVector w_angVec;
  std::cout << "Adding " << n_points <<
    " points with Weight() = 1 and Field('two') = 2\n";
  tree_pix.GenerateRandomPoints(angVec, n_points);
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    Stomp::WeightedAngularCoordinate tmp_ang(iter->Lambda(), iter->Eta(), 1.0,
					     Stomp::AngularCoordinate::Survey);
    tmp_ang.SetField("two", 2.0);
    w_angVec.push_back(tmp_ang);
    if (!tree_pix.AddPoint(tmp_ang))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }

  tree_pix.Ang(ang);

  // Quick global accounting check.
  std::cout << "\t" << tree_pix.NPoints() <<
    " points added.\n\tTotal Weight = " << tree_pix.Weight() <<
    "\n\tTotal Field('two') = " << tree_pix.FieldTotal("two") <<
    "\n\tTotal Field('three') = " << tree_pix.FieldTotal("three") << "\n";

  // Now, we check the annulus finding.
  Stomp::StompWatch stomp_watch;
  Stomp::AngularBin theta(0.05, 0.15);

  std::cout << "\nCalculating with FindPairs for Weight()...\n";
  stomp_watch.StartTimer();
  tree_pix.FindWeightedPairs(angVec, theta);
  stomp_watch.StopTimer();

  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", total Weight() = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  theta.Reset();

  std::cout << "\nCalculating with FindPairs for Field('two')...\n";
  stomp_watch.StartTimer();
  tree_pix.FindWeightedPairs(angVec, theta, "two");
  stomp_watch.StopTimer();

  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", FieldTotal('two') = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  theta.Reset();

  std::cout << "\nCalculating with FindPairs for Field('three')...\n";
  stomp_watch.StartTimer();
  tree_pix.FindWeightedPairs(angVec, theta, "three");
  stomp_watch.StopTimer();

  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", FieldTotal('three') = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  theta.Reset();

  std::cout <<
    "\nCalculating with FindPairs for Weight() x Field('two')...\n";
  stomp_watch.StartTimer();
  tree_pix.FindWeightedPairs(w_angVec, theta, "two");
  stomp_watch.StopTimer();

  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", FieldTotal('two') = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  theta.Reset();

  std::cout <<
    "\nCalculating with FindPairs for Field('two') x Field('two')...\n";
  stomp_watch.StartTimer();
  tree_pix.FindWeightedPairs(w_angVec, "two", theta, "two");
  stomp_watch.StopTimer();

  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", FieldTotal('two') = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  theta.Reset();

  std::cout <<
    "\nCalculating with FindPairs for Field('two') x Field('three')...\n";
  stomp_watch.StartTimer();
  tree_pix.FindWeightedPairs(w_angVec, "two", theta, "three");
  stomp_watch.StopTimer();

  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", FieldTotal('three') = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";
}

void StompTreePixelNeighborTests() {
  // Checking nearest neighbor finding routines.
  std::cout << "\n";
  std::cout << "*********************************************\n";
  std::cout << "*** StompTreePixel Nearest Neighbor Tests ***\n";
  std::cout << "*********************************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 100;
  uint32_t n_points = 50000;
  Stomp::TreePixel tree_pix(ang, Stomp::HPixResolution,
			    n_points_per_node);

  Stomp::AngularVector angVec;
  Stomp::WAngularVector w_angVec;
  std::cout << "Adding " << n_points << "...\n";
  tree_pix.GenerateRandomPoints(angVec, n_points);
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    Stomp::WeightedAngularCoordinate tmp_ang(iter->Lambda(), iter->Eta(), 1.0,
					     Stomp::AngularCoordinate::Survey);
    w_angVec.push_back(tmp_ang);
    if (!tree_pix.AddPoint(*iter))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
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
  double mean_neighbor_distance = 0.0;
  double mean_nodes_visited = 0.0;
  for (uint16_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_pix.NearestNeighborDistance(angVec[i], nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << " (0)\n" <<
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
    tree_pix._EdgeDistances(test_angVec[i], min_edge_distance,
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

DEFINE_bool(all_tree_pixel_tests, false, "Run all class unit tests.");
DEFINE_bool(tree_pixel_tests, false, "Run TreePixel tests");
DEFINE_bool(tree_pixel_pair_tests, false, "Run TreePixel pair tests");
DEFINE_bool(tree_pixel_coverage_tests, false, "Run TreePixel coverage tests");
DEFINE_bool(tree_pixel_field_pair_tests, false,
            "Run TreePixel Field pair tests");
DEFINE_bool(tree_pixel_neighbor_tests, false,
            "Run TreePixel nearest neighbor tests");

int main(int argc, char **argv) {
  void StompTreePixelTests();
  void StompTreePixelPairTests();
  void StompTreePixelCoverageTests();
  void StompTreePixelFieldPairTests();
  void StompTreePixelNeighborTests();

  std::string usage = "Usage: ";
  usage += argv[0];
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Check that the Stomp::TreePixel class is able to add points and
  // automatically generate sub-pixels.
  if (FLAGS_all_tree_pixel_tests || FLAGS_tree_pixel_tests)
    StompTreePixelTests();

  // Check the Stomp::TreePixel pair-finding routines.
  if (FLAGS_all_tree_pixel_tests || FLAGS_tree_pixel_pair_tests)
    StompTreePixelPairTests();

  // Check the Stomp::TreePixel Coverage method works as advertised.
  if (FLAGS_all_tree_pixel_tests || FLAGS_tree_pixel_coverage_tests)
    StompTreePixelCoverageTests();

  // Checking pair finding routines with Field values.
  if (FLAGS_all_tree_pixel_tests || FLAGS_tree_pixel_field_pair_tests)
    StompTreePixelFieldPairTests();

  // Checking nearest neighbor finding routines.
  if (FLAGS_all_tree_pixel_tests || FLAGS_tree_pixel_neighbor_tests)
    StompTreePixelNeighborTests();

  return 0;
}
