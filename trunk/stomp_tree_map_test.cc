#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_util.h"
#include "stomp_angular_coordinate.h"
#include "stomp_angular_bin.h"
#include "stomp_angular_correlation.h"
#include "stomp_pixel.h"
#include "stomp_map.h"
#include "stomp_tree_map.h"
#include "stomp_tree_map_test.h"

void TreeMapBasicTests() {
  std::cout << "\n";
  std::cout << "***************************\n";
  std::cout << "*** TreeMap Basic Tests ***\n";
  std::cout << "***************************\n";
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a non-standard resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint16_t resolution = 32;
  Stomp::TreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::TreeMap at " << resolution <<
    " resolution...\n";

  if (tree_map.AddPoint(ang)) {
    std::cout << "\tSuccessfully added a random position to the map.\n";
  } else {
    std::cout << "\tFailed to add a random position to the map.\n";
    exit(1);
  }

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::TreeMap.  We choose a large enough radius
  // to make sure that we're adding points that are well outside of pixel that
  // contained the initial point to make sure that we're properly adding new
  // pixels as we ingest points.
  double theta = 5.0;
  Stomp::Pixel tmp_pix(ang, resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  uint32_t n_points = 500000;
  Stomp::AngularVector angVec;
  stomp_map->GenerateRandomPoints(angVec, n_points);

  Stomp::StompWatch stomp_watch;
  bool added_point = false;
  stomp_watch.StartTimer();
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    added_point = tree_map.AddPoint(*iter);
    if (!added_point)
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }
  stomp_watch.StopTimer();
  std::cout << "\tAdded " << tree_map.NPoints() << "/" << n_points <<
    " points to this map (" << n_points/stomp_watch.ElapsedTime() <<
    " points/second)\n";

  std::cout << "\t" << tree_map.BaseNodes() << " nodes at " <<
    tree_map.Resolution() << " resolution used for this map; " <<
    tree_map.Nodes() << " total nodes.\n";

  // First we check to make sure that we're recovering the base nodes properly.
  // The Coverage method will raise alarms if this fails.  Otherwise, this can
  // pass without comment.
  Stomp::PixelVector superpix;
  tree_map.Coverage(superpix, tree_map.Resolution());

  // Now we verify that the superpixels for our TreeMap match those for the
  // Map that we built it from.
  tree_map.Coverage(superpix);
  tree_map.Coverage(superpix, tree_map.Resolution());
  tree_map.Coverage(superpix);
  std::cout << "\tMap covers " << superpix.size() << " superpixels:\n\t\t";
  for (Stomp::PixelIterator iter=superpix.begin();iter!=superpix.end();++iter)
    std::cout << iter->Superpixnum() << " ";
  std::cout << "\n";
  stomp_map->Coverage(superpix);
  std::cout << "\tInput map covers " << superpix.size() <<
    " superpixels:\n\t\t";
  for (Stomp::PixelIterator iter=superpix.begin();iter!=superpix.end();++iter)
    std::cout << iter->Superpixnum() << " ";
  std::cout << "\n";
}

void TreeMapPairTests() {
  std::cout << "\n";
  std::cout << "**************************\n";
  std::cout << "*** TreeMap Pair Tests ***\n";
  std::cout << "**************************\n";
  uint16_t n_points_per_node = 200;
  uint16_t resolution = 8;
  Stomp::TreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::TreeMap at " << resolution <<
    " resolution with " << n_points_per_node << " points per node...\n";
  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::TreeMap.  We choose a large enough radius
  // to make sure that we're adding points that are well outside of pixel that
  // contained the initial point to make sure that we're properly adding new
  // pixels as we ingest points.
  double lambda = 60.0;
  double eta = 0.0;
  Stomp::AngularCoordinate ang(lambda, eta, Stomp::AngularCoordinate::Survey);
  double theta_radius = 5.0;
  Stomp::Pixel tmp_pix(ang, 32);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta_radius, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);
  uint32_t n_points = 10000;
  Stomp::AngularVector angVec;
  stomp_map->GenerateRandomPoints(angVec, n_points);

  Stomp::PixelVector base_nodes;
  stomp_map->Coverage(base_nodes, resolution);
  std::cout << "\tInput map covers " << base_nodes.size() <<
    " basenodes:\n\t\t";
  for (Stomp::PixelIterator iter=base_nodes.begin();
       iter!=base_nodes.end();++iter)
    std::cout << iter->Superpixnum() << " ";
  std::cout << "\n";

  bool added_point = false;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    added_point = tree_map.AddPoint(*iter);
    if (!added_point)
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }
  std::cout << "\t" << tree_map.BaseNodes() << " base nodes at " <<
    tree_map.Resolution() << " resolution; " << tree_map.Nodes() <<
    " total nodes.\n";

  // We start by choosing a radius much larger than our area, which should mean
  // that we find all points in the map as pairs.
  std::cout << "\nEncompassing radius test:\n";
  uint32_t n_pairs = tree_map.FindPairs(ang, 4.0*theta_radius);
  std::cout << "\tFound " << n_pairs << " pairs (expect " << n_points << ")\n";

  // Next, we choose a smaller radius where we should begin to probe some of
  // the tree structure in the pixel.
  double theta = 0.1;
  uint32_t direct_pairs = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if (ang.AngularDistance(*iter) <= theta) direct_pairs++;
  }
  std::cout << "\nSmall circle test:\n";
  std::cout << "\tFound " << tree_map.FindPairs(ang, theta) << " pairs; " <<
    direct_pairs << " pairs from brute force.\n";

  // Now, we check the annulus finding.
  double theta_max = 0.5;
  double theta_min = 0.3;
  direct_pairs = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if ((ang.AngularDistance(*iter) <= theta_max) &&
	(ang.AngularDistance(*iter) >= theta_min)) direct_pairs++;
  }
  std::cout << "\nAnnulus test:\n";
  std::cout << "\tFound " << tree_map.FindPairs(ang, theta_min, theta_max) <<
    " pairs; " << direct_pairs << " pairs from brute force.\n";

  // Now we test to see if the pair finding works at the edges of our map.
  // First, the top of the map, where the angular distortion should be
  // maximum.
  double scale_factor = 1.0;
  ang.SetSurveyCoordinates(lambda+scale_factor*theta_radius,eta);
  while (stomp_map->FindLocation(ang)) {
    scale_factor += 0.1;
    ang.SetSurveyCoordinates(lambda+scale_factor*theta_radius,eta);
  }

  direct_pairs = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if ((ang.AngularDistance(*iter) <= theta_max) &&
	(ang.AngularDistance(*iter) >= theta_min)) direct_pairs++;
  }
  std::cout << "\nTop of the circle test (" <<
    ang.Lambda() << ", " << ang.Eta() << "):\n";
  std::cout << "\tFound " << tree_map.FindPairs(ang, theta_min, theta_max) <<
    " pairs; "<< direct_pairs << " pairs from brute force.\n";

  // Now, a location close to the right edge of the circle.
  scale_factor = 1.0;
  ang.SetSurveyCoordinates(lambda,eta+scale_factor*theta_radius);
  while (stomp_map->FindLocation(ang)) {
    scale_factor += 0.1;
    ang.SetSurveyCoordinates(lambda,eta+scale_factor*theta_radius);
  }

  direct_pairs = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if ((ang.AngularDistance(*iter) <= theta_max) &&
	(ang.AngularDistance(*iter) >= theta_min)) direct_pairs++;
  }
  std::cout << "\nSide of the circle test (" <<
    ang.Lambda() << ", " << ang.Eta() << "):\n";
  std::cout << "\tFound " << tree_map.FindPairs(ang, theta_min, theta_max) <<
    " pairs; "<< direct_pairs << " pairs from brute force.\n";

  Stomp::StompWatch stomp_watch;
  std::cout << "\nAngular Bin test:\n";
  // We initialize the AngularCorrelation object with the last parameter set
  // to "false".  This will prevent the constructor from assigning resolutions
  // to each of the angular bins, marking them all as bins to use the
  // pair-counting estimator on.
  Stomp::AngularCorrelation wtheta(0.01, 1.0, 5.0, false);
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(angVec, wtheta);
  stomp_watch.StopTimer();
  std::cout << "\tTime elapsed: " << stomp_watch.ElapsedTime() << " seconds.\n";

  std::cout << "Checking direct counts...\n";
  Stomp::AngularCorrelation wtheta_brute(0.01, 1.0, 5.0, false);
  Stomp::ThetaIterator theta_begin = wtheta_brute.Begin();
  Stomp::ThetaIterator theta_end = wtheta_brute.End();
  Stomp::ThetaIterator theta_iter;
  double costheta = 0.0;
  stomp_watch.StartTimer();
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
  stomp_watch.StopTimer();
  std::cout << "\tTime elapsed: " << stomp_watch.ElapsedTime() << " seconds.\n";

  std::cout << "Direct counting comparison:\n";
  for (Stomp::ThetaIterator iter=wtheta.Begin();iter!=wtheta.End();++iter) {
    double sin2theta = 0.5*(iter->Sin2ThetaMin()+iter->Sin2ThetaMax());
    theta_iter = wtheta_brute.Find(theta_begin, theta_end, sin2theta);
    std::cout << "\t" << iter->ThetaMin() << " - " << iter->ThetaMax() <<
      ": " << iter->Counter() << " pairs; " << theta_iter->Counter() <<
      " brute force pairs.\n";
  }
}

void TreeMapAreaTests() {
  std::cout << "\n";
  std::cout << "**************************\n";
  std::cout << "*** TreeMap Area Tests ***\n";
  std::cout << "**************************\n";
  // The goal here is to check that the estimate of the area associated with
  // the TreeMap improves as we add points.

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a coarse resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint16_t resolution = 4;
  Stomp::TreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::TreeMap at " << resolution <<
    " resolution...\n";

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::TreeMap.  We choose a large enough radius
  // to make sure that we're adding points that are well outside of pixel that
  // contained the initial point to make sure that we're properly adding new
  // pixels as we ingest points.
  double theta = 5.0;
  uint16_t annulus_resolution = 32;
  Stomp::Pixel tmp_pix(ang, annulus_resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Now we add an increasing number of points and track how the area of our
  // TreeMap compares with the area of the Map used to generate the random
  // points.  Since we used a finer resolution to make the Map than for our
  // TreeMap, we should start with a larger TreeMap area and then eventually
  // get closer to the Map area as we add points.
  for (uint32_t k=1;k<200;k*=2) {
    uint32_t n_points = 5000*k;

    Stomp::AngularVector angVec;
    stomp_map->GenerateRandomPoints(angVec, n_points);

    bool added_point = false;
    for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
      added_point = tree_map.AddPoint(*iter);
      if (!added_point)
	std::cout << "\t\tFailed to add point: " <<
	  iter->RA() << ", " << iter->DEC() << "\n";
    }

    std::cout << "\t" << tree_map.NPoints() << " points: " <<
      tree_map.BaseNodes() << " nodes at " << tree_map.Resolution() <<
      " resolution; " << tree_map.Nodes() << " total nodes.\n";
    std::cout << "\t\tInput map area: " << stomp_map->Area() <<
      "; TreeMap area: " << tree_map.Area() << "\n";
  }
}

void TreeMapRegionTests() {
  std::cout << "\n";
  std::cout << "****************************\n";
  std::cout << "*** TreeMap Region Tests ***\n";
  std::cout << "****************************\n";
  // The goal here is to check that the regionation code works, provided that
  // we have a dense sampling of our data area.

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a coarse resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint16_t resolution = 4;
  Stomp::TreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::TreeMap at " << resolution <<
    " resolution...\n";

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::TreeMap.  We choose a large enough radius
  // to make sure that we're adding points that are well outside of pixel that
  // contained the initial point to make sure that we're properly adding new
  // pixels as we ingest points.
  double theta = 5.0;
  uint16_t annulus_resolution = 32;
  Stomp::Pixel tmp_pix(ang, annulus_resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Now we add a high number of points to the TreeMap to flesh out its bounds.
  uint32_t n_points = 1000000;
  Stomp::AngularVector angVec;
  stomp_map->GenerateRandomPoints(angVec, n_points);

  bool added_point = false;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    added_point = tree_map.AddPoint(*iter);
    if (!added_point)
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }
  angVec.clear();

  std::cout << "\t" << tree_map.NPoints() << " points: " <<
    tree_map.Nodes() << " nodes at " << tree_map.Resolution() <<
    " resolution used for this map\n";
  std::cout << "\t\tInput map area: " << stomp_map->Area() <<
    "; TreeMap area: " << tree_map.Area() << "\n";

  // Before regionating our TreeMap, let's verify that the unmasked fraction
  // calculations that go into the Coverage maps are calculated correctly.
  Stomp::PixelVector map_coverage;
  stomp_map->Coverage(map_coverage);

  Stomp::PixelVector tree_coverage;
  tree_map.Coverage(tree_coverage);

  if (map_coverage.size() != tree_coverage.size()) {
    std::cout << "Disagreement between Map Coverage (" <<
      map_coverage.size() << " pixels) and TreeMap Coverage (" <<
      tree_coverage.size() << " pixels).\n";
  } else {
    std::cout << "Unmasked fractions:\n";
    for (uint32_t i=0;i<map_coverage.size();i++)
      std::cout << "\t" << map_coverage[i].Pixnum() << ": Map = " <<
	map_coverage[i].Weight() << ", TreeMap = " <<
	tree_coverage[i].Weight() << "\n";
  }

  // Now we regionate both maps into equal numbers of regions and compare
  // the areas in each of the regions.
  uint16_t n_regions = 10;
  uint16_t n_map_regions = stomp_map->InitializeRegions(n_regions, 4);
  uint16_t n_tree_regions = tree_map.InitializeRegions(n_regions);

  if (n_map_regions != n_tree_regions) {
    std::cout << "\tAttempted to split both maps into " << n_regions <<
      " regions;\n\t\tGot " << n_map_regions << " for the Map and " <<
      n_tree_regions << " for the TreeMap.\n";
  } else {
    std::cout << "\tGot " << n_map_regions << " regions for both maps:\n";
    for (uint16_t i=0;i<n_map_regions;i++) {
      std::cout << "\t\tRegion " << i << ": Map area = " <<
	stomp_map->RegionArea(i) << ", TreeMap area = " <<
	tree_map.RegionArea(i) << "\n";
    }
  }
}

void TreeMapFieldPairTests() {
  // Checking pair finding routines with Field values.
  std::cout << "\n";
  std::cout << "********************************\n";
  std::cout << "*** TreeMap Field Pair Tests ***\n";
  std::cout << "********************************\n";

  uint16_t resolution = 8;
  uint16_t n_points_per_node = 200;
  Stomp::TreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::TreeMap at " << resolution <<
    " resolution with " << n_points_per_node << " points per node...\n";

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::TreeMap.  We choose a large enough radius
  // to make sure that we're adding points that are well outside of pixel that
  // contained the initial point to make sure that we're properly adding new
  // pixels as we ingest points.
  double theta_bound = 5.0;
  uint16_t annulus_resolution = 32;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, annulus_resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta_bound, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  Stomp::PixelVector base_nodes;
  stomp_map->Coverage(base_nodes, resolution);
  std::cout << "\tInput map covers " << base_nodes.size() <<
    " basenodes:\n\t\t";
  for (Stomp::PixelIterator iter=base_nodes.begin();
       iter!=base_nodes.end();++iter)
    std::cout << iter->Superpixnum() << " ";
  std::cout << "\n";

  // Now we add a high number of points to the TreeMap to flesh out its bounds.
  uint32_t n_points = 100000;
  Stomp::AngularVector angVec;
  Stomp::WAngularVector w_angVec;
  std::cout << "Adding " << n_points <<
    " points with Weight() = 1 and Field('two') = 2\n";
  stomp_map->GenerateRandomPoints(angVec, n_points);

  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    Stomp::WeightedAngularCoordinate tmp_ang(iter->Lambda(), iter->Eta(), 1.0,
					     Stomp::AngularCoordinate::Survey);
    tmp_ang.SetField("two", 2.0);
    w_angVec.push_back(tmp_ang);

    if (!tree_map.AddPoint(tmp_ang))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }

  // Quick global accounting check.
  std::cout << "\t" << tree_map.BaseNodes() << " base nodes at " <<
    tree_map.Resolution() << " resolution; " << tree_map.Nodes() <<
    " total nodes.\n";
  std::cout << "\t" << tree_map.NPoints() <<
    " points added.\n\tTotal Weight = " << tree_map.Weight() <<
    "\n\tTotal Field('two') = " << tree_map.FieldTotal("two") <<
    "\n\tTotal Field('three') = " << tree_map.FieldTotal("three") << "\n";

  // Now, we check the annulus finding.
  Stomp::StompWatch stomp_watch;
  Stomp::AngularBin theta(0.05, 0.15);

  std::cout << "\nCalculating with FindPairs for Weight()...\n";
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(angVec, theta);
  stomp_watch.StopTimer();
  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", total Weight() = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  std::cout << "\nCalculating with FindPairs for Field('two')...\n";
  theta.Reset();
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(angVec, theta, "two");
  stomp_watch.StopTimer();
  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", FieldTotal('two') = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  std::cout << "\nCalculating with FindPairs for Field('three')...\n";
  theta.Reset();
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(angVec, theta, "three");
  stomp_watch.StopTimer();
  std::cout << "\tTotal pairs = " << theta.Counter() <<
    ", FieldTotal('three') = " << theta.Weight() << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  std::cout <<
    "\nCalculating with FindPairs for Weight() x Field('two')...\n";
  theta.Reset();
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(w_angVec, theta, "two");
  stomp_watch.StopTimer();
  std::cout << "\tTotal pairs = " << theta.Counter() <<
    "\n\tWeight()xFieldTotal('two') = " << theta.Weight() <<
    "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime() << "s\n";

  std::cout <<
    "\nCalculating with FindPairs for Field('two') x Field('two')...\n";
  theta.Reset();
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(w_angVec, "two", theta, "two");
  stomp_watch.StopTimer();
  std::cout << "\tTotal pairs = " << theta.Counter() <<
    "\n\tFieldTotal('two')xFieldTotal('two') = " << theta.Weight() <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime() << "s\n";

  std::cout <<
    "\nCalculating with FindPairs for Field('two') x Field('three')...\n";
  theta.Reset();
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(w_angVec, "two", theta, "three");
  stomp_watch.StopTimer();
  std::cout << "\tTotal pairs = " << theta.Counter() <<
    "\n\tFieldTotal('two')xFieldTotal('three') = " << theta.Weight() <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime() << "s\n";

  // One last test to see how the access time scales with the number of fields
  tree_map.Clear();
  std::cout << "\nAdding Field('three')...\n";
  for (Stomp::WAngularIterator iter=w_angVec.begin();
       iter!=w_angVec.end();++iter) {
    iter->SetField("three", 3.0);

    if (!tree_map.AddPoint(*iter))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }

  std::cout <<
    "\tCalculating with FindPairs for Field('two') x Field('three')...\n";
  theta.Reset();
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(w_angVec, "two", theta, "three");
  stomp_watch.StopTimer();
  std::cout << "\t\tTotal pairs = " << theta.Counter() <<
    "\n\t\tFieldTotal('two')xFieldTotal('three') = " << theta.Weight() <<
    "\n\t\t\tTime elapsed = " << stomp_watch.ElapsedTime() << "s\n";

  tree_map.Clear();
    std::cout << "\nAdding Field('four'), Field('five') and Field('six')...\n";
  for (Stomp::WAngularIterator iter=w_angVec.begin();
       iter!=w_angVec.end();++iter) {
    iter->SetField("four", 4.0);
    iter->SetField("five", 5.0);
    iter->SetField("six", 6.0);

    if (!tree_map.AddPoint(*iter))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }

  std::cout <<
    "\tCalculating with FindPairs for Field('two') x Field('three')...\n";
  theta.Reset();
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(w_angVec, "two", theta, "three");
  stomp_watch.StopTimer();
  std::cout << "\t\tTotal pairs = " << theta.Counter() <<
    "\n\t\tFieldTotal('two')xFieldTotal('three') = " << theta.Weight() <<
    "\n\t\t\tTime elapsed = " << stomp_watch.ElapsedTime() << "s\n";

  tree_map.Clear();
    std::cout <<
      "\nAdding Field('seven'), Field('eight'), Field('nine') & Field('ten')\n";
  for (Stomp::WAngularIterator iter=w_angVec.begin();
       iter!=w_angVec.end();++iter) {
    iter->SetField("seven", 7.0);
    iter->SetField("eight", 8.0);
    iter->SetField("nine", 9.0);
    iter->SetField("ten", 10.0);

    if (!tree_map.AddPoint(*iter))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }

  std::cout <<
    "\tCalculating with FindPairs for Field('two') x Field('ten')...\n";
  theta.Reset();
  stomp_watch.StartTimer();
  tree_map.FindWeightedPairs(w_angVec, "two", theta, "ten");
  stomp_watch.StopTimer();
  std::cout << "\t\tTotal pairs = " << theta.Counter() <<
    "\n\t\tFieldTotal('two')xFieldTotal('ten') = " << theta.Weight() <<
    "\n\t\t\tTime elapsed = " << stomp_watch.ElapsedTime() << "s\n";
}

void TreeMapNeighborTests() {
  // Checking nearest neighbor finding routines.
  std::cout << "\n";
  std::cout << "**************************************\n";
  std::cout << "*** TreeMap Nearest Neighbor Tests ***\n";
  std::cout << "**************************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a coarse resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint16_t resolution = 4;
  Stomp::TreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::TreeMap at " << resolution <<
    " resolution...\n";

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::TreeMap.  We choose a large enough radius
  // to make sure that we're adding points that are well outside of pixel that
  // contained the initial point to make sure that we're properly adding new
  // pixels as we ingest points.
  double theta_bound = 5.0;
  uint16_t annulus_resolution = 32;
  Stomp::Pixel tmp_pix(ang, annulus_resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta_bound, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Now we add a high number of points to the TreeMap to flesh out its bounds.
  uint32_t n_points = 100000;
  Stomp::AngularVector angVec;
  Stomp::WAngularVector w_angVec;
  std::cout << "Adding " << n_points << " points\n";
  stomp_map->GenerateRandomPoints(angVec, n_points);

  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    Stomp::WeightedAngularCoordinate tmp_ang(iter->Lambda(), iter->Eta(), 1.0,
					     Stomp::AngularCoordinate::Survey);
    if (!tree_map.AddPoint(tmp_ang))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
  }

  // Quick global accounting check.
  std::cout << "\t" << tree_map.NPoints() << " points added; " <<
    1.0*tree_map.NPoints()/tree_map.Area() << " points/sq. degree.\n";

  // Generate a smaller set of test points.
  uint32_t n_test_points = 1000;
  Stomp::AngularVector test_angVec;
  stomp_map->GenerateRandomPoints(test_angVec, n_test_points);

  // Now, we check the annulus finding.
  Stomp::StompWatch stomp_watch;
  Stomp::WeightedAngularCoordinate test_point;

  std::cout <<
    "\nFinding nearest neighbor distances using " << n_test_points <<
    " points in the tree...\n";
  uint16_t total_nodes = tree_map.Nodes();
  uint16_t nodes_visited = 0;
  double mean_neighbor_distance = 0.0;
  double mean_nodes_visited = 0.0;
  stomp_watch.StartTimer();
  for (uint32_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_map.NearestNeighborDistance(angVec[i], nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";

  int k = 10;
  std::cout << "\nFinding " << k <<
    "th nearest neighbor distance using " << n_test_points <<
    " points in the tree...\n";
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  mean_nodes_visited = 0.0;
  for (uint32_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_map.KNearestNeighborDistance(angVec[i], k, nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to " << k << "th nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";


  std::cout <<
    "\nFinding nearest neighbor distances using " << n_test_points <<
    " points in the map...\n";
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  mean_nodes_visited = 0.0;
  for (uint32_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_map.NearestNeighborDistance(test_angVec[i], nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";

  std::cout << "\nFinding " << k << "th nearest neighbor distance using " <<
    n_test_points << " points in the map...\n";
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  mean_nodes_visited = 0.0;
  for (uint32_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_map.KNearestNeighborDistance(test_angVec[i], k, nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to " << k << "th nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";


  std::cout <<
    "\nFinding nearest neighbor distances using " << n_test_points <<
    " points outside the map...\n";
  for (uint32_t i=0;i<n_test_points;i++)
    test_angVec[i].SetSurveyCoordinates(test_angVec[i].Lambda()-3.0*theta_bound,
					0.0);
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  mean_nodes_visited = 0.0;
  for (uint32_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_map.NearestNeighborDistance(test_angVec[i], nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";

  std::cout << "\nFinding " << k << "th nearest neighbor distance using " <<
    n_test_points << " points outside the map...\n";
  stomp_watch.StartTimer();
  mean_neighbor_distance = 0.0;
  mean_nodes_visited = 0.0;
  for (uint32_t i=0;i<n_test_points;i++) {
    mean_neighbor_distance +=
      tree_map.KNearestNeighborDistance(test_angVec[i], k, nodes_visited);
    mean_nodes_visited += static_cast<double>(nodes_visited);
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to " << k << "th nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n\t\tMean nodes visited = " <<
    mean_nodes_visited/n_test_points << "/" << total_nodes <<
    "\n\t\tTime elapsed = " << stomp_watch.ElapsedTime()/n_test_points << "s\n";
}

int main(int argc, char **argv) {
  void TreeMapBasicTests();
  void TreeMapPairTests();
  void TreeMapAreaTests();
  void TreeMapRegionTests();
  void TreeMapFieldPairTests();
  void TreeMapNeighborTests();

  std::string usage = "Usage: ";
  usage += argv[0];
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Check that the Stomp::TreeMap class is able to add points and
  // automatically generate sub-pixels.
  if (FLAGS_all_tree_map_tests || FLAGS_tree_map_basic_tests)
    TreeMapBasicTests();

  // Check the Stomp::TreeMap pair-finding routines.
  if (FLAGS_all_tree_map_tests || FLAGS_tree_map_pair_tests)
    TreeMapPairTests();

  // Check that the Stomp::TreeMap class area calculations improve as more
  // points are added to the map.
  if (FLAGS_all_tree_map_tests || FLAGS_tree_map_area_tests)
    TreeMapAreaTests();

  // Check that the Stomp::TreeMap class is able to regionate properly.
  if (FLAGS_all_tree_map_tests || FLAGS_tree_map_region_tests)
    TreeMapRegionTests();

  // Checking pair finding routines with Field values.
  if (FLAGS_all_tree_map_tests || FLAGS_tree_map_field_pair_tests)
    TreeMapFieldPairTests();

  // Checking nearest neighbor routines.
  if (FLAGS_all_tree_map_tests || FLAGS_tree_map_neighbor_tests)
    TreeMapNeighborTests();

  return 0;
}
