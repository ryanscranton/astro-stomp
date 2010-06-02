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
#include "stomp_itree_map.h"

void IndexedTreeMapBasicTests() {
  std::cout << "\n";
  std::cout << "**********************************\n";
  std::cout << "*** IndexedTreeMap Basic Tests ***\n";
  std::cout << "**********************************\n";
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a non-standard resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint32_t resolution = 32;
  Stomp::IndexedTreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::IndexedTreeMap at " << resolution <<
    " resolution...\n";

  if (tree_map.AddPoint(ang, 0)) {
    std::cout << "\tSuccessfully added a random position to the map.\n";
  } else {
    std::cout << "\tFailed to add a random position to the map.\n";
    exit(1);
  }

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::IndexedTreeMap.  We choose a large
  // enough radius to make sure that we're adding points that are well outside
  // of pixel that contained the initial point to make sure that we're properly
  // adding new pixels as we ingest points.
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
  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    added_point = tree_map.AddPoint(*iter, idx);
    if (!added_point)
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
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

  // Now we verify that the superpixels for our IndexedTreeMap match those
  // for the Map that we built it from.
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

void IndexedTreeMapPairTests() {
  std::cout << "\n";
  std::cout << "*********************************\n";
  std::cout << "*** IndexedTreeMap Pair Tests ***\n";
  std::cout << "*********************************\n";
  uint16_t n_points_per_node = 50;
  uint32_t resolution = 4;
  Stomp::IndexedTreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::IndexedTreeMap at " << resolution <<
    " resolution with " << n_points_per_node << " points per node...\n";
  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::IndexedTreeMap.  We choose a large
  // enough radius to make sure that we're adding points that are well outside
  // of pixel that contained the initial point to make sure that we're properly
  // adding new pixels as we ingest points.
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

  bool added_point = false;
  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    added_point = tree_map.AddPoint(*iter, idx);
    if (!added_point)
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
  }
  std::cout << "\t" << tree_map.BaseNodes() << " base nodes at " <<
    tree_map.Resolution() << " resolution; " << tree_map.Nodes() <<
    " total nodes.\n";

  // We start by choosing a radius much larger than our area, which should mean
  // that we find all points in the map as pairs.
  std::cout << "\nEncompassing radius test:\n";
  Stomp::IndexVector pair_indices;
  tree_map.FindPairs(ang, 4.0*theta_radius, pair_indices);
  std::cout << "\tFound " << pair_indices.size() <<
    " pairs (expect " << n_points << ")\n";

  // Next, we choose a smaller radius where we should begin to probe some of
  // the tree structure in the pixel.
  double theta = 0.1;
  uint32_t direct_pairs = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    if (ang.AngularDistance(*iter) <= theta) direct_pairs++;
  }
  std::cout << "\nSmall circle test:\n";
  tree_map.FindPairs(ang, theta, pair_indices); 
  std::cout << "\tFound " << pair_indices.size() << " pairs; " <<
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
  tree_map.FindPairs(ang, theta_min, theta_max, pair_indices);
  std::cout << "\tFound " << pair_indices.size() <<
    " pairs; " << direct_pairs << " pairs from brute force.\n";

  // Now we test to see if the pair finding works at the edges of our map.
  // First, the top of the map, where the angular distortion should be
  // maximum.
  double scale_factor = 1.0;
  ang.SetSurveyCoordinates(lambda+scale_factor*theta_radius,eta);
  while (stomp_map->Contains(ang)) {
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
  tree_map.FindPairs(ang, theta_min, theta_max, pair_indices);
  std::cout << "\tFound " << pair_indices.size() << " pairs; " <<
    direct_pairs << " pairs from brute force.\n";

  // Now, a location close to the right edge of the circle.
  scale_factor = 1.0;
  ang.SetSurveyCoordinates(lambda,eta+scale_factor*theta_radius);
  while (stomp_map->Contains(ang)) {
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
  tree_map.FindPairs(ang, theta_min, theta_max, pair_indices);
  std::cout << "\tFound " << pair_indices.size() << " pairs; " <<
    direct_pairs << " pairs from brute force.\n";
}

void IndexedTreeMapAreaTests() {
  std::cout << "\n";
  std::cout << "*********************************\n";
  std::cout << "*** IndexedTreeMap Area Tests ***\n";
  std::cout << "*********************************\n";
  // The goal here is to check that the estimate of the area associated with
  // the IndexedTreeMap improves as we add points.

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a coarse resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint32_t resolution = 4;
  Stomp::IndexedTreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::IndexedTreeMap at " << resolution <<
    " resolution...\n";

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::IndexedTreeMap.  We choose a large
  // enough radius to make sure that we're adding points that are well outside
  // of pixel that contained the initial point to make sure that we're properly
  // adding new pixels as we ingest points.
  double theta = 5.0;
  uint32_t annulus_resolution = 32;
  Stomp::Pixel tmp_pix(ang, annulus_resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Now we add an increasing number of points and track how the area of our
  // IndexedTreeMap compares with the area of the Map used to generate the
  // random points.  Since we used a finer resolution to make the Map than
  // for our IndexedTreeMap, we should start with a larger IndexedTreeMap area
  // and then eventually get closer to the Map area as we add points.
  uint32_t idx = 0;
  for (uint32_t k=1;k<200;k*=2) {
    uint32_t n_points = 5000*k;

    Stomp::AngularVector angVec;
    stomp_map->GenerateRandomPoints(angVec, n_points);

    bool added_point = false;
    for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
      added_point = tree_map.AddPoint(*iter, idx);
      if (!added_point)
	std::cout << "\t\tFailed to add point: " <<
	  iter->RA() << ", " << iter->DEC() << "\n";
      idx++;
    }

    std::cout << "\t" << tree_map.NPoints() << " points: " <<
      tree_map.BaseNodes() << " nodes at " << tree_map.Resolution() <<
      " resolution; " << tree_map.Nodes() << " total nodes.\n";
    std::cout << "\t\tInput map area: " << stomp_map->Area() <<
      "; IndexedTreeMap area: " << tree_map.Area() << "\n";
  }
}

void IndexedTreeMapRegionTests() {
  std::cout << "\n";
  std::cout << "***********************************\n";
  std::cout << "*** IndexedTreeMap Region Tests ***\n";
  std::cout << "***********************************\n";
  // The goal here is to check that the regionation code works, provided that
  // we have a dense sampling of our data area.

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a coarse resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint32_t resolution = 4;
  Stomp::IndexedTreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::IndexedTreeMap at " << resolution <<
    " resolution...\n";

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::IndexedTreeMap.  We choose a large
  // enough radius to make sure that we're adding points that are well outside
  // of pixel that contained the initial point to make sure that we're properly
  // adding new pixels as we ingest points.
  double theta = 5.0;
  uint32_t annulus_resolution = 32;
  Stomp::Pixel tmp_pix(ang, annulus_resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Now we add a high number of points to the IndexedTreeMap to flesh out
  // its bounds.
  uint32_t n_points = 1000000;
  Stomp::AngularVector angVec;
  stomp_map->GenerateRandomPoints(angVec, n_points);

  bool added_point = false;
  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    added_point = tree_map.AddPoint(*iter, idx);
    if (!added_point)
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
  }
  angVec.clear();

  std::cout << "\t" << tree_map.NPoints() << " points: " <<
    tree_map.Nodes() << " nodes at " << tree_map.Resolution() <<
    " resolution used for this map\n";
  std::cout << "\t\tInput map area: " << stomp_map->Area() <<
    "; IndexedTreeMap area: " << tree_map.Area() << "\n";

  // Before regionating our IndexedTreeMap, let's verify that the unmasked
  // fraction calculations that go into the Coverage maps are calculated
  // correctly.
  Stomp::PixelVector map_coverage;
  stomp_map->Coverage(map_coverage);

  Stomp::PixelVector tree_coverage;
  tree_map.Coverage(tree_coverage);

  if (map_coverage.size() != tree_coverage.size()) {
    std::cout << "Disagreement between Map Coverage (" <<
      map_coverage.size() << " pixels) and IndexedTreeMap Coverage (" <<
      tree_coverage.size() << " pixels).\n";
  } else {
    std::cout << "Unmasked fractions:\n";
    for (uint32_t i=0;i<map_coverage.size();i++)
      std::cout << "\t" << map_coverage[i].Pixnum() << ": Map = " <<
	map_coverage[i].Weight() << ", IndexedTreeMap = " <<
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
      n_tree_regions << " for the IndexedTreeMap.\n";
  } else {
    std::cout << "\tGot " << n_map_regions << " regions for both maps:\n";
    for (uint16_t i=0;i<n_map_regions;i++) {
      std::cout << "\t\tRegion " << i << ": Map area = " <<
	stomp_map->RegionArea(i) << ", IndexedTreeMap area = " <<
	tree_map.RegionArea(i) << "\n";
    }
  }
}

void IndexedTreeMapNeighborTests() {
  // Checking nearest neighbor finding routines.
  std::cout << "\n";
  std::cout << "**************************************\n";
  std::cout << "*** IndexedTreeMap Nearest Neighbor Tests ***\n";
  std::cout << "**************************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a coarse resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint32_t resolution = 8;
  Stomp::IndexedTreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::IndexedTreeMap at " << resolution <<
    " resolution...\n";

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::IndexedTreeMap.  We choose a large
  // enough radius to make sure that we're adding points that are well outside
  // of pixel that contained the initial point to make sure that we're properly
  // adding new pixels as we ingest points.
  double theta_bound = 5.0;
  uint32_t annulus_resolution = 32;
  Stomp::Pixel tmp_pix(ang, annulus_resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta_bound, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Now we add a high number of points to the IndexedTreeMap to flesh out
  // its bounds.
  uint32_t n_points = 100000;
  Stomp::AngularVector angVec;
  Stomp::IAngularVector i_angVec;
  std::cout << "Adding " << n_points << " points\n";
  stomp_map->GenerateRandomPoints(angVec, n_points);

  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    Stomp::IndexedAngularCoordinate tmp_ang(iter->UnitSphereX(),
					    iter->UnitSphereY(),
					    iter->UnitSphereZ(), idx);
    if (!tree_map.AddPoint(tmp_ang))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
  }

  // Quick global accounting check.
  std::cout << "\t" << tree_map.NPoints() << " points added; " <<
    1.0*tree_map.NPoints()/tree_map.Area() << " points/sq. degree.\n";
  std::cout << "\t" << tree_map.BaseNodes() << " base nodes at " <<
    tree_map.Resolution() << " resolution; " << tree_map.Nodes() <<
    " total nodes.\n";

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
  uint16_t failed_matches = 0;
  double mean_neighbor_distance = 0.0;
  double mean_nodes_visited = 0.0;
  Stomp::IndexedAngularCoordinate nearest_neighbor;
  stomp_watch.StartTimer();
  for (uint32_t i=0;i<n_test_points;i++) {
    nodes_visited = tree_map.FindNearestNeighbor(angVec[i], nearest_neighbor);
    mean_neighbor_distance += angVec[i].AngularDistance(nearest_neighbor);
    mean_nodes_visited += static_cast<double>(nodes_visited);
    if (i != nearest_neighbor.Index()) failed_matches++;
  }
  stomp_watch.StopTimer();

  std::cout << "\tMean distance to nearest neighbor = " <<
    mean_neighbor_distance/n_test_points << "\n" <<
    "\t\tFailed matches: " << failed_matches << "\n" <<
    "\t\tMean nodes visited = " << mean_nodes_visited/n_test_points <<
    "/" << total_nodes << "\n\t\tTime elapsed = " <<
    stomp_watch.ElapsedTime()/n_test_points << "s\n";

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

void IndexedTreeMapMatchTests() {
  // Checking closest match finding routines.
  std::cout << "\n";
  std::cout << "***********************************\n";
  std::cout << "*** IndexedTreeMap Closest Match Tests ***\n";
  std::cout << "***********************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  uint16_t n_points_per_node = 200;
  // Make the map at a coarse resolution so that we can test the
  // Coverage and NodeMap methods later on.
  uint32_t resolution = 8;
  Stomp::IndexedTreeMap tree_map(resolution, n_points_per_node);
  std::cout << "Building Stomp::IndexedTreeMap at " << resolution <<
    " resolution...\n";

  // Now we generate a bunch of random points around our original point and
  // attempt to add them to the Stomp::IndexedTreeMap.  We choose a large
  // enough radius to make sure that we're adding points that are well outside
  // of pixel that contained the initial point to make sure that we're properly
  // adding new pixels as we ingest points.
  double theta_bound = 5.0;
  uint32_t annulus_resolution = 32;
  Stomp::Pixel tmp_pix(ang, annulus_resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta_bound, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Now we add a high number of points to the IndexedTreeMap to flesh out
  // its bounds.
  uint32_t n_points = 100000;
  Stomp::AngularVector angVec;
  Stomp::WAngularVector w_angVec;
  std::cout << "Adding " << n_points << " points\n";
  stomp_map->GenerateRandomPoints(angVec, n_points);

  uint32_t idx = 0;
  for (Stomp::AngularIterator iter=angVec.begin();iter!=angVec.end();++iter) {
    Stomp::IndexedAngularCoordinate tmp_ang(iter->UnitSphereX(),
					    iter->UnitSphereY(),
					    iter->UnitSphereZ(), idx);
    if (!tree_map.AddPoint(tmp_ang))
      std::cout << "\t\tFailed to add point: " <<
	iter->RA() << ", " << iter->DEC() << "\n";
    idx++;
  }

  // Quick global accounting check.
  std::cout << "\t" << tree_map.NPoints() << " points added; " <<
    1.0*tree_map.NPoints()/tree_map.Area() << " points/sq. degree.\n";
  std::cout << "\t" << tree_map.BaseNodes() << " base nodes at " <<
    tree_map.Resolution() << " resolution; " << tree_map.Nodes() <<
    " total nodes.\n";

  // Generate a smaller set of test points.
  uint32_t n_test_points = 1000;
  Stomp::AngularVector test_angVec;
  stomp_map->GenerateRandomPoints(test_angVec, n_test_points);

  // Now, we check the annulus finding.
  Stomp::StompWatch stomp_watch;
  Stomp::IndexedAngularCoordinate test_point;

  double match_radius = 3.0;
  std::cout <<
    "\nFinding closest matches using " << n_test_points <<
    " points in the tree and " << match_radius << " arcsecond radius...\n";
  stomp_watch.StartTimer();
  uint16_t n_match = 0;
  double mean_distance = 0.0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_map.ClosestMatch(angVec[i], match_radius/3600.0, test_point) &&
	i == test_point.Index()) {
      n_match++;
      mean_distance += angVec[i].AngularDistance(test_point);
    }
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points <<
    " points; Average distance = " << 3600.0*mean_distance/n_match <<
    " arcseconds.\n" << "\tTime elapsed = " <<
    stomp_watch.ElapsedTime()/n_test_points << "s\n";

  match_radius = 0.5;
  std::cout <<
    "\nFinding closest matches using " << n_test_points <<
    " points in the tree and " << match_radius << " arcsecond radius..\n";
  stomp_watch.StartTimer();
  n_match = 0;
  mean_distance = 0.0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_map.ClosestMatch(angVec[i], match_radius/3600.0, test_point)) {
      n_match++;
      mean_distance += angVec[i].AngularDistance(test_point);
    }
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points <<
    " points; Average distance = " << 3600.0*mean_distance/n_match <<
    " arcseconds.\n" << "\tTime elapsed = " <<
    stomp_watch.ElapsedTime()/n_test_points << "s\n";

  match_radius = 3.0;
  std::cout << "\nFinding matches using " << n_test_points <<
    " points in the pixel with " << match_radius << " arcsecond radius...\n";
  stomp_watch.StartTimer();
  n_match = 0;
  mean_distance = 0.0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_map.ClosestMatch(test_angVec[i],
			      match_radius/3600.0, test_point)) {
      n_match++;
      mean_distance += test_angVec[i].AngularDistance(test_point);
      if (i < 5) {
	std::cout << "(" << test_angVec[i].Lambda() <<
	  "," << test_angVec[i].Eta() <<
	  ") -> (" << test_point.Lambda() << "," << test_point.Eta() <<
	  "): " << test_angVec[i].AngularDistance(test_point) << "\n";
      }
    }
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points <<
    " points; Average distance = " << 3600.0*mean_distance/n_match <<
    " arcseconds.\n" << "\tTime elapsed = " <<
    stomp_watch.ElapsedTime()/n_test_points << "s\n";

  match_radius = 0.5;
  std::cout << "\nFinding matches using " << n_test_points <<
    " points in the pixel with " << match_radius << " arcsecond radius...\n";
  stomp_watch.StartTimer();
  n_match = 0;
  mean_distance = 0.0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_map.ClosestMatch(test_angVec[i],
			      match_radius/3600.0, test_point)) {
      n_match++;
      mean_distance += test_angVec[i].AngularDistance(test_point);
    }
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points <<
    " points; Average distance = " << 3600.0*mean_distance/n_match <<
    " arcseconds.\n" << "\tTime elapsed = " <<
    stomp_watch.ElapsedTime()/n_test_points << "s\n";

  match_radius = 3.0;
  std::cout << "\nFinding matches using " << n_test_points <<
    " points outside the pixel and " << match_radius <<
    " arcsecond radius...\n";
  for (uint16_t i=0;i<n_test_points;i++)
    test_angVec[i].SetSurveyCoordinates(test_angVec[i].Lambda()-40.0,
					test_angVec[i].Eta());

  stomp_watch.StartTimer();
  n_match = 0;
  mean_distance = 0.0;
  for (uint16_t i=0;i<n_test_points;i++) {
    if (tree_map.ClosestMatch(test_angVec[i],
			      match_radius/3600.0, test_point)) {
      n_match++;
      mean_distance += test_angVec[i].AngularDistance(test_point);
    }
  }
  stomp_watch.StopTimer();

  std::cout << "\tFound " << n_match << "/" << n_test_points <<
    " points; Average distance = " << 3600.0*mean_distance/n_match <<
    " arcseconds.\n" << "\tTime elapsed = " <<
    stomp_watch.ElapsedTime()/n_test_points << "s\n";
}

// Define our command line flags here so we can use these flags later.
DEFINE_bool(all_itree_map_tests, false, "Run all class unit tests.");
DEFINE_bool(itree_map_basic_tests, false, "Run IndexedTreeMap basic tests");
DEFINE_bool(itree_map_pair_tests, false, "Run IndexedTreeMap pair tests");
DEFINE_bool(itree_map_area_tests, false, "Run IndexedTreeMap area tests");
DEFINE_bool(itree_map_region_tests, false, "Run IndexedTreeMap region tests");
DEFINE_bool(itree_map_neighbor_tests, false,
            "Run IndexedTreeMap nearest neighbor tests");
DEFINE_bool(itree_map_match_tests, false,
            "Run IndexedTreeMap closest match tests");

void IndexedTreeMapUnitTests(bool run_all_tests) {
  void IndexedTreeMapBasicTests();
  void IndexedTreeMapPairTests();
  void IndexedTreeMapAreaTests();
  void IndexedTreeMapRegionTests();
  void IndexedTreeMapNeighborTests();
  void IndexedTreeMapMatchTests();

  if (run_all_tests) FLAGS_all_itree_map_tests = true;

  // Check that the Stomp::IndexedTreeMap class is able to add points and
  // automatically generate sub-pixels.
  if (FLAGS_all_itree_map_tests || FLAGS_itree_map_basic_tests)
    IndexedTreeMapBasicTests();

  // Check the Stomp::IndexedTreeMap pair-finding routines.
  if (FLAGS_all_itree_map_tests || FLAGS_itree_map_pair_tests)
    IndexedTreeMapPairTests();

  // Check that the Stomp::IndexedTreeMap class area calculations improve
  // as more points are added to the map.
  if (FLAGS_all_itree_map_tests || FLAGS_itree_map_area_tests)
    IndexedTreeMapAreaTests();

  // Check that the Stomp::IndexedTreeMap class is able to regionate properly.
  if (FLAGS_all_itree_map_tests || FLAGS_itree_map_region_tests)
    IndexedTreeMapRegionTests();

  // Checking nearest neighbor routines.
  if (FLAGS_all_itree_map_tests || FLAGS_itree_map_neighbor_tests)
    IndexedTreeMapNeighborTests();

  // Checking closest match routines.
  if (FLAGS_all_itree_map_tests || FLAGS_itree_map_match_tests)
    IndexedTreeMapMatchTests();
}
