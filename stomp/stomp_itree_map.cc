// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the IndexedTreeMap class.  The core work of pair finding
// and K nearest neighbor searches is done in the IndexedTreePixel class.
// However, due to the pixelization scheme used in STOMP, there is a maximum
// pixel size that does not span the entire sphere.  Hence, a vector of
// IndexedTreePixels is necessary to describe an arbitrary collection of points.
// IndexedTreeMap manages that vector of TreePixels, adding them as necessary
// based on the input points.

#include "stomp_core.h"
#include "stomp_itree_map.h"
#include "stomp_map.h"
#include "stomp_angular_bin.h"
#include "stomp_util.h"

namespace Stomp {

IndexedTreeMap::IndexedTreeMap(uint32_t input_resolution,
			       uint16_t maximum_points) {
  resolution_ = input_resolution;
  maximum_points_ = maximum_points;
  point_count_ = 0;
  modified_ = false;
  area_ = 0.0;
  ClearRegions();
}

IndexedTreeMap::IndexedTreeMap(const std::string& input_file,
			       uint32_t input_resolution,
			       uint16_t maximum_points,
			       AngularCoordinate::Sphere sphere,
			       bool verbose, uint8_t theta_column,
			       uint8_t phi_column, int8_t index_column) {
  resolution_ = input_resolution;
  maximum_points_ = maximum_points;
  point_count_ = 0;
  modified_ = false;
  area_ = 0.0;
  ClearRegions();

  Read(input_file, sphere, verbose, theta_column, phi_column, index_column);
}

IndexedTreeMap::~IndexedTreeMap() {
  Clear();
  resolution_ = 0;
  maximum_points_ = 0;
  point_count_ = 0;
  modified_ = false;
  area_ = 0.0;
  ClearRegions();
}

void IndexedTreeMap::FindPairs(AngularCoordinate& ang, AngularBin& theta,
			       IAngularVector& i_angVec) {
  if (!i_angVec.empty()) i_angVec.clear();

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  center_pix.BoundingRadius(ang, theta.ThetaMax(), pix);

  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    ITreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end()) {
      IAngularVector i_ang;
      tree_map_[pix_iter->Pixnum()]->FindPairs(ang, theta, i_ang);
      for (IAngularIterator ang_iter=i_ang.begin();
	   ang_iter!=i_ang.end();++ang_iter) {
	i_angVec.push_back(*ang_iter);
      }
    }
  }
}

void IndexedTreeMap::FindPairs(AngularCoordinate& ang, AngularBin& theta,
			       IndexVector& pair_indices) {
  if (!pair_indices.empty()) pair_indices.clear();

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  center_pix.BoundingRadius(ang, theta.ThetaMax(), pix);

  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    ITreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end()) {
      IndexVector indices;
      tree_map_[pix_iter->Pixnum()]->FindPairs(ang, theta, indices);
      for (IndexIterator idx_iter=indices.begin();
	   idx_iter!=indices.end();++idx_iter) {
	pair_indices.push_back(*idx_iter);
      }
    }
  }
}

void IndexedTreeMap::FindPairs(AngularCoordinate& ang,
			       double theta_min, double theta_max,
			       IAngularVector& i_angVec) {
  AngularBin theta(theta_min, theta_max);
  FindPairs(ang, theta, i_angVec);
}

void IndexedTreeMap::FindPairs(AngularCoordinate& ang,
			       double theta_min, double theta_max,
			       IndexVector& pair_indices) {
  AngularBin theta(theta_min, theta_max);
  FindPairs(ang, theta, pair_indices);
}

void IndexedTreeMap::FindPairs(AngularCoordinate& ang, double theta_max,
			       IAngularVector& i_angVec) {
  AngularBin theta(0.0, theta_max);
  FindPairs(ang, theta, i_angVec);
}

void IndexedTreeMap::FindPairs(AngularCoordinate& ang, double theta_max,
			       IndexVector& pair_indices) {
  AngularBin theta(0.0, theta_max);
  FindPairs(ang, theta, pair_indices);
}

uint16_t IndexedTreeMap::FindKNearestNeighbors(AngularCoordinate& ang,
					       uint8_t n_neighbors,
					       IAngularVector& neighbor_ang) {
  IndexedTreeNeighbor neighbors(ang, n_neighbors);

  _NeighborRecursion(ang, neighbors);

  neighbors.NearestNeighbors(neighbor_ang, false);

  return neighbors.NodesVisited();
}

uint16_t IndexedTreeMap::FindNearestNeighbor(
  AngularCoordinate& ang, IndexedAngularCoordinate& neighbor_ang) {

  IAngularVector angVec;

  uint16_t nodes_visited = FindKNearestNeighbors(ang, 1, angVec);

  neighbor_ang = angVec[0];

  return nodes_visited;
}

double IndexedTreeMap::KNearestNeighborDistance(AngularCoordinate& ang,
						uint8_t n_neighbors,
						uint16_t& nodes_visited) {

  IndexedTreeNeighbor neighbors(ang, n_neighbors);

  _NeighborRecursion(ang, neighbors);

  nodes_visited = neighbors.NodesVisited();

  return neighbors.MaxAngularDistance();
}

double IndexedTreeMap::NearestNeighborDistance(AngularCoordinate& ang,
					       uint16_t& nodes_visited) {
  return KNearestNeighborDistance(ang, 1, nodes_visited);
}

bool IndexedTreeMap::ClosestMatch(AngularCoordinate& ang,
				  double max_distance,
				  IndexedAngularCoordinate& match_ang) {
  IndexedTreeNeighbor neighbors(ang, 1, max_distance);

  _MatchRecursion(ang, neighbors);

  bool found_match = false;
  if (neighbors.Neighbors() == neighbors.MaxNeighbors() &&
      neighbors.MaxAngularDistance() < max_distance) {
    found_match = true;

    IAngularVector neighbor_ang;
    neighbors.NearestNeighbors(neighbor_ang, false);
    match_ang = neighbor_ang[0];
  }

  return found_match;
}

void IndexedTreeMap::_NeighborRecursion(AngularCoordinate& ang,
					IndexedTreeNeighbor& neighbors) {

  // First we need to find out if the input point is within our map area.
  Pixel center_pix(ang, resolution_);
  ITreeDictIterator iter = tree_map_.find(center_pix.Pixnum());

  // If a node containing this point exists, then start finding neighbors there.
  if (iter != tree_map_.end())
    tree_map_[center_pix.Pixnum()]->_NeighborRecursion(ang, neighbors);

  // That should give us back a TreeNeighbor object that contains a workable
  // set of neighbors and a search radius for possible matches.  Now we just
  // need to iterate over those nodes that didn't contain the input point
  // to verify that there can't be any points in their sub-nodes which might
  // be closer to the input point.
  //
  // There's also the possibility that the input point is completely outside
  // our tree.  In that case (where the number of neighbors in the
  // TreeNeighbor object is less than the maximum), we want to check
  // all nodes.
  PixelVector pix;
  if (neighbors.Neighbors() == neighbors.MaxNeighbors()) {
    // We've got a starting list of neighbors, so we only have to look at
    // nodes within our current range.
    center_pix.BoundingRadius(ang, neighbors.MaxAngularDistance(), pix);
  } else {
    // The point is outside of the map area, so we have to check all of the
    // nodes.
    Coverage(pix, resolution_);
  }

  // Now we construct a priority queue so that we're search the nodes closest
  // to the input point first.
  IPixelQueue pix_queue;
  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    ITreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end() && !pix_iter->Contains(ang)) {
      double min_edge_distance, max_edge_distance;
      tree_map_[pix_iter->Pixnum()]->EdgeDistances(ang, min_edge_distance,
						   max_edge_distance);
      DistanceIPixelPair dist_pair(min_edge_distance,
				   tree_map_[pix_iter->Pixnum()]);
      pix_queue.push(dist_pair);
    }
  }

  // And iterate over that queue to check for neighbors.
  while (!pix_queue.empty()) {
    double pix_distance = pix_queue.top().first;
    IndexedTreePixel* pix_iter = pix_queue.top().second;
    if (pix_distance < neighbors.MaxDistance()) {
      pix_iter->_NeighborRecursion(ang, neighbors);
    }
    pix_queue.pop();
  }
}

void IndexedTreeMap::_MatchRecursion(AngularCoordinate& ang,
				     IndexedTreeNeighbor& neighbors) {

  // First we need to find out if the input point is within our map area.
  Pixel center_pix(ang, resolution_);
  ITreeDictIterator iter = tree_map_.find(center_pix.Pixnum());

  // If a node containing this point exists, then start finding matches there.
  if (iter != tree_map_.end())
    tree_map_[center_pix.Pixnum()]->_NeighborRecursion(ang, neighbors);

  // There's also a possibility that the matching point is just on the other
  // side of a pixel boundary.  To see if that's possible, check the edge
  // distance to our reference pixel.
  double min_edge_distance, max_edge_distance;
  center_pix.EdgeDistances(ang, min_edge_distance, max_edge_distance);
  if (min_edge_distance < neighbors.MaxDistance()) {
    // We're near enough to a boundary that we need to check the neighboring
    // pixels.
    PixelVector pix;
    center_pix.BoundingRadius(ang, neighbors.MaxAngularDistance(), pix);

    // Now we construct a priority queue so that we're searching the nodes
    // closest to the input point first.
    IPixelQueue pix_queue;
    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      ITreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end() && !pix_iter->Contains(ang)) {
	tree_map_[pix_iter->Pixnum()]->EdgeDistances(ang, min_edge_distance,
						     max_edge_distance);
	DistanceIPixelPair dist_pair(min_edge_distance,
				     tree_map_[pix_iter->Pixnum()]);
	pix_queue.push(dist_pair);
      }
    }

    // And iterate over that queue to check for neighbors.
    while (!pix_queue.empty()) {
      double pix_distance = pix_queue.top().first;
      IndexedTreePixel* pix_iter = pix_queue.top().second;
      if (pix_distance < neighbors.MaxDistance()) {
	pix_iter->_NeighborRecursion(ang, neighbors);
      }
      pix_queue.pop();
    }
  }
}

bool IndexedTreeMap::AddPoint(IndexedAngularCoordinate* ang) {
  Pixel pix;
  pix.SetResolution(resolution_);
  pix.SetPixnumFromAng(*ang);

  ITreeDictIterator iter = tree_map_.find(pix.Pixnum());
  if (iter == tree_map_.end()) {
    // If we didn't find the pixnum key in the map, then we need to add this
    // pixnum to the map and re-do the search.
    tree_map_.insert(
      std::pair<uint32_t, IndexedTreePixel *>(
	pix.Pixnum(), new IndexedTreePixel(pix.PixelX(), pix.PixelY(),
					   resolution_, maximum_points_)));
    iter = tree_map_.find(pix.Pixnum());
    if (iter == tree_map_.end()) {
      std::cout << "Stomp::IndexedTreeMap::AddPoint - " <<
	"Creating new IndexedTreeMap node failed. Exiting.\n";
      exit(2);
    }
  }
  bool added_point = (*iter).second->AddPoint(ang);
  if (added_point) point_count_++;

  modified_ = true;

  return added_point;
}

bool IndexedTreeMap::AddPoint(IndexedAngularCoordinate& i_ang) {
  IndexedAngularCoordinate* ang_copy =
    new IndexedAngularCoordinate(i_ang.UnitSphereX(), i_ang.UnitSphereY(),
				  i_ang.UnitSphereZ(), i_ang.Index());
  return AddPoint(ang_copy);
}

bool IndexedTreeMap::AddPoint(AngularCoordinate& ang, uint32_t index) {
  IndexedAngularCoordinate* i_ang =
    new IndexedAngularCoordinate(ang.UnitSphereX(), ang.UnitSphereY(),
				 ang.UnitSphereZ(), index);
  return AddPoint(i_ang);
}

bool IndexedTreeMap::Read(const std::string& input_file,
			  AngularCoordinate::Sphere sphere, bool verbose,
			  uint8_t theta_column, uint8_t phi_column,
			  int8_t index_column) {
  bool io_success = false;

  uint32_t n_lines = 0;
  uint32_t check_lines = 128;
  if (theta_column != phi_column) {
    std::ifstream input_file_str(input_file.c_str());

    uint8_t index_idx = static_cast<uint8_t>(index_column);

    if (input_file_str) {
      io_success = true;
      if (verbose) std::cout << "Stomp::IndexedTreeMap::Read - " <<
		     "Reading from " << input_file << "...\n";
      while (!input_file_str.eof()) {
	if (!input_file_str.eof()) {
	  // This should read each line into a buffer, convert that buffer into
	  // a string and then break that string into a vector of strings.  We
	  // should then be able to access theta and phi by converting the
	  // appropriate elements of that vector to doubles.
	  char line_buffer[1000];
	  std::vector<std::string> line_elements;

	  input_file_str.getline(line_buffer, 1000);
	  std::string line_string(line_buffer);
	  Tokenize(line_string, line_elements, " ");
	  n_lines++;
	  if (verbose && n_lines == check_lines) {
	    std::cout << "\tRead " << check_lines << " lines...\n";
	    check_lines *= 2;
	  }

	  if ((line_elements.size() > theta_column) &&
	      (line_elements.size() > phi_column)) {
	    double theta = strtod(line_elements[theta_column].c_str(), NULL);
	    double phi = strtod(line_elements[phi_column].c_str(), NULL);
	    uint32_t index = n_lines - 1;
	    if ((index_column > -1) &&
		(line_elements.size() > index_idx)) {
	      index = static_cast<uint32_t>(
		strtoul(line_elements[index_idx].c_str(), NULL, 10));
	    }

	    IndexedAngularCoordinate* i_ang =
	      new IndexedAngularCoordinate(theta, phi, index, sphere);
	    if (!AddPoint(i_ang)) io_success = false;
	  }
	}
      }
      input_file_str.close();
    } else {
      std::cout << "Stomp::IndexedTreeMap::Read - " << input_file <<
	" does not exist!\n";
    }
  }

  if (verbose && io_success)
    std::cout << "Stomp::IndexedTreeMap::Read - Read " << n_lines-1 <<
      " lines from " << input_file << "; loaded " <<
      NPoints() << " into tree...\n";

  return io_success;
}

void IndexedTreeMap::Coverage(PixelVector& superpix, uint32_t resolution,
			      bool calculate_fraction) {
  if (!superpix.empty()) superpix.clear();

  if (resolution > resolution_) {
    std::cout << "Stomp::IndexedTreeMap::Coverage - " <<
      "WARNING: Requested resolution is higher than " <<
      "the map resolution!\nReseting to map resolution...\n";
    resolution = resolution_;
  }

  // We need to make a vector of pixels that cover the current IndexedTreeMap
  // area.  If the requested resolution is the same as our base node
  // resolution, then this is simple.
  if (resolution_ == resolution) {
    superpix.reserve(tree_map_.size());
    for (ITreeDictIterator iter=tree_map_.begin();
	 iter!=tree_map_.end();++iter) {
      Pixel pix(iter->second->PixelX(), iter->second->PixelY(),
		iter->second->Resolution(), iter->second->Coverage());
      superpix.push_back(pix);
    }
  } else {
    // If that's not the case, then we need to do some work.  First, the case
    // where the requested resolution is coarser than our base nodes.
    if (resolution < resolution_) {
      // We need to find the unique superpixels covered by the map and output
      // those.  We can use a temporary ITreeDict to do this quickly.
      ITreeDict tmp_map;
      for (ITreeDictIterator iter=tree_map_.begin();
	   iter!=tree_map_.end();++iter) {
	IndexedTreePixel* pix =
	  new IndexedTreePixel(iter->second->PixelX(), iter->second->PixelY(),
			       iter->second->Resolution(), 0);
	pix->SetToSuperPix(resolution);
	if (tmp_map.find(pix->Pixnum()) == tmp_map.end()) {
	  tmp_map[pix->Pixnum()] = pix;
	}
      }

      superpix.reserve(tmp_map.size());
      for (ITreeDictIterator iter=tmp_map.begin();
	   iter!=tmp_map.end();++iter) {
	Pixel pix(iter->second->PixelX(), iter->second->PixelY(),
		  iter->second->Resolution(), 1.0);
	if (calculate_fraction) pix.SetWeight(FindUnmaskedFraction(pix));
	superpix.push_back(pix);
	delete iter->second;
      }
    } else {
      // If the requested map is at higher resolution, then we iterate over
      // our map, finding the sub-pixels at the requested resolution.
      for (ITreeDictIterator iter=tree_map_.begin();
	   iter!=tree_map_.end();++iter) {
	PixelVector sub_pix;
	iter->second->SubPix(resolution, sub_pix);

	for (PixelIterator sub_iter=sub_pix.begin();
	     sub_iter!=sub_pix.end();++sub_iter) {
	  // For each of the pixels in the superpixel, we check its status
	  // against the current map.  This is faster than finding the unmasked
	  // fraction directly and immediately tells us which pixels we can
	  // eliminate and which of those we do keep require further
	  // calculations to find the unmasked fraction.
	  int8_t unmasked_status = FindUnmaskedStatus(*sub_iter);
	  if (unmasked_status != 0) {
	    if (calculate_fraction) {
	      sub_iter->SetWeight(FindUnmaskedFraction(*sub_iter));
	    } else {
	      sub_iter->SetWeight(1.0);
	    }
	    superpix.push_back(*sub_iter);
	  }
	}
      }
    }
  }

  // Sort them into the expected order before sending them back since we don't
  // know a priori what the order is coming out of the map object.
  sort(superpix.begin(), superpix.end(), Pixel::SuperPixelBasedOrder);
}

bool IndexedTreeMap::Covering(Map& stomp_map, uint32_t maximum_pixels) {
  if (!stomp_map.Empty()) stomp_map.Clear();

  PixelVector pix;
  Coverage(pix);

  bool met_pixel_requirement;
  if (pix.size() > maximum_pixels) {
    // If the number of requested pixels is smaller than the number of
    // superpixels in the map, then we'd never be able to meet that requirement.
    // In this case, we set the output Map to the coarsest possible case
    // and return false.
    met_pixel_requirement = false;

    stomp_map.Initialize(pix);
  } else {
    // Ok, in this case, we can definitely produce a map that has at most
    // maximum_pixels in it.
    met_pixel_requirement = true;

    // To possibly save ourselves some effort, we start by making the Map
    // equivalent of our base level nodes.
    NodeMap(stomp_map);

    if (maximum_pixels < stomp_map.Size()) {
      // If our Map is still too large, we can use the Map's
      // Covering method to do our work for us.
      Map tmp_map;
      met_pixel_requirement = stomp_map.Covering(tmp_map, maximum_pixels);
      if (tmp_map.Size() < maximum_pixels) {
	stomp_map = tmp_map;
      }
    }
  }

  return met_pixel_requirement;
}

double IndexedTreeMap::FindUnmaskedFraction(Pixel& pix) {
  double unmasked_fraction = 0.0;

  if (pix.Resolution() >= resolution_) {
    // If our input pixel is the size of our base-node or smaller, then we
    // can use each node's Coverage method to find the unmasked fraction,
    // provideded that a matching node can be found.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    ITreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end()) {
      unmasked_fraction = iter->second->Coverage(pix);
    }
  } else {
    // If that's not the case, then we need to find the subpixels of the input
    // pixel that match our base node resolution and iterate over them.
    double pixel_fraction =
      static_cast<double> (pix.Resolution()*pix.Resolution())/
      (resolution_*resolution_);

    PixelVector sub_pix;
    pix.SubPix(resolution_, sub_pix);

    for (PixelIterator sub_iter=sub_pix.begin();
	 sub_iter!=sub_pix.end();++sub_iter) {
      ITreeDictIterator iter = tree_map_.find(sub_iter->Pixnum());
      if (iter != tree_map_.end()) {
	unmasked_fraction += pixel_fraction*iter->second->Coverage(pix);
      }
    }
  }

  return unmasked_fraction;
}

int8_t IndexedTreeMap::FindUnmaskedStatus(Pixel& pix) {
  int8_t unmasked_status = 0;

  // Since we don't have a strong notion of the exact geometry of our map,
  // the return values won't be as exact as they are in the Map or ScalarMap
  // classes.  The important thing, though, is returning a non-zero value if
  // we have reason to believe that there is data in the input pixel.
  if (pix.Resolution() >= resolution_) {
    // If our input pixel is the size of our base-node or smaller, then we
    // just need to find the containing node.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    ITreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end()) {
      unmasked_status = 1;
    }
  } else {
    // If that's not the case, then we need to find the subpixels of the input
    // pixel that match our base node resolution and iterate over them.
    PixelVector sub_pix;
    pix.SubPix(resolution_, sub_pix);
    PixelIterator sub_iter=sub_pix.begin();

    while (sub_iter!=sub_pix.end() && unmasked_status == 0) {
      ITreeDictIterator iter = tree_map_.find(sub_iter->Pixnum());
      if (iter != tree_map_.end()) unmasked_status = -1;

      ++sub_iter;
    }
  }

  return unmasked_status;
}

void IndexedTreeMap::NodeMap(Map& stomp_map) {
  if (!stomp_map.Empty()) stomp_map.Clear();

  PixelVector pix;
  pix.reserve(tree_map_.size());
  for (ITreeDictIterator iter=tree_map_.begin();iter!=tree_map_.end();++iter) {
    Pixel tmp_pix(iter->second->PixelX(), iter->second->PixelY(),
		   iter->second->Resolution(), 1.0);
    pix.push_back(tmp_pix);
  }

  stomp_map.Initialize(pix);
}

uint32_t IndexedTreeMap::Resolution() {
  return resolution_;
}

uint16_t IndexedTreeMap::PixelCapacity() {
  return maximum_points_;
}

void IndexedTreeMap::SetResolution(uint32_t resolution) {
  Clear();
  resolution_ = resolution;
}

void IndexedTreeMap::SetPixelCapacity(int pixel_capacity) {
  Clear();
  maximum_points_ = pixel_capacity;
}

uint32_t IndexedTreeMap::NPoints(uint32_t k) {
  return (k == MaxPixnum ? point_count_ :
	  (tree_map_.find(k) != tree_map_.end() ?
	   tree_map_[k]->NPoints() : 0));
}

uint32_t IndexedTreeMap::NPoints(Pixel& pix) {
  uint32_t total_points = 0;

  // First we check to see if the input pixel is larger than our base-nodes.
  if (pix.Resolution() < resolution_) {
    // If this is the case, then we need to iterate over all of the sub-pixels
    // at the base-node resolution to see if there is any overlap and return
    // the sum.  We can just take the total since we know that the sub-pixel
    // covers the entire node.
    PixelVector pixVec;
    pix.SubPix(resolution_, pixVec);
    for (PixelIterator pix_iter=pixVec.begin();
	 pix_iter!=pixVec.end();++pix_iter) {
      ITreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end())
	total_points += tree_map_[pix_iter->Pixnum()]->NPoints();
    }
  } else {
    // If the input pixel is the same size as our nodes or smaller, then we
    // look for the appropriate node and return the results from that node.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    ITreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end())
      total_points = tree_map_[tmp_pix.Pixnum()]->NPoints(pix);
  }

  return total_points;
}

void IndexedTreeMap::Points(IAngularVector& i_ang) {
  if (!i_ang.empty()) i_ang.clear();

  // Fairly simple, just iterate over all of the base-nodes and return the
  // aggregate vector.  We can speed things up a bit since we know the total
  // number of points in the map.
  i_ang.reserve(point_count_);

  for (ITreeDictIterator iter=tree_map_.begin();
       iter!=tree_map_.end();++iter) {
    IAngularVector tmp_ang;
    iter->second->Points(tmp_ang);
    for (IAngularIterator ang_iter=tmp_ang.begin();
	 ang_iter!=tmp_ang.end();++ang_iter) i_ang.push_back(*ang_iter);
  }
}

void IndexedTreeMap::Points(IAngularVector& i_ang, Pixel& pix) {
  if (!i_ang.empty()) i_ang.clear();

  // First we check to see if the input pixel is larger than our base-nodes.
  if (pix.Resolution() < resolution_) {
    // If this is the case, then we need to iterate over all of the sub-pixels
    // at the base-node resolution to see if there is any overlap and return
    // the aggregate.  We can just take all of the points in each node since
    // we know that the sub-pixel covers the entire node.
    PixelVector pixVec;
    pix.SubPix(resolution_, pixVec);
    for (PixelIterator pix_iter=pixVec.begin();
	 pix_iter!=pixVec.end();++pix_iter) {
      ITreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end()) {
	IAngularVector tmp_ang;

	iter->second->Points(tmp_ang);
	for (IAngularIterator ang_iter=tmp_ang.begin();
	     ang_iter!=tmp_ang.end();++ang_iter) i_ang.push_back(*ang_iter);
      }
    }
  } else {
    // If the input pixel is the same size as our nodes or smaller, then we
    // look for the appropriate node and return the results from that node.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    ITreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end()) {
      IAngularVector tmp_ang;

      iter->second->Points(tmp_ang, pix);
      for (IAngularIterator ang_iter=tmp_ang.begin();
	   ang_iter!=tmp_ang.end();++ang_iter) i_ang.push_back(*ang_iter);
    }
  }
}

void IndexedTreeMap::Indices(Pixel& pix, IndexVector& indices) {
  if (!indices.empty()) indices.clear();

  // First we check to see if the input pixel is larger than our base-nodes.
  if (pix.Resolution() < resolution_) {
    // If this is the case, then we need to iterate over all of the sub-pixels
    // at the base-node resolution to see if there is any overlap and return
    // the sum.  We can just take the total since we know that the sub-pixel
    // covers the entire node.
    PixelVector pixVec;
    pix.SubPix(resolution_, pixVec);
    for (PixelIterator pix_iter=pixVec.begin();
	 pix_iter!=pixVec.end();++pix_iter) {
      ITreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end()) {
	IndexVector tmp_indices;
	tree_map_[pix_iter->Pixnum()]->Indices(pix, tmp_indices);
	for (IndexIterator idx_iter=tmp_indices.begin();
	     idx_iter!=tmp_indices.end();++idx_iter)
	  indices.push_back(*idx_iter);
      }
    }
  } else {
    // If the input pixel is the same size as our nodes or smaller, then we
    // look for the appropriate node and return the results for that node.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    ITreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end())
      tree_map_[tmp_pix.Pixnum()]->Indices(pix, indices);
  }
}

uint16_t IndexedTreeMap::BaseNodes() {
  return tree_map_.size();
}

uint16_t IndexedTreeMap::Nodes() {
  uint16_t total_nodes = 0;
  for (ITreeDictIterator iter=tree_map_.begin();
       iter!=tree_map_.end();++iter) total_nodes += iter->second->Nodes();

  return total_nodes;
}

uint32_t IndexedTreeMap::Size() {
  return point_count_;
}

double IndexedTreeMap::Area() {
  if (modified_) CalculateArea();

  return area_;
}

void IndexedTreeMap::CalculateArea() {
  area_ = 0.0;
  for (ITreeDictIterator iter=tree_map_.begin();
       iter!=tree_map_.end();++iter) {
    area_ += iter->second->Coverage()*iter->second->Area();
  }

  modified_ = false;
}

uint32_t IndexedTreeMap::MinResolution() {
  return resolution_;
}

uint32_t IndexedTreeMap::MaxResolution() {
  return resolution_;
}

uint8_t IndexedTreeMap::MinLevel() {
  return Pixel::ResolutionToLevel(resolution_);
}

uint8_t IndexedTreeMap::MaxLevel() {
  return Pixel::ResolutionToLevel(resolution_);
}

bool IndexedTreeMap::Empty() {
  return (tree_map_.empty() ? true : false);
}

void IndexedTreeMap::Clear() {
  if (!tree_map_.empty()) {
    for (ITreeDictIterator iter=tree_map_.begin();
	 iter!=tree_map_.end();++iter) {
      iter->second->Clear();
      delete iter->second;
    }
    tree_map_.clear();
    area_ = 0.0;
    point_count_ = 0;
    modified_ = false;
  }
  ClearRegions();
}

} // end namespace Stomp
