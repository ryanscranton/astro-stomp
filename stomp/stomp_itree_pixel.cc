// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains a variant on the TreePixel class.  Again, we are
// using the hierarchical pixelization as the scaffolding for a spatial tree
// structure.  However, in this case, the tree is populated by
// IndexedAngularCoordinate objects rather than their weighted cousins.  For
// integer data types, we're not going to be interested in aggregate statistics
// like with the floating point data, but rather returning vectors of objects
// or indices.

#include "stomp_core.h"
#include "stomp_itree_pixel.h"
#include "stomp_angular_bin.h"

namespace Stomp {

IndexedTreePixel::IndexedTreePixel() {
  SetWeight(0.0);
  maximum_points_ = 0;
  point_count_ = 0;
  InitializeCorners();
}

IndexedTreePixel::IndexedTreePixel(const uint32_t input_resolution,
				   const uint32_t input_pixnum,
				   const uint16_t maximum_points) {
  SetResolution(input_resolution);

  uint32_t tmp_y = input_pixnum/(Stomp::Nx0*Resolution());
  uint32_t tmp_x = input_pixnum - Stomp::Nx0*Resolution()*tmp_y;

  SetPixnumFromXY(tmp_x, tmp_y);
  SetWeight(0.0);
  maximum_points_ = maximum_points;
  point_count_ = 0;
  initialized_subpixels_ = false;
  InitializeCorners();
}

IndexedTreePixel:: IndexedTreePixel(const uint32_t input_x,
				    const uint32_t input_y,
				    const uint32_t input_resolution,
				    const uint16_t maximum_points) {
  SetResolution(input_resolution);
  SetPixnumFromXY(input_x, input_y);
  SetWeight(0.0);
  maximum_points_ = maximum_points;
  point_count_ = 0;
  initialized_subpixels_ = false;
  InitializeCorners();
}

IndexedTreePixel::IndexedTreePixel(AngularCoordinate& ang,
				   const uint32_t input_resolution,
				   const uint16_t maximum_points) {
  SetResolution(input_resolution);
  SetPixnumFromAng(ang);
  SetWeight(0.0);
  maximum_points_ = maximum_points;
  point_count_ = 0;
  initialized_subpixels_ = false;
  InitializeCorners();
}

IndexedTreePixel::~IndexedTreePixel() {
  ang_.clear();
  subpix_.clear();
  maximum_points_ = 0;
  point_count_ = 0;
  initialized_subpixels_ = false;
}

bool IndexedTreePixel::_InitializeSubPixels() {
  initialized_subpixels_ = false;
  // If we're already at the maximum resolution, then we shouldn't be trying
  // to generate sub-pixels.
  if (Resolution() < Stomp::MaxPixelResolution) {
    PixelVector tmp_pix;
    SubPix(Resolution()*2, tmp_pix);
    subpix_.reserve(4);

    // Provided we passed that test, we create a vector of sub-pixels.
    for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter) {
      IndexedTreePixel* tree_pix =
	new IndexedTreePixel(iter->PixelX(), iter->PixelY(),
			     iter->Resolution(), maximum_points_);
      subpix_.push_back(tree_pix);
    }
    initialized_subpixels_ = true;

    // Now we iterate over all of the AngularCoordinates in the current pixel.
    // Provided that we find a home for all of them, we return true.  If any
    // of them fail to fit into a sub-pixel, then we return false.
    bool transferred_point_to_subpixels = false;
    for (uint32_t i=0;i<ang_.size();++i) {
      transferred_point_to_subpixels = false;
      for (uint32_t j=0;j<subpix_.size();++j) {
	if (subpix_[j]->AddPoint(ang_[i])) {
	  j = subpix_.size();
	  transferred_point_to_subpixels = true;
	}
      }
      if (!transferred_point_to_subpixels) initialized_subpixels_ = false;
    }
    ang_.clear();
  }

  return initialized_subpixels_;
}

void IndexedTreePixel::FindPairs(AngularCoordinate& ang, AngularBin& theta,
				 IAngularVector& i_angVec) {
  if (!i_angVec.empty()) i_angVec.clear();

  IAngularPtrVector i_ang;

  // If we have AngularCoordinates in this pixel, then this is just a
  // matter of iterating through them and finding how many satisfy the
  // angular bounds.
  if (!ang_.empty()) {
    if (theta.ThetaMax() < 90.0) {
      for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	if (theta.WithinCosBounds((*iter)->DotProduct(ang)))
	  i_ang.push_back(*iter);
      }
    } else {
      for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	if (theta.WithinBounds((*iter)->AngularDistance(ang)))
	  i_ang.push_back(*iter);
      }
    }
  } else {
    // If the current pixel doesn't contain any points, then we need to see
    // if either the current pixel is either fully or partially contained in
    // the annulus.  For the former case, we can just send back the total
    // number of points in this pixel.  In the latter case, we pass things
    // along to the sub-pixels.  If neither of those things are true, then
    // we're done and we send back zero.
    int8_t intersects_annulus = IntersectsAnnulus(ang, theta);

    if (intersects_annulus == 1) {
      // Fully contained in the annulus, so we get a copy of all of the pointers
      // in the sub-nodes.
      for (ITreePtrIterator iter=subpix_.begin();
	   iter!=subpix_.end();++iter) {
	_PointPtrs(i_ang);
      }
    } else {
      if (intersects_annulus == -1) {
	// Partial intersection with the annulus, so to pass the bounds along
	// to the sub-nodes.
	for (ITreePtrIterator iter=subpix_.begin();
	     iter!=subpix_.end();++iter) {
	  (*iter)->_FindPairs(ang, theta, i_ang);
	}
      }
    }
  }

  // At this point, we have a vector of pointers to our IndexedAngularCoordinate
  // objects.  We don't want to pass those outside of the tree and risk someone
  // breaking our tree by deleting those pointers, so we make a vector of
  // copies and pass that back.
  if (!i_ang.empty()) {
    i_angVec.reserve(i_ang.size());

    for (IAngularPtrIterator iter=i_ang.begin();iter!=i_ang.end();++iter) {
      i_angVec.push_back(IndexedAngularCoordinate((*iter)->UnitSphereX(),
						  (*iter)->UnitSphereY(),
						  (*iter)->UnitSphereZ(),
						  (*iter)->Index()));
    }
  }
}

void IndexedTreePixel::FindPairs(AngularCoordinate& ang, AngularBin& theta,
				 IndexVector& pair_indices) {
  if (!pair_indices.empty()) pair_indices.clear();

  IAngularVector i_ang;
  FindPairs(ang, theta, i_ang);

  if (!i_ang.empty()) {
    pair_indices.reserve(i_ang.size());

    for (IAngularIterator iter=i_ang.begin();iter!=i_ang.end();++iter) {
      pair_indices.push_back(iter->Index());
    }
  }
}


void IndexedTreePixel::FindPairs(AngularCoordinate& ang,
				 double theta_min, double theta_max,
				 IAngularVector& i_angVec) {
  AngularBin theta(theta_min, theta_max);
  FindPairs(ang, theta, i_angVec);
}

void IndexedTreePixel::FindPairs(AngularCoordinate& ang,
				 double theta_min, double theta_max,
				 IndexVector& pair_indices) {
  AngularBin theta(theta_min, theta_max);
  FindPairs(ang, theta, pair_indices);
}

void IndexedTreePixel::FindPairs(AngularCoordinate& ang, double theta_max,
				 IAngularVector& i_angVec) {
  AngularBin theta(0.0, theta_max);
  FindPairs(ang, theta, i_angVec);
}

void IndexedTreePixel::FindPairs(AngularCoordinate& ang, double theta_max,
				 IndexVector& pair_indices) {
  AngularBin theta(0.0, theta_max);
  FindPairs(ang, theta, pair_indices);
}

void IndexedTreePixel::_FindPairs(AngularCoordinate& ang, AngularBin& theta,
				  IAngularPtrVector& i_ang) {
  // If we have AngularCoordinates in this pixel, then this is just a
  // matter of iterating through them and finding how many satisfy the
  // angular bounds.
  if (!ang_.empty()) {
    if (theta.ThetaMax() < 90.0) {
      for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	if (theta.WithinCosBounds((*iter)->DotProduct(ang)))
	  i_ang.push_back(*iter);
      }
    } else {
      for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	if (theta.WithinBounds((*iter)->AngularDistance(ang)))
	  i_ang.push_back(*iter);
      }
    }
  } else {
    // If the current pixel doesn't contain any points, then we need to see
    // if either the current pixel is either fully or partially contained in
    // the annulus.  For the former case, we can just send back the total
    // number of points in this pixel.  In the latter case, we pass things
    // along to the sub-pixels.  If neither of those things are true, then
    // we're done and we send back zero.
    int8_t intersects_annulus = IntersectsAnnulus(ang, theta);

    if (intersects_annulus == 1) {
      // Fully contained in the annulus, so we get a copy of all of the pointers
      // in the sub-nodes.
      for (ITreePtrIterator iter=subpix_.begin();
	   iter!=subpix_.end();++iter) {
	_PointPtrs(i_ang);
      }
    } else {
      if (intersects_annulus == -1) {
	// Partial intersection with the annulus, so to pass the bounds along
	// to the sub-nodes.
	for (ITreePtrIterator iter=subpix_.begin();
	     iter!=subpix_.end();++iter) {
	  (*iter)->_FindPairs(ang, theta, i_ang);
	}
      }
    }
  }
}

void IndexedTreePixel::_PointPtrs(IAngularPtrVector& i_ang) {
  // Add any pointers contained in this node to the input vector or pass the
  // request along to the sub-nodes.
  if (!ang_.empty()) {
    for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      i_ang.push_back(*iter);
    }
  } else {
    for (ITreePtrIterator iter=subpix_.begin();
	 iter!=subpix_.end();++iter) {
      _PointPtrs(i_ang);
    }
  }
}

uint16_t IndexedTreePixel::FindNearestNeighbor(
  AngularCoordinate& ang, IndexedAngularCoordinate& neighbor_ang) {

  IAngularVector angVec;

  uint16_t nodes_visited = FindKNearestNeighbors(ang, 1, angVec);

  neighbor_ang = angVec[0];

  return nodes_visited;
}

double IndexedTreePixel::KNearestNeighborDistance(AngularCoordinate& ang,
						  uint8_t n_neighbors,
						  uint16_t& nodes_visited) {

  IndexedTreeNeighbor neighbors(ang, n_neighbors);

  _NeighborRecursion(ang, neighbors);

  nodes_visited = neighbors.NodesVisited();

  return neighbors.MaxAngularDistance();
}

double IndexedTreePixel::NearestNeighborDistance(AngularCoordinate& ang,
						 uint16_t& nodes_visited) {
  return KNearestNeighborDistance(ang, 1, nodes_visited);
}

bool IndexedTreePixel::ClosestMatch(AngularCoordinate& ang, double max_distance,
				    IndexedAngularCoordinate& match_ang) {
  IndexedTreeNeighbor neighbors(ang, 1, max_distance);

  _NeighborRecursion(ang, neighbors);

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

void IndexedTreePixel::_NeighborRecursion(AngularCoordinate& ang,
					  IndexedTreeNeighbor& neighbors) {

  neighbors.AddNode();

  if (!ang_.empty()) {
    // We have no sub-nodes in this tree, so we'll just iterate over the
    // points here and take the nearest N neighbors.
    for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter)
      neighbors.TestPoint(*iter);
  } else {
    // This node is the root node for our tree, so we first find the sub-node
    // that contains the point and start recursing there.
    //
    // While we iterate through the nodes, we'll also calculate the edge
    // distances for those nodes that don't contain the point and store them
    // in a priority queue.  This will let us do a follow-up check on nodes in
    // the most productive order.
    IPixelQueue pix_queue;
    for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
      if ((*iter)->Contains(ang)) {
	(*iter)->_NeighborRecursion(ang, neighbors);
      } else {
	double min_edge_distance, max_edge_distance;
	(*iter)->EdgeDistances(ang, min_edge_distance, max_edge_distance);
	DistanceIPixelPair dist_pair(min_edge_distance, (*iter));
	pix_queue.push(dist_pair);
      }
    }

    // That should give us back a TreeNeighbor object that contains a workable
    // set of neighbors and a search radius for possible matches.  Now we just
    // need to iterate over those sub-nodes that didn't contain the input point
    // to verify that there can't be any points in their sub-nodes which might
    // be closer to the input point.
    //
    // There's also the possibility that the input point is completely outside
    // our tree.  In that case (where the number of neighbors in the
    // TreeNeighbor object is less than the maximum), we want to check
    // all nodes.
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

void IndexedTreePixel::InitializeCorners() {
  unit_sphere_x_ = -1.0*sin(Lambda()*Stomp::DegToRad);
  unit_sphere_y_ = cos(Lambda()*Stomp::DegToRad)*
    cos(Eta()*Stomp::DegToRad+Stomp::EtaPole);
  unit_sphere_z_ = cos(Lambda()*Stomp::DegToRad)*
    sin(Eta()*Stomp::DegToRad+Stomp::EtaPole);

  unit_sphere_x_ul_ = -1.0*sin(LambdaMax()*Stomp::DegToRad);
  unit_sphere_y_ul_ = cos(LambdaMax()*Stomp::DegToRad)*
    cos(EtaMin()*Stomp::DegToRad+Stomp::EtaPole);
  unit_sphere_z_ul_ = cos(LambdaMax()*Stomp::DegToRad)*
    sin(EtaMin()*Stomp::DegToRad+Stomp::EtaPole);

  unit_sphere_x_ur_ = -1.0*sin(LambdaMax()*Stomp::DegToRad);
  unit_sphere_y_ur_ = cos(LambdaMax()*Stomp::DegToRad)*
    cos(EtaMax()*Stomp::DegToRad+Stomp::EtaPole);
  unit_sphere_z_ur_ = cos(LambdaMax()*Stomp::DegToRad)*
    sin(EtaMax()*Stomp::DegToRad+Stomp::EtaPole);

  unit_sphere_x_ll_ = -1.0*sin(LambdaMin()*Stomp::DegToRad);
  unit_sphere_y_ll_ = cos(LambdaMin()*Stomp::DegToRad)*
    cos(EtaMin()*Stomp::DegToRad+Stomp::EtaPole);
  unit_sphere_z_ll_ = cos(LambdaMin()*Stomp::DegToRad)*
    sin(EtaMin()*Stomp::DegToRad+Stomp::EtaPole);

  unit_sphere_x_lr_ = -1.0*sin(LambdaMin()*Stomp::DegToRad);
  unit_sphere_y_lr_ = cos(LambdaMin()*Stomp::DegToRad)*
    cos(EtaMax()*Stomp::DegToRad+Stomp::EtaPole);
  unit_sphere_z_lr_ = cos(LambdaMin()*Stomp::DegToRad)*
    sin(EtaMax()*Stomp::DegToRad+Stomp::EtaPole);
}

bool IndexedTreePixel::AddPoint(IndexedAngularCoordinate* ang) {
  bool added_to_pixel = false;
  if (Contains(*ang)) {
    if ((point_count_ < maximum_points_) ||
	(Resolution() == Stomp::MaxPixelResolution)) {
      if (point_count_ == 0) ang_.reserve(maximum_points_);
      ang_.push_back(ang);
      added_to_pixel = true;
    } else {
      if (!initialized_subpixels_) {
	if (!_InitializeSubPixels()) {
	  std::cout << "Failed to initialize sub-pixels.  Exiting.\n";
	  exit(2);
	}
      }
      for (uint32_t i=0;i<subpix_.size();++i) {
	if (subpix_[i]->Contains(*ang)) {
	  added_to_pixel = subpix_[i]->AddPoint(ang);
	  i = subpix_.size();
	}
      }
    }
  } else {
    added_to_pixel = false;
  }

  if (added_to_pixel) point_count_++;

  return added_to_pixel;
}

bool IndexedTreePixel::AddPoint(IndexedAngularCoordinate& i_ang) {
  IndexedAngularCoordinate* ang_copy =
    new IndexedAngularCoordinate(i_ang.UnitSphereX(), i_ang.UnitSphereY(),
				 i_ang.UnitSphereZ(), i_ang.Index());
  return AddPoint(ang_copy);
}

bool IndexedTreePixel::AddPoint(AngularCoordinate& ang, uint32_t index) {
  IndexedAngularCoordinate* i_ang =
    new IndexedAngularCoordinate(ang.UnitSphereX(), ang.UnitSphereY(),
				 ang.UnitSphereZ(), index);
  return AddPoint(i_ang);
}

uint32_t IndexedTreePixel::NPoints() {
  return point_count_;
}

uint32_t IndexedTreePixel::NPoints(Pixel& pix) {
  uint32_t total_points = 0;

  // First check to see if the input pixel contains the current pixel.
  if (pix.Contains(Resolution(), PixelX(), PixelY())) {
    // If so, then it also contains all of the points in the current pixel.
    total_points = NPoints();
  } else {
    // If not, then either the input pixel doesn't overlap the current one or
    // it's a sub-pixel of the current one.
    if (Contains(pix)) {
      // If we contain the input pixel, then we either iterate over the
      // sub-nodes to this pixel or iterate over the points contained in this
      // pixel.
      if (initialized_subpixels_) {
	for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	  total_points += (*iter)->NPoints(pix);
	}
      } else {
	for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	  if (pix.Contains(*(*iter))) total_points++;
	}
      }
    }
  }

  return total_points;
}

void IndexedTreePixel::Indices(Pixel& pix, IndexVector& indices) {
  if (!indices.empty()) indices.clear();
  IAngularPtrVector i_ang;

  // First check to see if the input pixel contains the current pixel.
  if (pix.Contains(Resolution(), PixelX(), PixelY())) {
    // If so, then it also contains all of the points in the current pixel.
    _PointPtrs(i_ang);
  } else {
    // If not, then either the input pixel doesn't overlap the current one or
    // it's a sub-pixel of the current one.
    if (Contains(pix)) {
      // If we contain the input pixel, then we either iterate over the
      // sub-nodes to this pixel or iterate over the points contained in this
      // pixel.
      if (initialized_subpixels_) {
	for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	  (*iter)->Indices(pix, indices);
	}
      } else {
	for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	  if (pix.Contains(*(*iter))) i_ang.push_back(*iter);
	}
      }
    }
  }

  if (!i_ang.empty()) {
    for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      indices.push_back((*iter)->Index());
    }
  }
}

double IndexedTreePixel::Coverage() {
  double total_coverage = 0.0;

  // First check to see if the current pixel contains any sub-pixels.
  if (!initialized_subpixels_) {
    // If there are no sub-pixels, then we have no further information and
    // must assume that the whole pixel is covered with data, provided that
    // there's at least one point here.
    if (point_count_ > 0) total_coverage = 1.0;
  } else {
    // If we have sub-pixels, then we want to recursively probe the tree
    // structure to find out how many of the sub-pixels contain data.  Since
    // each sub-pixel contributes 1/4 the area of the current pixel, we
    // scale their results accordingly.
    for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
      total_coverage += 0.25*(*iter)->Coverage();
    }
  }

  return total_coverage;
}

double IndexedTreePixel::Coverage(Pixel& pix) {
  double total_coverage = 0.0;

  // First check to see if the input pixel contains the current pixel.
  if (pix.Contains(Resolution(), PixelX(), PixelY())) {
    // If so, then we return the coverage for this pixel, normalized by the
    // relative areas between the two pixels.
    total_coverage = Coverage()*Area()/pix.Area();
  } else {
    // If the input pixel doesn't contain the current pixel, then we need
    // to verify that the converse is true.  Otherwise, the Coverage() is 0.
    if (Contains(pix)) {
      // If there are no sub-pixels, then all we can say is that the input
      // pixel is completely covered by the current pixel.
      if (!initialized_subpixels_) {
	if (point_count_ > 0) total_coverage = 1.0;
      } else {
	// If we have sub-pixels, then we want to find the one that contains
	// the input pixel and recurse down the tree until we find either the
	// pixel itself or the last node that contains it.
	for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	  if ((*iter)->Contains(pix))
	    total_coverage = (*iter)->Coverage(pix);
	}
      }
    }
  }

  return total_coverage;
}

void IndexedTreePixel::Points(IAngularVector& i_ang) {
  if (!i_ang.empty()) i_ang.clear();
  i_ang.reserve(point_count_);

  // If we haven't initialized any sub-nodes, then this is just a matter of
  // creating a copy of all of the points in the current pixel.
  if (!initialized_subpixels_) {
    for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      IndexedAngularCoordinate tmp_ang = *(*iter);

      i_ang.push_back(tmp_ang);
    }
  } else {
    // If not, then we need to iterate through our sub-nodes and return an
    // aggregate list.
    for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
      IAngularVector tmp_ang;
      (*iter)->Points(tmp_ang);
      for (IAngularIterator ang_iter=tmp_ang.begin();
	   ang_iter!=tmp_ang.end();++ang_iter) i_ang.push_back(*ang_iter);
    }
  }
}

void IndexedTreePixel::Points(IAngularVector& i_ang, Pixel& pix) {
  if (!i_ang.empty()) i_ang.clear();
  i_ang.reserve(point_count_);

  // First, we need to check to verify that the input pixel is either contained
  // in the current pixel or contains it.
  if (Contains(pix) || pix.Contains(Resolution(), PixelX(), PixelY())) {
    // If we haven't initialized any sub-nodes, then this is just a matter of
    // creating a copy of all of the points in the current pixel that are
    // contained in the input pixel.
    if (!initialized_subpixels_) {
      for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	if (pix.Contains(*(*iter))) {
	  IndexedAngularCoordinate tmp_ang((*iter)->UnitSphereX(),
					   (*iter)->UnitSphereY(),
					   (*iter)->UnitSphereZ(),
					   (*iter)->Index());
	  i_ang.push_back(tmp_ang);
	}
      }
    } else {
      // If not, then we need to iterate through our sub-nodes and return an
      // aggregate list.
      for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	IAngularVector tmp_ang;
	(*iter)->Points(tmp_ang, pix);
	for (IAngularIterator ang_iter=tmp_ang.begin();
	     ang_iter!=tmp_ang.end();++ang_iter) i_ang.push_back(*ang_iter);
      }
    }
  }
}

uint16_t IndexedTreePixel::Nodes() {
  uint16_t n_nodes = 1;

  for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter)
    (*iter)->_AddSubNodes(n_nodes);

  return n_nodes;
}

void IndexedTreePixel::_AddSubNodes(uint16_t& n_nodes) {
  n_nodes++;

  for (ITreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter)
    (*iter)->_AddSubNodes(n_nodes);
}

void IndexedTreePixel::SetPixelCapacity(uint16_t maximum_points) {
  maximum_points_ = maximum_points;
}

uint16_t IndexedTreePixel::PixelCapacity() {
  return maximum_points_;
}

IAngularPtrIterator IndexedTreePixel::PointsBegin() {
  return ang_.begin();
}

IAngularPtrIterator IndexedTreePixel::PointsEnd() {
  return ang_.end();
}

ITreePtrIterator IndexedTreePixel::NodesBegin() {
  return subpix_.begin();
}

ITreePtrIterator IndexedTreePixel::NodesEnd() {
  return subpix_.end();
}

bool IndexedTreePixel::HasPoints() {
  return (ang_.empty() ? false : true);
}

bool IndexedTreePixel::HasNodes() {
  return (subpix_.empty() ? false : true);
}

void IndexedTreePixel::Clear() {
  if (!ang_.empty())
    for (IAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter)
      delete *iter;
  ang_.clear();
  if (!subpix_.empty())
    for (uint32_t i=0;i<subpix_.size();i++) {
      subpix_[i]->Clear();
      delete subpix_[i];
    }
  subpix_.clear();
}

double IndexedTreePixel::UnitSphereX() {
  return unit_sphere_x_;
}

double IndexedTreePixel::UnitSphereY() {
  return unit_sphere_y_;
}

double IndexedTreePixel::UnitSphereZ() {
  return unit_sphere_z_;
}

double IndexedTreePixel::UnitSphereX_UL() {
  return unit_sphere_x_ul_;
}

double IndexedTreePixel::UnitSphereY_UL() {
  return unit_sphere_y_ul_;
}

double IndexedTreePixel::UnitSphereZ_UL() {
  return unit_sphere_z_ul_;
}

double IndexedTreePixel::UnitSphereX_UR() {
  return unit_sphere_x_ur_;
}

double IndexedTreePixel::UnitSphereY_UR() {
  return unit_sphere_y_ur_;
}

double IndexedTreePixel::UnitSphereZ_UR() {
  return unit_sphere_z_ur_;
}

double IndexedTreePixel::UnitSphereX_LL() {
  return unit_sphere_x_ll_;
}

double IndexedTreePixel::UnitSphereY_LL() {
  return unit_sphere_y_ll_;
}

double IndexedTreePixel::UnitSphereZ_LL() {
  return unit_sphere_z_ll_;
}

double IndexedTreePixel::UnitSphereX_LR() {
  return unit_sphere_x_lr_;
}

double IndexedTreePixel::UnitSphereY_LR() {
  return unit_sphere_y_lr_;
}

double IndexedTreePixel::UnitSphereZ_LR() {
  return unit_sphere_z_lr_;
}

void IndexedTreePixel::WithinAnnulus(AngularBin& theta, PixelVector& pix,
				     bool check_full_pixel) {
  if (!pix.empty()) pix.clear();

  uint32_t y_min;
  uint32_t y_max;
  std::vector<uint32_t> x_min;
  std::vector<uint32_t> x_max;

  XYBounds(theta.ThetaMax(), x_min, x_max, y_min, y_max, true);

  uint32_t nx = Nx0*Resolution();
  uint32_t nx_pix;

  for (uint32_t y=y_min,n=0;y<=y_max;y++,n++) {
    if ((x_max[n] < x_min[n]) && (x_min[n] > nx/2)) {
      nx_pix = nx - x_min[n] + x_max[n] + 1;
    } else {
      nx_pix = x_max[n] - x_min[n] + 1;
    }
    if (nx_pix > nx) nx_pix = nx;
    for (uint32_t m=0,x=x_min[n];m<nx_pix;m++,x++) {
      if (x == nx) x = 0;
      IndexedTreePixel tree_pix(x, y, Resolution());
      bool within_bounds =
	theta.WithinCosBounds(UnitSphereX()*tree_pix.UnitSphereX() +
			      UnitSphereY()*tree_pix.UnitSphereY() +
			      UnitSphereZ()*tree_pix.UnitSphereZ());
      if (check_full_pixel && within_bounds) {
	if (theta.WithinCosBounds(UnitSphereX()*tree_pix.UnitSphereX_UL() +
				  UnitSphereY()*tree_pix.UnitSphereY_UL() +
				  UnitSphereZ()*tree_pix.UnitSphereZ_UL()) &&
	    theta.WithinCosBounds(UnitSphereX()*tree_pix.UnitSphereX_UR() +
				  UnitSphereY()*tree_pix.UnitSphereY_UR() +
				  UnitSphereZ()*tree_pix.UnitSphereZ_UR()) &&
	    theta.WithinCosBounds(UnitSphereX()*tree_pix.UnitSphereX_LL() +
				  UnitSphereY()*tree_pix.UnitSphereY_LL() +
				  UnitSphereZ()*tree_pix.UnitSphereZ_LL()) &&
	    theta.WithinCosBounds(UnitSphereX()*tree_pix.UnitSphereX_LR() +
				  UnitSphereY()*tree_pix.UnitSphereY_LR() +
				  UnitSphereZ()*tree_pix.UnitSphereZ_LR())) {
	  within_bounds = true;
	} else {
	  within_bounds = false;
	}
      }
      if (within_bounds) pix.push_back(Pixel(tree_pix.PixelX(),
					     tree_pix.PixelY(),
					     tree_pix.Resolution(), 1.0));
    }
  }
}

IndexedTreeNeighbor::IndexedTreeNeighbor(AngularCoordinate& reference_ang,
					 uint8_t n_neighbor) {
  reference_ang_ = reference_ang;
  n_neighbors_ = n_neighbor;
  max_distance_ = 100.0;
  n_nodes_visited_ = 0;
}

IndexedTreeNeighbor::IndexedTreeNeighbor(AngularCoordinate& reference_ang,
					 uint8_t n_neighbor,
					 double max_distance) {
  reference_ang_ = reference_ang;
  n_neighbors_ = n_neighbor;
  max_distance_ = sin(DegToRad*max_distance)*sin(DegToRad*max_distance);
  n_nodes_visited_ = 0;
}

IndexedTreeNeighbor::~IndexedTreeNeighbor() {
  n_neighbors_ = 0;
  max_distance_ = 100.0;
  n_nodes_visited_ = 0;
}

void IndexedTreeNeighbor::NearestNeighbors(IAngularVector& i_ang,
					   bool save_neighbors) {
  if (!i_ang.empty()) i_ang.clear();

  std::vector<DistanceIPointPair> backup_copy;

  while (!ang_queue_.empty()) {
    DistanceIPointPair dist_pair = ang_queue_.top();
    ang_queue_.pop();

    IndexedAngularCoordinate tmp_ang(dist_pair.second->UnitSphereX(),
				     dist_pair.second->UnitSphereY(),
				     dist_pair.second->UnitSphereZ(),
				     dist_pair.second->Index());

    i_ang.push_back(tmp_ang);
    backup_copy.push_back(dist_pair);
  }

  if (save_neighbors) {
    for (uint8_t i=0;i<backup_copy.size();i++) {
      ang_queue_.push(backup_copy[i]);
    }
  }
}

uint8_t IndexedTreeNeighbor::Neighbors() {
  return ang_queue_.size();
}

uint8_t IndexedTreeNeighbor::MaxNeighbors() {
  return n_neighbors_;
}

bool IndexedTreeNeighbor::TestPoint(IndexedAngularCoordinate* test_ang) {
  bool kept_point = false;

  double costheta = reference_ang_.DotProduct(test_ang);

  double sin2theta = 1.0 - costheta*costheta;

  if (sin2theta < max_distance_ || Neighbors() < MaxNeighbors()) {
    // If the new point is closer than the most distant point in our queue or
    // we're still filling our heap, then we keep it.
    kept_point = true;

    // Throw away the current most distant point if we're at capacity.
    if (Neighbors() == MaxNeighbors()) ang_queue_.pop();

    // Create a new pair for the test point and add it to the queue.
    DistanceIPointPair dist_pair(sin2theta, test_ang);
    ang_queue_.push(dist_pair);

    // And reset our maximum distance using the new top of the heap.
    max_distance_ = ang_queue_.top().first;
  }

  return kept_point;
}

double IndexedTreeNeighbor::MaxDistance() {
  return max_distance_;
}

double IndexedTreeNeighbor::MaxAngularDistance() {
  // Numerical precision can sometimes make the max_distance_ negative.
  return RadToDeg*asin(sqrt(fabs(max_distance_)));
}

uint16_t IndexedTreeNeighbor::NodesVisited() {
  return n_nodes_visited_;
}

void IndexedTreeNeighbor::AddNode() {
  n_nodes_visited_++;
}

} // end namespace Stomp
