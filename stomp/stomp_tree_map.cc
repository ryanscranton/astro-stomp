// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the TreeMap class.  The core work of pair finding
// and K nearest neighbor searches is done in the TreePixel class.  However,
// due to the pixelization scheme used in STOMP, there is a maximum pixel size
// that does not span the entire sphere.  Hence, a vector of TreePixels is
// necessary to describe an arbitrary collection of points.  TreeMap manages
// that vector of TreePixels, adding them as necessary based on the input
// points.

#include "stomp_core.h"
#include "stomp_tree_map.h"
#include "stomp_map.h"
#include "stomp_angular_bin.h"
#include "stomp_radial_bin.h"
#include "stomp_angular_correlation.h"
#include "stomp_util.h"

namespace Stomp {

TreeMap::TreeMap(uint32_t input_resolution, uint16_t maximum_points) {
  resolution_ = input_resolution;
  maximum_points_ = maximum_points;
  weight_ = 0.0;
  point_count_ = 0;
  modified_ = false;
  area_ = 0.0;
  ClearRegions();
}

TreeMap::TreeMap(const std::string& input_file, uint32_t input_resolution,
		 uint16_t maximum_points, AngularCoordinate::Sphere sphere,
		 bool verbose, uint8_t theta_column, uint8_t phi_column,
		 int8_t weight_column) {
  resolution_ = input_resolution;
  maximum_points_ = maximum_points;
  weight_ = 0.0;
  point_count_ = 0;
  modified_ = false;
  area_ = 0.0;
  ClearRegions();

  Read(input_file, sphere, verbose, theta_column, phi_column, weight_column);
}

TreeMap::TreeMap(const std::string& input_file, FieldColumnDict& field_columns,
		 uint32_t input_resolution, uint16_t maximum_points,
		 AngularCoordinate::Sphere sphere, bool verbose,
		 uint8_t theta_column, uint8_t phi_column,
		 int8_t weight_column) {
  resolution_ = input_resolution;
  maximum_points_ = maximum_points;
  weight_ = 0.0;
  point_count_ = 0;
  modified_ = false;
  area_ = 0.0;
  ClearRegions();

  Read(input_file, field_columns, sphere, verbose,
       theta_column, phi_column, weight_column);
}

TreeMap::~TreeMap() {
  Clear();
  resolution_ = 0;
  maximum_points_ = 0;
  weight_ = 0.0;
  point_count_ = 0;
  modified_ = false;
  area_ = 0.0;
  ClearRegions();
}

uint32_t TreeMap::FindPairs(AngularCoordinate& ang, AngularBin& theta) {
  uint32_t pair_count = 0;

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  center_pix.BoundingRadius(ang, theta.ThetaMax(), pix);

  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end())
      pair_count += tree_map_[pix_iter->Pixnum()]->FindPairs(ang, theta);
  }
  return pair_count;
}

uint32_t TreeMap::FindPairs(AngularCoordinate& ang,
				 double theta_min, double theta_max) {
  AngularBin theta(theta_min, theta_max);
  return FindPairs(ang, theta);
}

uint32_t TreeMap::FindPairs(AngularCoordinate& ang, double theta_max) {
  AngularBin theta(0.0, theta_max);
  return FindPairs(ang, theta);
}

void TreeMap::FindPairs(AngularVector& ang, AngularBin& theta) {
  uint32_t n_pairs = 0;

  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    n_pairs = FindPairs(*ang_iter, theta);
  }
}

void TreeMap::FindPairs(AngularVector& ang, AngularCorrelation& wtheta) {
  uint32_t n_pairs = 0;

  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (AngularIterator ang_iter=ang.begin();
	 ang_iter!=ang.end();++ang_iter) {
      n_pairs = FindPairs(*ang_iter, *theta_iter);
    }
  }
}

double TreeMap::FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta) {
  double total_weight = 0.0;

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  center_pix.BoundingRadius(ang, theta.ThetaMax(), pix);

  // Now we iterate through the possibilities and add their contributions to
  // the total.
  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end())
      total_weight +=
	tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(ang, theta);
  }
  return total_weight;
}

double TreeMap::FindWeightedPairs(AngularCoordinate& ang,
				  double theta_min, double theta_max) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(ang, theta);
}

double TreeMap::FindWeightedPairs(AngularCoordinate& ang, double theta_max) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(ang, theta);
}

void TreeMap::FindWeightedPairs(AngularVector& ang, AngularBin& theta) {
  double total_weight = 0.0;

  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, theta);
  }
}

void TreeMap::FindWeightedPairs(AngularVector& ang,
				AngularCorrelation& wtheta) {
  double total_weight = 0.0;

  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (AngularIterator ang_iter=ang.begin();
	 ang_iter!=ang.end();++ang_iter) {
      total_weight = FindWeightedPairs(*ang_iter, *theta_iter);
    }
  }
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  AngularBin& theta) {
  double total_weight = 0.0;

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  center_pix.BoundingRadius(w_ang, theta.ThetaMax(), pix);

  // Now we iterate through the possibilities and add their contributions to
  // the total.
  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end())
      total_weight +=
	tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(w_ang, theta);
  }
  return total_weight;
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  double theta_min, double theta_max) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(w_ang, theta);
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  double theta_max) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(w_ang, theta);
}

void TreeMap::FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta) {
  double total_weight = 0.0;

  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, theta);
  }
}

void TreeMap::FindWeightedPairs(WAngularVector& w_ang,
				AngularCorrelation& wtheta) {
  double total_weight = 0.0;

  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (WAngularIterator ang_iter=w_ang.begin();
	 ang_iter!=w_ang.end();++ang_iter) {
      total_weight = FindWeightedPairs(*ang_iter, *theta_iter);
    }
  }
}

double TreeMap::FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
				  const std::string& field_name) {
  double total_weight = 0.0;

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  center_pix.BoundingRadius(ang, theta.ThetaMax(), pix);

  // Now we iterate through the possibilities and add their contributions to
  // the total.
  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end())
      total_weight +=
	tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(ang, theta,
							 field_name);
  }
  return total_weight;
}

double TreeMap::FindWeightedPairs(AngularCoordinate& ang,
				  double theta_min, double theta_max,
				  const std::string& field_name) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(ang, theta, field_name);
}

double TreeMap::FindWeightedPairs(AngularCoordinate& ang, double theta_max,
				  const std::string& field_name) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(ang, theta, field_name);
}

void TreeMap::FindWeightedPairs(AngularVector& ang, AngularBin& theta,
				const std::string& field_name) {
  double total_weight = 0.0;

  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, theta, field_name);
  }
}

void TreeMap::FindWeightedPairs(AngularVector& ang,
				AngularCorrelation& wtheta,
				const std::string& field_name) {
  double total_weight = 0.0;

  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (AngularIterator ang_iter=ang.begin();
	 ang_iter!=ang.end();++ang_iter) {
      total_weight = FindWeightedPairs(*ang_iter, *theta_iter, field_name);
    }
  }
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  AngularBin& theta,
				  const std::string& field_name) {
  double total_weight = 0.0;

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  center_pix.BoundingRadius(w_ang, theta.ThetaMax(), pix);

  // Now we iterate through the possibilities and add their contributions to
  // the total.
  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end())
      total_weight +=
	tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(w_ang, theta,
							 field_name);
  }
  return total_weight;
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  double theta_min, double theta_max,
				  const std::string& field_name) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(w_ang, theta, field_name);
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  double theta_max,
				  const std::string& field_name) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(w_ang, theta, field_name);
}

void TreeMap::FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta,
				const std::string& field_name) {
  double total_weight = 0.0;

  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, theta, field_name);
  }
}

void TreeMap::FindWeightedPairs(WAngularVector& w_ang,
				AngularCorrelation& wtheta,
				const std::string& field_name) {
  double total_weight = 0.0;

  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (WAngularIterator ang_iter=w_ang.begin();
	 ang_iter!=w_ang.end();++ang_iter) {
      total_weight = FindWeightedPairs(*ang_iter, *theta_iter, field_name);
    }
  }
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  const std::string& ang_field_name,
				  AngularBin& theta,
				  const std::string& field_name) {
  double total_weight = 0.0;

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  center_pix.BoundingRadius(w_ang, theta.ThetaMax(), pix);

  // Now we iterate through the possibilities and add their contributions to
  // the total.
  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end())
      total_weight +=
	tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(w_ang, ang_field_name,
							 theta, field_name);
  }
  return total_weight;
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  const std::string& ang_field_name,
				  double theta_min, double theta_max,
				  const std::string& field_name) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(w_ang, ang_field_name, theta, field_name);
}

double TreeMap::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				  const std::string& ang_field_name,
				  double theta_max,
				  const std::string& field_name) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(w_ang, ang_field_name, theta, field_name);
}

void TreeMap::FindWeightedPairs(WAngularVector& w_ang,
				const std::string& ang_field_name,
				AngularBin& theta,
				const std::string& field_name) {
  double total_weight = 0.0;

  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, ang_field_name,
				     theta, field_name);
  }
}

void TreeMap::FindWeightedPairs(WAngularVector& w_ang,
				const std::string& ang_field_name,
				AngularCorrelation& wtheta,
				const std::string& field_name) {
  double total_weight = 0.0;

  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (WAngularIterator ang_iter=w_ang.begin();
	 ang_iter!=w_ang.end();++ang_iter) {
      total_weight = FindWeightedPairs(*ang_iter, ang_field_name,
				       *theta_iter, field_name);
    }
  }
}

void TreeMap::FindPairsWithRegions(AngularVector& ang, AngularBin& theta) {
  if (!RegionsInitialized()) {
    std::cout <<
      "Stomp::TreeMap::FindPairsWithRegions - " <<
      "Must initialize regions before calling FindPairsWithRegions\n" <<
      "\tExiting...\n";
    exit(2);
  }

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  uint32_t n_pair = 0;
  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    center_pix.BoundingRadius(*ang_iter, theta.ThetaMax(), pix);
    uint16_t region = FindRegion(center_pix);

    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      if (region == FindRegion(*pix_iter)) {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  n_pair = tree_map_[pix_iter->Pixnum()]->FindPairs(*ang_iter,
							    theta, region);
      } else {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  n_pair = tree_map_[pix_iter->Pixnum()]->FindPairs(*ang_iter, theta);
      }
    }
  }
}

void TreeMap::FindPairsWithRegions(AngularVector& ang,
				   AngularCorrelation& wtheta) {
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter)
    FindPairsWithRegions(ang, *theta_iter);
}

void TreeMap::FindWeightedPairsWithRegions(AngularVector& ang,
					   AngularBin& theta) {
  if (!RegionsInitialized()) {
    std::cout <<
      "Stomp::TreeMap::FindWeightedPairsWithRegions - " <<
      "Must initialize regions before calling FindPairsWithRegions\n" <<
      "\tExiting...\n";
    exit(2);
  }

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  double total_weight = 0.0;
  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    center_pix.BoundingRadius(*ang_iter, theta.ThetaMax(), pix);
    uint16_t region = FindRegion(center_pix);

    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      if (region == FindRegion(*pix_iter)) {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter,
							     theta, region);
      } else {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter, theta);
      }
    }
  }
}

void TreeMap::FindWeightedPairsWithRegions(AngularVector& ang,
					   AngularCorrelation& wtheta) {
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter)
    FindWeightedPairsWithRegions(ang, *theta_iter);
}

void TreeMap::FindWeightedPairsWithRegions(WAngularVector& w_ang,
					   AngularBin& theta) {
  if (!RegionsInitialized()) {
    std::cout <<
      "Stomp::TreeMap::FindWeightedPairsWithRegions - " <<
      "Must initialize regions before calling FindPairsWithRegions\n" <<
      "\tExiting...\n";
    exit(2);
  }

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  double total_weight = 0.0;
  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    center_pix.BoundingRadius(*ang_iter, theta.ThetaMax(), pix);
    uint16_t region = FindRegion(center_pix);

    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      if (region == FindRegion(*pix_iter)) {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter,
							     theta, region);
      } else {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter, theta);
      }
    }
  }
}

void TreeMap::FindWeightedPairsWithRegions(CosmoVector& c_ang,
					   RadialBin& radius) {
  if (!RegionsInitialized()) {
    std::cout <<
      "Stomp::TreeMap::FindWeightedPairsWithRegions - " <<
      "Must initialize regions before calling FindPairsWithRegions\n" <<
      "\tExiting...\n";
    exit(2);
  }

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  double total_weight = 0.0;
  for (CosmoIterator ang_iter=c_ang.begin();
       ang_iter!=c_ang.end();++ang_iter) {
    radius.SetRedshift(ang_iter->Redshift());
    center_pix.BoundingRadius(*ang_iter, radius.ThetaMax(), pix);
    uint16_t region = FindRegion(center_pix);

    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      if (region == FindRegion(*pix_iter)) {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter,
							     radius, region);
      } else {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter, radius);
      }
    }
  }
}

void TreeMap::FindWeightedPairsWithRegions(WAngularVector& w_ang,
					   AngularCorrelation& wtheta) {
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter)
    FindWeightedPairsWithRegions(w_ang, *theta_iter);
}

void TreeMap::FindWeightedPairsWithRegions(AngularVector& ang,
					   AngularBin& theta,
					   const std::string& field_name) {
  if (!RegionsInitialized()) {
    std::cout <<
      "Stomp::TreeMap::FindWeightedPairsWithRegions - " <<
      "Must initialize regions before calling FindPairsWithRegions\n" <<
      "\tExiting...\n";
    exit(2);
  }

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  double total_weight = 0.0;
  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    center_pix.BoundingRadius(*ang_iter, theta.ThetaMax(), pix);
    uint16_t region = FindRegion(center_pix);

    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      if (region == FindRegion(*pix_iter)) {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter, theta,
							     field_name,region);
      } else {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter, theta,
							     field_name);
      }
    }
  }
}

void TreeMap::FindWeightedPairsWithRegions(AngularVector& ang,
					   AngularCorrelation& wtheta,
					   const std::string& field_name) {
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter)
    FindWeightedPairsWithRegions(ang, *theta_iter, field_name);
}

void TreeMap::FindWeightedPairsWithRegions(WAngularVector& w_ang,
					   AngularBin& theta,
					   const std::string& field_name) {
  if (!RegionsInitialized()) {
    std::cout <<
      "Stomp::TreeMap::FindWeightedPairsWithRegions - " <<
      "Must initialize regions before calling FindPairsWithRegions\n" <<
      "\tExiting...\n";
    exit(2);
  }

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  double total_weight = 0.0;
  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    center_pix.BoundingRadius(*ang_iter, theta.ThetaMax(), pix);
    uint16_t region = FindRegion(center_pix);

    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      if (region == FindRegion(*pix_iter)) {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter, theta,
							     field_name,region);
      } else {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter, theta,
							     field_name);
      }
    }
  }
}

void TreeMap::FindWeightedPairsWithRegions(WAngularVector& w_ang,
					   AngularCorrelation& wtheta,
					   const std::string& field_name) {
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter)
    FindWeightedPairsWithRegions(w_ang, *theta_iter, field_name);
}

void TreeMap::FindWeightedPairsWithRegions(WAngularVector& w_ang,
					   const std::string& ang_field_name,
					   AngularBin& theta,
					   const std::string& field_name) {
  if (!RegionsInitialized()) {
    std::cout <<
      "Stomp::TreeMap::FindWeightedPairsWithRegions - " <<
      "Must initialize regions before calling FindPairsWithRegions\n" <<
      "\tExiting...\n";
    exit(2);
  }

  // First we need to find out which pixels this angular bin possibly touches.
  Pixel center_pix;
  center_pix.SetResolution(resolution_);
  PixelVector pix;
  double total_weight = 0.0;
  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    center_pix.BoundingRadius(*ang_iter, theta.ThetaMax(), pix);
    uint16_t region = FindRegion(center_pix);

    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      if (region == FindRegion(*pix_iter)) {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter,
							     ang_field_name,
							     theta, field_name,
							     region);
      } else {
	TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
	if (iter != tree_map_.end())
	  total_weight =
	    tree_map_[pix_iter->Pixnum()]->FindWeightedPairs(*ang_iter,
							     ang_field_name,
							     theta, field_name);
      }
    }
  }
}

void TreeMap::FindWeightedPairsWithRegions(WAngularVector& w_ang,
					   const std::string& ang_field_name,
					   AngularCorrelation& wtheta,
					   const std::string& field_name) {
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter)
    FindWeightedPairsWithRegions(w_ang, ang_field_name,
				 *theta_iter, field_name);
}

uint16_t TreeMap::FindKNearestNeighbors(AngularCoordinate& ang,
					uint8_t n_neighbors,
					WAngularVector& neighbor_ang) {
  TreeNeighbor neighbors(ang, n_neighbors);

  _NeighborRecursion(ang, neighbors);

  neighbors.NearestNeighbors(neighbor_ang, false);

  return neighbors.NodesVisited();
}

uint16_t TreeMap::FindNearestNeighbor(AngularCoordinate& ang,
				    WeightedAngularCoordinate& neighbor_ang) {
  WAngularVector angVec;

  uint16_t nodes_visited = FindKNearestNeighbors(ang, 1, angVec);

  neighbor_ang = angVec[0];

  return nodes_visited;
}

double TreeMap::KNearestNeighborDistance(AngularCoordinate& ang,
					 uint8_t n_neighbors,
					 uint16_t& nodes_visited) {

  TreeNeighbor neighbors(ang, n_neighbors);

  _NeighborRecursion(ang, neighbors);

  nodes_visited = neighbors.NodesVisited();

  return neighbors.MaxAngularDistance();
}

double TreeMap::NearestNeighborDistance(AngularCoordinate& ang,
					uint16_t& nodes_visited) {
  return KNearestNeighborDistance(ang, 1, nodes_visited);
}

bool TreeMap::ClosestMatch(AngularCoordinate& ang,
			   double max_distance,
			   WeightedAngularCoordinate& match_ang) {
  TreeNeighbor neighbors(ang, 1, max_distance);

  _MatchRecursion(ang, neighbors);

  bool found_match = false;
  if (neighbors.Neighbors() == neighbors.MaxNeighbors() &&
      neighbors.MaxAngularDistance() < max_distance) {
    found_match = true;

    WAngularVector neighbor_ang;
    neighbors.NearestNeighbors(neighbor_ang, false);
    match_ang = neighbor_ang[0];
  }

  return found_match;
}

void TreeMap::_NeighborRecursion(AngularCoordinate& ang,
				 TreeNeighbor& neighbors) {

  // First we need to find out if the input point is within our map area.
  Pixel center_pix(ang, resolution_);
  TreeDictIterator iter = tree_map_.find(center_pix.Pixnum());

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
  PixelQueue pix_queue;
  for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
    TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
    if (iter != tree_map_.end() && !pix_iter->Contains(ang)) {
      double min_edge_distance, max_edge_distance;
      tree_map_[pix_iter->Pixnum()]->EdgeDistances(ang, min_edge_distance,
						   max_edge_distance);
      DistancePixelPair dist_pair(min_edge_distance,
				  tree_map_[pix_iter->Pixnum()]);
      pix_queue.push(dist_pair);
    }
  }

  // And iterate over that queue to check for neighbors.
  while (!pix_queue.empty()) {
    double pix_distance = pix_queue.top().first;
    TreePixel* pix_iter = pix_queue.top().second;
    if (pix_distance < neighbors.MaxDistance()) {
      pix_iter->_NeighborRecursion(ang, neighbors);
    }
    pix_queue.pop();
  }
}

void TreeMap::_MatchRecursion(AngularCoordinate& ang,
				 TreeNeighbor& neighbors) {

  // First we need to find out if the input point is within our map area.
  Pixel center_pix(ang, resolution_);
  TreeDictIterator iter = tree_map_.find(center_pix.Pixnum());

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
    PixelQueue pix_queue;
    for (PixelIterator pix_iter=pix.begin();pix_iter!=pix.end();++pix_iter) {
      TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end() && !pix_iter->Contains(ang)) {
	tree_map_[pix_iter->Pixnum()]->EdgeDistances(ang, min_edge_distance,
						     max_edge_distance);
	DistancePixelPair dist_pair(min_edge_distance,
				    tree_map_[pix_iter->Pixnum()]);
	pix_queue.push(dist_pair);
      }
    }

    // And iterate over that queue to check for neighbors.
    while (!pix_queue.empty()) {
      double pix_distance = pix_queue.top().first;
      TreePixel* pix_iter = pix_queue.top().second;
      if (pix_distance < neighbors.MaxDistance()) {
      pix_iter->_NeighborRecursion(ang, neighbors);
      }
      pix_queue.pop();
    }
  }
}

bool TreeMap::AddPoint(WeightedAngularCoordinate* ang) {
  Pixel pix;
  pix.SetResolution(resolution_);
  pix.SetPixnumFromAng(*ang);

  TreeDictIterator iter = tree_map_.find(pix.Pixnum());
  if (iter == tree_map_.end()) {
    // If we didn't find the pixnum key in the map, then we need to add this
    // pixnum to the map and re-do the search.
    tree_map_.insert(std::pair<uint32_t,
		     TreePixel *>(pix.Pixnum(),
				  new TreePixel(pix.PixelX(), pix.PixelY(),
						resolution_, maximum_points_)));
    iter = tree_map_.find(pix.Pixnum());
    if (iter == tree_map_.end()) {
      std::cout << "Stomp::TreeMap::AddPoint - " <<
	"Creating new TreeMap node failed. Exiting.\n";
      exit(2);
    }
  }
  bool added_point = (*iter).second->AddPoint(ang);
  if (added_point) {
    point_count_++;
    weight_ += ang->Weight();
    if (ang->HasFields()) {
      for (FieldIterator iter=ang->FieldBegin();iter!=ang->FieldEnd();++iter) {
	if (field_total_.find(iter->first) != field_total_.end()) {
	  field_total_[iter->first] += iter->second;
	} else {
	  field_total_[iter->first] = iter->second;
	}
      }
    }
  }

  modified_ = true;

  return added_point;
}

bool TreeMap::AddPoint(WeightedAngularCoordinate& w_ang) {
  WeightedAngularCoordinate* ang_copy =
    new WeightedAngularCoordinate(w_ang.UnitSphereX(), w_ang.UnitSphereY(),
				  w_ang.UnitSphereZ(), w_ang.Weight());
  ang_copy->CopyFields(w_ang);
  return AddPoint(ang_copy);
}

bool TreeMap::AddPoint(AngularCoordinate& ang, double object_weight) {
  WeightedAngularCoordinate* w_ang =
    new WeightedAngularCoordinate(ang.UnitSphereX(), ang.UnitSphereY(),
				  ang.UnitSphereZ(), object_weight);
  return AddPoint(w_ang);
}

bool TreeMap::Read(const std::string& input_file,
		   AngularCoordinate::Sphere sphere, bool verbose,
		   uint8_t theta_column, uint8_t phi_column,
		   int8_t weight_column) {
  bool io_success = false;

  uint32_t n_lines = 0;
  uint32_t check_lines = 128;
  if (theta_column != phi_column) {
    std::ifstream input_file_str(input_file.c_str());

    uint8_t weight_idx = static_cast<uint8_t>(weight_column);

    if (input_file_str) {
      io_success = true;
      if (verbose) std::cout << "Stomp::TreeMap::Read - " <<
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
	    double weight = 1.0;
	    if ((weight_column > -1) &&
		(line_elements.size() > weight_idx)) {
	      weight = strtod(line_elements[weight_idx].c_str(), NULL);
	    }

	    WeightedAngularCoordinate* w_ang =
	      new WeightedAngularCoordinate(theta, phi, weight, sphere);
	    if (!AddPoint(w_ang)) io_success = false;
	  }
	}
      }
      input_file_str.close();
    } else {
      std::cout << "Stomp::TreeMap::Read - " << input_file <<
	" does not exist!\n";
    }
  }

  if (verbose && io_success)
    std::cout << "Stomp::TreeMap::Read - Read " << n_lines-1 <<
      " lines from " << input_file << "; loaded " << NPoints() <<
      " into tree...\n";

  return io_success;
}

bool TreeMap::Read(const std::string& input_file,
		   FieldColumnDict& field_columns,
		   AngularCoordinate::Sphere sphere, bool verbose,
		   uint8_t theta_column, uint8_t phi_column,
		   int8_t weight_column) {
  bool io_success = false;

  uint32_t n_lines = 0;
  uint32_t check_lines = 128;
  if (theta_column != phi_column) {
    std::ifstream input_file_str(input_file.c_str());

    uint8_t weight_idx = static_cast<uint8_t>(weight_column);

    if (input_file_str) {
      io_success = true;
      if (verbose) std::cout << "Stomp::TreeMap::Read - Reading from " <<
		     input_file << "...\n";
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
	    double weight = 1.0;
	    if ((weight_column > -1) &&
		(line_elements.size() > weight_idx)) {
	      weight = strtod(line_elements[weight_idx].c_str(), NULL);
	    }

	    FieldDict fields;
	    for (FieldColumnIterator iter=field_columns.begin();
		 iter!=field_columns.end();++iter) {
	      fields[iter->first] = 0.0;
	      if (line_elements.size() > iter->second) {
		fields[iter->first] =
		  strtod(line_elements[iter->second].c_str(), NULL);
	      }
	    }

	    WeightedAngularCoordinate* w_ang =
	      new WeightedAngularCoordinate(theta, phi, weight, fields, sphere);
	    if (!AddPoint(w_ang)) io_success = false;
	  }
	}
      }
      input_file_str.close();
    } else {
      std::cout << "Stomp::TreeMap::Read - " << input_file <<
	" does not exist!\n";
    }
  }

  if (verbose && io_success)
    std::cout << "Stomp::TreeMap::Read - Read " << n_lines-1 <<
      " lines from " << input_file << "; loaded " << NPoints() <<
      " into tree...\n";

  return io_success;
}

void TreeMap::Coverage(PixelVector& superpix, uint32_t resolution,
		       bool calculate_fraction) {
  if (!superpix.empty()) superpix.clear();

  if (resolution > resolution_) {
    std::cout << "Stomp::TreeMap::Coverage - " <<
      "WARNING: Requested resolution is higher than " <<
      "the map resolution!\nReseting to map resolution...\n";
    resolution = resolution_;
  }

  // We need to make a vector of pixels that cover the current TreeMap area.
  // If the requested resolution is the same as our base node resolution, then
  // this is simple.
  if (resolution_ == resolution) {
    superpix.reserve(tree_map_.size());
    for (TreeDictIterator iter=tree_map_.begin();iter!=tree_map_.end();++iter) {
      Pixel pix(iter->second->PixelX(), iter->second->PixelY(),
		iter->second->Resolution(), iter->second->Coverage());
      superpix.push_back(pix);
    }
  } else {
    // If that's not the case, then we need to do some work.  First, the case
    // where the requested resolution is coarser than our base nodes.
    if (resolution < resolution_) {
      // We need to find the unique superpixels covered by the map and output
      // those.  We can use a temporary TreeDict to do this quickly.
      TreeDict tmp_map;
      for (TreeDictIterator iter=tree_map_.begin();
	   iter!=tree_map_.end();++iter) {
	TreePixel* pix =
	  new TreePixel(iter->second->PixelX(), iter->second->PixelY(),
			iter->second->Resolution(), 0);
	pix->SetToSuperPix(resolution);
	if (tmp_map.find(pix->Pixnum()) == tmp_map.end()) {
	  tmp_map[pix->Pixnum()] = pix;
	}
      }

      superpix.reserve(tmp_map.size());
      for (TreeDictIterator iter=tmp_map.begin();
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
      for (TreeDictIterator iter=tree_map_.begin();
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

bool TreeMap::Covering(Map& stomp_map, uint32_t maximum_pixels) {
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

double TreeMap::FindUnmaskedFraction(Pixel& pix) {
  double unmasked_fraction = 0.0;

  if (pix.Resolution() >= resolution_) {
    // If our input pixel is the size of our base-node or smaller, then we
    // can use each node's Coverage method to find the unmasked fraction,
    // provideded that a matching node can be found.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    TreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
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
      TreeDictIterator iter = tree_map_.find(sub_iter->Pixnum());
      if (iter != tree_map_.end()) {
	unmasked_fraction += pixel_fraction*iter->second->Coverage(pix);
      }
    }
  }

  return unmasked_fraction;
}

int8_t TreeMap::FindUnmaskedStatus(Pixel& pix) {
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

    TreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
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
      TreeDictIterator iter = tree_map_.find(sub_iter->Pixnum());
      if (iter != tree_map_.end()) unmasked_status = -1;

      ++sub_iter;
    }
  }

  return unmasked_status;
}

void TreeMap::NodeMap(Map& stomp_map) {
  if (!stomp_map.Empty()) stomp_map.Clear();

  PixelVector pix;
  pix.reserve(tree_map_.size());
  for (TreeDictIterator iter=tree_map_.begin();iter!=tree_map_.end();++iter) {
    Pixel tmp_pix(iter->second->PixelX(), iter->second->PixelY(),
		   iter->second->Resolution(), 1.0);
    pix.push_back(tmp_pix);
  }

  stomp_map.Initialize(pix);
}

uint32_t TreeMap::Resolution() {
  return resolution_;
}

uint16_t TreeMap::PixelCapacity() {
  return maximum_points_;
}

void TreeMap::SetResolution(uint32_t resolution) {
  Clear();
  resolution_ = resolution;
}

void TreeMap::SetPixelCapacity(int pixel_capacity) {
  Clear();
  maximum_points_ = pixel_capacity;
}

uint32_t TreeMap::NPoints(uint32_t k) {
  return (k == MaxPixnum ? point_count_ :
	  (tree_map_.find(k) != tree_map_.end() ?
	   tree_map_[k]->NPoints() : 0));
}

uint32_t TreeMap::NPoints(Pixel& pix) {
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
      TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end())
	total_points += tree_map_[pix_iter->Pixnum()]->NPoints();
    }
  } else {
    // If the input pixel is the same size as our nodes or smaller, then we
    // look for the appropriate node and return the results from that node.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    TreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end())
      total_points = tree_map_[tmp_pix.Pixnum()]->NPoints(pix);
  }

  return total_points;
}

void TreeMap::Points(WAngularVector& w_ang) {
  if (!w_ang.empty()) w_ang.clear();

  // Fairly simple, just iterate over all of the base-nodes and return the
  // aggregate vector.  We can speed things up a bit since we know the total
  // number of points in the map.
  w_ang.reserve(point_count_);

  for (TreeDictIterator iter=tree_map_.begin();
       iter!=tree_map_.end();++iter) {
    WAngularVector tmp_ang;
    iter->second->Points(tmp_ang);
    for (WAngularIterator ang_iter=tmp_ang.begin();
	 ang_iter!=tmp_ang.end();++ang_iter) w_ang.push_back(*ang_iter);
  }
}

void TreeMap::Points(WAngularVector& w_ang, Pixel& pix) {
  if (!w_ang.empty()) w_ang.clear();

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
      TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end()) {
	WAngularVector tmp_ang;

	iter->second->Points(tmp_ang);
	for (WAngularIterator ang_iter=tmp_ang.begin();
	     ang_iter!=tmp_ang.end();++ang_iter) w_ang.push_back(*ang_iter);
      }
    }
  } else {
    // If the input pixel is the same size as our nodes or smaller, then we
    // look for the appropriate node and return the results from that node.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    TreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end()) {
      WAngularVector tmp_ang;

      iter->second->Points(tmp_ang, pix);
      for (WAngularIterator ang_iter=tmp_ang.begin();
	   ang_iter!=tmp_ang.end();++ang_iter) w_ang.push_back(*ang_iter);
    }
  }
}

double TreeMap::Weight(uint32_t k) {
  return (k == Stomp::MaxPixnum ? weight_ :
	  (tree_map_.find(k) != tree_map_.end() ?
	   tree_map_[k]->Weight() : 0.0));
}

double TreeMap::Weight(Pixel& pix) {
  double total_weight = 0;

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
      TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end())
	total_weight += tree_map_[pix_iter->Pixnum()]->Weight();
    }
  } else {
    // If the input pixel is the same size as our nodes or smaller, then we
    // look for the appropriate node and return the results for that node.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    TreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end())
      total_weight = tree_map_[tmp_pix.Pixnum()]->PixelWeight(pix);
  }

  return total_weight;
}

double TreeMap::FieldTotal(const std::string& field_name, uint32_t k) {
  return (k == MaxPixnum ?
	  (field_total_.find(field_name) != field_total_.end() ?
	   field_total_[field_name] : 0.0) :
	  (tree_map_.find(k) != tree_map_.end() ?
	   tree_map_[k]->FieldTotal(field_name) : 0.0));
}

double TreeMap::FieldTotal(const std::string& field_name, Pixel& pix) {
  double field_total = 0;

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
      TreeDictIterator iter = tree_map_.find(pix_iter->Pixnum());
      if (iter != tree_map_.end())
	field_total += tree_map_[pix_iter->Pixnum()]->FieldTotal(field_name);
    }
  } else {
    // If the input pixel is the same size as our nodes or smaller, then we
    // look for the appropriate node and return the results for that node.
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    TreeDictIterator iter = tree_map_.find(tmp_pix.Pixnum());
    if (iter != tree_map_.end())
      field_total = tree_map_[tmp_pix.Pixnum()]->FieldTotal(field_name, pix);
  }

  return field_total;
}

uint16_t TreeMap::NField() {
  return field_total_.size();
}

bool TreeMap::HasFields() {
  return (field_total_.size() > 0 ? true : false);
}

void TreeMap::FieldNames(std::vector<std::string>& field_names) {
  field_names.clear();
  for (FieldIterator iter=field_total_.begin();
       iter!=field_total_.end();++iter) field_names.push_back(iter->first);
}

uint16_t TreeMap::BaseNodes() {
  return tree_map_.size();
}

uint16_t TreeMap::Nodes() {
  uint16_t total_nodes = 0;
  for (TreeDictIterator iter=tree_map_.begin();
       iter!=tree_map_.end();++iter) total_nodes += iter->second->Nodes();

  return total_nodes;
}

uint32_t TreeMap::Size() {
  return point_count_;
}

double TreeMap::Area() {
  if (modified_) CalculateArea();

  return area_;
}

void TreeMap::CalculateArea() {
  area_ = 0.0;
  for (TreeDictIterator iter=tree_map_.begin();
       iter!=tree_map_.end();++iter) {
    area_ += iter->second->Coverage()*iter->second->Area();
  }

  modified_ = false;
}

uint32_t TreeMap::MinResolution() {
  return resolution_;
}

uint32_t TreeMap::MaxResolution() {
  return resolution_;
}

uint8_t TreeMap::MinLevel() {
  return Pixel::ResolutionToLevel(resolution_);
}

uint8_t TreeMap::MaxLevel() {
  return Pixel::ResolutionToLevel(resolution_);
}

bool TreeMap::Empty() {
  return (tree_map_.empty() ? true : false);
}

void TreeMap::Clear() {
  if (!tree_map_.empty()) {
    for (TreeDictIterator iter=tree_map_.begin();
	 iter!=tree_map_.end();++iter) {
      iter->second->Clear();
      delete iter->second;
    }
    tree_map_.clear();
    field_total_.clear();
    weight_ = 0.0;
    area_ = 0.0;
    point_count_ = 0;
    modified_ = false;
  }
  ClearRegions();
}

} // end namespace Stomp
