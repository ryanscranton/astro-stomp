// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains a variant on the Pixel class.  The goal here is to
// use the hierarchical nature of the Pixel class to form the basis for a
// spatial quad tree structure.  Hence, a given TreePixel object will have a
// number of points associated with it, where they have been stored in such a
// way that pair finding and K nearest neighbor searches will run in ln(N) time.

#include "stomp_core.h"
#include "stomp_tree_pixel.h"
#include "stomp_angular_bin.h"
#include "stomp_radial_bin.h"
#include "stomp_angular_correlation.h"

namespace Stomp {

TreePixel::TreePixel() {
  SetWeight(0.0);
  maximum_points_ = 0;
  point_count_ = 0;
  InitializeCorners();
}

TreePixel::TreePixel(const uint32_t input_resolution,
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

TreePixel:: TreePixel(const uint32_t input_x,
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

TreePixel::TreePixel(AngularCoordinate& ang,
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

TreePixel::~TreePixel() {
  ang_.clear();
  subpix_.clear();
  maximum_points_ = 0;
  point_count_ = 0;
  initialized_subpixels_ = false;
}

bool TreePixel::_InitializeSubPixels() {
  initialized_subpixels_ = false;
  // If we're already at the maximum resolution, then we shouldn't be trying
  // to generate sub-pixels.
  if (Resolution() < Stomp::MaxPixelResolution) {
    PixelVector tmp_pix;
    SubPix(Resolution()*2, tmp_pix);
    subpix_.reserve(4);

    // Provided we passed that test, we create a vector of sub-pixels.
    for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter) {
      TreePixel* tree_pix = new TreePixel(iter->PixelX(), iter->PixelY(),
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

uint32_t TreePixel::DirectPairCount(AngularCoordinate& ang,
				    AngularBin& theta,
				    int16_t region) {
  uint32_t pair_count = 0;
  if (theta.ThetaMax() < 90.0) {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinCosBounds((*iter)->DotProduct(ang))) pair_count++;
    }
  } else {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinBounds((*iter)->AngularDistance(ang))) pair_count++;
    }
  }

  theta.AddToCounter(pair_count, region);
  return pair_count;
}

uint32_t TreePixel::FindPairs(AngularCoordinate& ang, AngularBin& theta,
			      int16_t region) {
  uint32_t pair_count = 0;

  // If we have AngularCoordinates in this pixel, then this is just a
  // matter of iterating through them and finding how many satisfy the
  // angular bounds.
  if (!ang_.empty()) {
    pair_count = DirectPairCount(ang, theta, region);
  } else {
    // If the current pixel doesn't contain any points, then we need to see
    // if either the current pixel is either fully or partially contained in
    // the annulus.  For the former case, we can just send back the total
    // number of points in this pixel.  In the latter case, we pass things
    // along to the sub-pixels.  If neither of those things are true, then
    // we're done and we send back zero.
    int8_t intersects_annulus = IntersectsAnnulus(ang, theta);

    if (intersects_annulus == 1) {
      // Fully contained in the annulus.
      pair_count = point_count_;
      theta.AddToCounter(point_count_, region);
    } else {
      if (intersects_annulus == -1) {
      // Partial intersection with the annulus.
	for (TreePtrIterator iter=subpix_.begin();
	     iter!=subpix_.end();++iter) {
	  pair_count += (*iter)->FindPairs(ang, theta, region);
	}
      } else {
	// Completely outside the annulus.
	pair_count = 0;
      }
    }
  }
  return pair_count;
}

uint32_t TreePixel::FindPairs(AngularCoordinate& ang,
			      double theta_min, double theta_max) {
  AngularBin theta(theta_min, theta_max);
  return FindPairs(ang, theta);
}

uint32_t TreePixel::FindPairs(AngularCoordinate& ang, double theta_max) {
  AngularBin theta(0.0, theta_max);
  return FindPairs(ang, theta);
}

double TreePixel::DirectWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
				      int16_t region) {
  double total_weight = 0.0;
  uint32_t n_pairs = 0;

  if (theta.ThetaMax() < 90.0) {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinCosBounds((*iter)->DotProduct(ang))) {
	total_weight += (*iter)->Weight();
	n_pairs++;
      }
    }
  } else {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinBounds((*iter)->AngularDistance(ang))) {
	total_weight += (*iter)->Weight();
	n_pairs++;
      }
    }
  }

  theta.AddToWeight(total_weight, region);
  theta.AddToCounter(n_pairs, region);

  return total_weight;
}

double TreePixel::FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
				    int16_t region) {
  double total_weight = 0.0;
  // If we have AngularCoordinates in this pixel, then this is just a
  // matter of iterating through them and finding how many satisfy the
  // angular bounds.
  if (!ang_.empty()) {
    total_weight = DirectWeightedPairs(ang, theta, region);
  } else {
    // If the current pixel doesn't contain any points, then we need to see
    // if either the current pixel is either fully or partially contained in
    // the annulus.  For the former case, we can just send back the total
    // number of points in this pixel.  In the latter case, we pass things
    // along to the sub-pixels.  If neither of those things are true, then
    // we're done and we send back zero.
    int8_t intersects_annulus = IntersectsAnnulus(ang, theta);
    if (intersects_annulus == 1) {
      // Fully contained in the annulus.
      total_weight = Weight();
      theta.AddToWeight(Weight(), region);
      theta.AddToCounter(point_count_, region);
    } else {
      if (intersects_annulus == -1) {
      // Partial intersection with the annulus.
	for (TreePtrIterator iter=subpix_.begin();
	     iter!=subpix_.end();++iter) {
	  total_weight += (*iter)->FindWeightedPairs(ang, theta, region);
	}
      } else {
	// Completely outside the annulus.
	total_weight = 0.0;
      }
    }
  }
  return total_weight;
}

double TreePixel::FindWeightedPairs(AngularCoordinate& ang,
				    double theta_min, double theta_max) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(ang, theta);
}

double TreePixel::FindWeightedPairs(AngularCoordinate& ang, double theta_max) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(ang, theta);
}

double TreePixel::DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
				      AngularBin& theta, int16_t region) {
  double total_weight = 0.0;
  uint32_t n_pairs = 0;

  if (theta.ThetaMax() < 90.0) {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinCosBounds((*iter)->DotProduct(w_ang))) {
	total_weight += (*iter)->Weight();
	n_pairs++;
      }
    }
  } else {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinBounds((*iter)->AngularDistance(w_ang))) {
	total_weight += (*iter)->Weight();
	n_pairs++;
      }
    }
  }

  total_weight *= w_ang.Weight();

  theta.AddToWeight(total_weight, region);
  theta.AddToCounter(n_pairs, region);

  return total_weight;
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    AngularBin& theta, int16_t region) {
  double total_weight = 0.0;
  // If we have AngularCoordinates in this pixel, then this is just a
  // matter of iterating through them and finding how many satisfy the
  // angular bounds.
  if (!ang_.empty()) {
    total_weight = DirectWeightedPairs(w_ang, theta, region);
  } else {
    // If the current pixel doesn't contain any points, then we need to see
    // if either the current pixel is either fully or partially contained in
    // the annulus.  For the former case, we can just send back the total
    // number of points in this pixel.  In the latter case, we pass things
    // along to the sub-pixels.  If neither of those things are true, then
    // we're done and we send back zero.
    int8_t intersects_annulus = IntersectsAnnulus(w_ang, theta);
    if (intersects_annulus == 1) {
      // Fully contained in the annulus.
      total_weight = Weight()*w_ang.Weight();
      theta.AddToWeight(Weight(), region);
      theta.AddToCounter(point_count_, region);
    } else {
      if (intersects_annulus == -1) {
      // Partial intersection with the annulus.
	for (TreePtrIterator iter=subpix_.begin();
	     iter!=subpix_.end();++iter) {
	  total_weight += (*iter)->FindWeightedPairs(w_ang, theta, region);
	}
      } else {
	// Completely outside the annulus.
	total_weight = 0.0;
      }
    }
  }
  return total_weight;
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    double theta_min, double theta_max) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(w_ang, theta);
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    double theta_max) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(w_ang, theta);
}

void TreePixel::FindPairs(AngularVector& ang, AngularBin& theta,
			  int16_t region) {
  uint32_t n_pairs = 0;
  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    n_pairs = FindPairs(*ang_iter, theta, region);
  }
}

void TreePixel::FindPairs(AngularVector& ang, AngularCorrelation& wtheta,
			  int16_t region) {
  uint32_t n_pairs = 0;
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (AngularIterator ang_iter=ang.begin();
	 ang_iter!=ang.end();++ang_iter) {
      n_pairs = FindPairs(*ang_iter, *theta_iter, region);
    }
  }
}

void TreePixel::FindWeightedPairs(AngularVector& ang, AngularBin& theta,
				  int16_t region) {
  double total_weight = 0.0;
  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, theta, region);
  }
}

void TreePixel::FindWeightedPairs(AngularVector& ang,
				  AngularCorrelation& wtheta,
				  int16_t region) {
  double total_weight = 0.0;
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (AngularIterator ang_iter=ang.begin();
	 ang_iter!=ang.end();++ang_iter) {
      total_weight = FindWeightedPairs(*ang_iter, *theta_iter, region);
    }
  }
}

void TreePixel::FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta,
				  int16_t region) {
  double total_weight = 0.0;
  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, theta, region);
  }
}

void TreePixel::FindWeightedPairs(WAngularVector& w_ang,
				  AngularCorrelation& wtheta, int16_t region) {
  double total_weight = 0.0;
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (WAngularIterator ang_iter=w_ang.begin();
	 ang_iter!=w_ang.end();++ang_iter) {
      total_weight = FindWeightedPairs(*ang_iter, *theta_iter, region);
    }
  }
}

double TreePixel::DirectWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
				      const std::string& field_name,
				      int16_t region) {
  double total_weight = 0.0;
  uint32_t n_pairs = 0;

  if (theta.ThetaMax() < 90.0) {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinCosBounds((*iter)->DotProduct(ang))) {
	total_weight += (*iter)->Field(field_name);
	n_pairs++;
      }
    }
  } else {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinBounds((*iter)->AngularDistance(ang))) {
	total_weight += (*iter)->Field(field_name);
	n_pairs++;
      }
    }
  }

  theta.AddToWeight(total_weight, region);
  theta.AddToCounter(n_pairs, region);

  return total_weight;
}

double TreePixel::FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
				    const std::string& field_name,
				    int16_t region) {
  double total_weight = 0.0;
  // If we have AngularCoordinates in this pixel, then this is just a
  // matter of iterating through them and finding how many satisfy the
  // angular bounds.
  if (!ang_.empty()) {
    total_weight = DirectWeightedPairs(ang, theta, field_name, region);
  } else {
    // If the current pixel doesn't contain any points, then we need to see
    // if either the current pixel is either fully or partially contained in
    // the annulus.  For the former case, we can just send back the total
    // number of points in this pixel.  In the latter case, we pass things
    // along to the sub-pixels.  If neither of those things are true, then
    // we're done and we send back zero.
    int8_t intersects_annulus = IntersectsAnnulus(ang, theta);
    if (intersects_annulus == 1) {
      // Fully contained in the annulus.
      total_weight = FieldTotal(field_name);
      theta.AddToWeight(total_weight, region);
      theta.AddToCounter(point_count_, region);
    } else {
      if (intersects_annulus == -1) {
      // Partial intersection with the annulus.
	for (TreePtrIterator iter=subpix_.begin();
	     iter!=subpix_.end();++iter) {
	  total_weight +=
	    (*iter)->FindWeightedPairs(ang, theta, field_name, region);
	}
      } else {
	// Completely outside the annulus.
	total_weight = 0.0;
      }
    }
  }
  return total_weight;
}

double TreePixel::FindWeightedPairs(AngularCoordinate& ang,
				    double theta_min, double theta_max,
				    const std::string& field_name) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(ang, theta, field_name);
}

double TreePixel::FindWeightedPairs(AngularCoordinate& ang, double theta_max,
				    const std::string& field_name) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(ang, theta, field_name);
}

void TreePixel::FindWeightedPairs(AngularVector& ang, AngularBin& theta,
				  const std::string& field_name,
				  int16_t region) {
  double total_weight = 0.0;
  for (AngularIterator ang_iter=ang.begin();ang_iter!=ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, theta, field_name, region);
  }
}

void TreePixel::FindWeightedPairs(AngularVector& ang,
				  AngularCorrelation& wtheta,
				  const std::string& field_name, int16_t region) {
  double total_weight = 0.0;
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (AngularIterator ang_iter=ang.begin();
	 ang_iter!=ang.end();++ang_iter) {
      total_weight =
	FindWeightedPairs(*ang_iter, *theta_iter, field_name, region);
    }
  }
}

double TreePixel::DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
				      AngularBin& theta,
				      const std::string& field_name,
				      int16_t region) {
  double total_weight = 0.0;
  uint32_t n_pairs = 0;

  if (theta.ThetaMax() < 90.0) {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinCosBounds((*iter)->DotProduct(w_ang))) {
	total_weight += (*iter)->Field(field_name);
	n_pairs++;
      }
    }
  } else {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinBounds((*iter)->AngularDistance(w_ang))) {
	total_weight += (*iter)->Field(field_name);
	n_pairs++;
      }
    }
  }

  total_weight *= w_ang.Weight();

  theta.AddToWeight(total_weight, region);
  theta.AddToCounter(n_pairs, region);

  return total_weight;
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    AngularBin& theta,
				    const std::string& field_name,
				    int16_t region) {
  double total_weight = 0.0;
  // If we have AngularCoordinates in this pixel, then this is just a
  // matter of iterating through them and finding how many satisfy the
  // angular bounds.
  if (!ang_.empty()) {
    total_weight = DirectWeightedPairs(w_ang, theta, field_name, region);
  } else {
    // If the current pixel doesn't contain any points, then we need to see
    // if either the current pixel is either fully or partially contained in
    // the annulus.  For the former case, we can just send back the total
    // number of points in this pixel.  In the latter case, we pass things
    // along to the sub-pixels.  If neither of those things are true, then
    // we're done and we send back zero.
    int8_t intersects_annulus = IntersectsAnnulus(w_ang, theta);
    if (intersects_annulus == 1) {
      // Fully contained in the annulus.
      total_weight = FieldTotal(field_name)*w_ang.Weight();
      theta.AddToWeight(total_weight, region);
      theta.AddToCounter(point_count_, region);
    } else {
      if (intersects_annulus == -1) {
      // Partial intersection with the annulus.
	for (TreePtrIterator iter=subpix_.begin();
	     iter!=subpix_.end();++iter) {
	  total_weight +=
	    (*iter)->FindWeightedPairs(w_ang, theta, field_name, region);
	}
      } else {
	// Completely outside the annulus.
	total_weight = 0.0;
      }
    }
  }
  return total_weight;
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    double theta_min, double theta_max,
				    const std::string& field_name) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(w_ang, theta, field_name);
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    double theta_max,
				    const std::string& field_name) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(w_ang, theta, field_name);
}

void TreePixel::FindWeightedPairs(WAngularVector& w_ang,
				  AngularBin& theta,
				  const std::string& field_name,
				  int16_t region) {
  double total_weight = 0.0;
  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, theta, field_name, region);
  }
}

void TreePixel::FindWeightedPairs(WAngularVector& w_ang,
				  AngularCorrelation& wtheta,
				  const std::string& field_name,
				  int16_t region) {
  double total_weight = 0.0;
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (WAngularIterator ang_iter=w_ang.begin();
	 ang_iter!=w_ang.end();++ang_iter) {
      total_weight =
	FindWeightedPairs(*ang_iter, *theta_iter, field_name, region);
    }
  }
}

double TreePixel::DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
				      const std::string& ang_field_name,
				      AngularBin& theta,
				      const std::string& field_name,
				      int16_t region) {
  double total_weight = 0.0;
  uint32_t n_pairs = 0;

  if (theta.ThetaMax() < 90.0) {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinCosBounds((*iter)->DotProduct(w_ang))) {
	total_weight += (*iter)->Field(field_name);
	n_pairs++;
      }
    }
  } else {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      if (theta.WithinBounds((*iter)->AngularDistance(w_ang))) {
	total_weight += (*iter)->Field(field_name);
	n_pairs++;
      }
    }
  }

  total_weight *= w_ang.Field(ang_field_name);

  theta.AddToWeight(total_weight, region);
  theta.AddToCounter(n_pairs, region);

  return total_weight;
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    const std::string& ang_field_name,
				    AngularBin& theta,
				    const std::string& field_name,
				    int16_t region) {
  double total_weight = 0.0;
  // If we have AngularCoordinates in this pixel, then this is just a
  // matter of iterating through them and finding how many satisfy the
  // angular bounds.
  if (!ang_.empty()) {
    total_weight =
      DirectWeightedPairs(w_ang, ang_field_name, theta, field_name, region);
  } else {
    // If the current pixel doesn't contain any points, then we need to see
    // if either the current pixel is either fully or partially contained in
    // the annulus.  For the former case, we can just send back the total
    // number of points in this pixel.  In the latter case, we pass things
    // along to the sub-pixels.  If neither of those things are true, then
    // we're done and we send back zero.
    int8_t intersects_annulus = IntersectsAnnulus(w_ang, theta);
    if (intersects_annulus == 1) {
      // Fully contained in the annulus.
      total_weight = FieldTotal(field_name)*w_ang.Field(ang_field_name);
      theta.AddToWeight(total_weight, region);
      theta.AddToCounter(point_count_, region);
    } else {
      if (intersects_annulus == -1) {
      // Partial intersection with the annulus.
	for (TreePtrIterator iter=subpix_.begin();
	     iter!=subpix_.end();++iter) {
	  total_weight +=
	    (*iter)->FindWeightedPairs(w_ang, ang_field_name, theta,
				       field_name, region);
	}
      } else {
	// Completely outside the annulus.
	total_weight = 0.0;
      }
    }
  }
  return total_weight;
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    const std::string& ang_field_name,
				    double theta_min, double theta_max,
				    const std::string& field_name) {
  AngularBin theta(theta_min, theta_max);
  return FindWeightedPairs(w_ang, ang_field_name, theta, field_name);
}

double TreePixel::FindWeightedPairs(WeightedAngularCoordinate& w_ang,
				    const std::string& ang_field_name,
				    double theta_max,
				    const std::string& field_name) {
  AngularBin theta(0.0, theta_max);
  return FindWeightedPairs(w_ang, ang_field_name, theta, field_name);
}

void TreePixel::FindWeightedPairs(WAngularVector& w_ang,
				  const std::string& ang_field_name,
				  AngularBin& theta,
				  const std::string& field_name,
				  int16_t region) {
  double total_weight = 0.0;
  for (WAngularIterator ang_iter=w_ang.begin();
       ang_iter!=w_ang.end();++ang_iter) {
    total_weight = FindWeightedPairs(*ang_iter, ang_field_name, theta,
				     field_name, region);
  }
}

void TreePixel::FindWeightedPairs(WAngularVector& w_ang,
				  const std::string& ang_field_name,
				  AngularCorrelation& wtheta,
				  const std::string& field_name,
				  int16_t region) {
  double total_weight = 0.0;
  for (ThetaIterator theta_iter=wtheta.Begin(0);
       theta_iter!=wtheta.End(0);++theta_iter) {
    for (WAngularIterator ang_iter=w_ang.begin();
	 ang_iter!=w_ang.end();++ang_iter) {
      total_weight =
	FindWeightedPairs(*ang_iter, ang_field_name, *theta_iter,
			  field_name, region);
    }
  }
}

uint16_t TreePixel::FindKNearestNeighbors(AngularCoordinate& ang,
					  uint8_t n_neighbors,
					  WAngularVector& neighbor_ang) {
  TreeNeighbor neighbors(ang, n_neighbors);

  _NeighborRecursion(ang, neighbors);

  neighbors.NearestNeighbors(neighbor_ang, false);

  return neighbors.NodesVisited();
}

uint16_t TreePixel::FindNearestNeighbor(AngularCoordinate& ang,
					WeightedAngularCoordinate& nbr_ang) {
  WAngularVector angVec;

  uint16_t nodes_visited = FindKNearestNeighbors(ang, 1, angVec);

  nbr_ang = angVec[0];

  return nodes_visited;
}

double TreePixel::KNearestNeighborDistance(AngularCoordinate& ang,
					   uint8_t n_neighbors,
					   uint16_t& nodes_visited) {

  TreeNeighbor neighbors(ang, n_neighbors);

  _NeighborRecursion(ang, neighbors);

  nodes_visited = neighbors.NodesVisited();

  return neighbors.MaxAngularDistance();
}

double TreePixel::NearestNeighborDistance(AngularCoordinate& ang,
					  uint16_t& nodes_visited) {
  return KNearestNeighborDistance(ang, 1, nodes_visited);
}

bool TreePixel::ClosestMatch(AngularCoordinate& ang, double max_distance,
			     WeightedAngularCoordinate& match_ang) {
  TreeNeighbor neighbors(ang, 1, max_distance);

  _NeighborRecursion(ang, neighbors);

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

void TreePixel::_NeighborRecursion(AngularCoordinate& ang,
				   TreeNeighbor& neighbors) {

  neighbors.AddNode();

  if (!ang_.empty()) {
    // We have no sub-nodes in this tree, so we'll just iterate over the
    // points here and take the nearest N neighbors.
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter)
      neighbors.TestPoint(*iter);
  } else {
    // This node is the root node for our tree, so we first find the sub-node
    // that contains the point and start recursing there.
    //
    // While we iterate through the nodes, we'll also calculate the edge
    // distances for those nodes that don't contain the point and store them
    // in a priority queue.  This will let us do a follow-up check on nodes in
    // the most productive order.
    PixelQueue pix_queue;
    for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
      if ((*iter)->Contains(ang)) {
	(*iter)->_NeighborRecursion(ang, neighbors);
      } else {
	double min_edge_distance, max_edge_distance;
	(*iter)->EdgeDistances(ang, min_edge_distance, max_edge_distance);
	DistancePixelPair dist_pair(min_edge_distance, (*iter));
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
      TreePixel* pix_iter = pix_queue.top().second;
      if (pix_distance < neighbors.MaxDistance()) {
	pix_iter->_NeighborRecursion(ang, neighbors);
      }
      pix_queue.pop();
    }
  }
}

void TreePixel::InitializeCorners() {
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

bool TreePixel::AddPoint(WeightedAngularCoordinate* ang) {
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
	  std::cout << "Stomp::TreePixel::AddPoint - " <<
	    "Failed to initialize sub-pixels.  Exiting.\n";
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

  if (added_to_pixel) {
    AddToWeight(ang->Weight());
    if (ang->HasFields()) {
      for (FieldIterator iter=ang->FieldBegin();iter!=ang->FieldEnd();++iter) {
	if (field_total_.find(iter->first) != field_total_.end()) {
	  field_total_[iter->first] += iter->second;
	} else {
	  field_total_[iter->first] = iter->second;
	}
      }
    }
    point_count_++;
  }

  return added_to_pixel;
}

bool TreePixel::AddPoint(WeightedAngularCoordinate& w_ang) {
  WeightedAngularCoordinate* ang_copy =
    new WeightedAngularCoordinate(w_ang.UnitSphereX(), w_ang.UnitSphereY(),
				  w_ang.UnitSphereZ(), w_ang.Weight());
  ang_copy->CopyFields(w_ang);
  return AddPoint(ang_copy);
}

bool TreePixel::AddPoint(AngularCoordinate& ang, double object_weight) {
  WeightedAngularCoordinate* w_ang =
    new WeightedAngularCoordinate(ang.UnitSphereX(), ang.UnitSphereY(),
				  ang.UnitSphereZ(), object_weight);
  return AddPoint(w_ang);
}

uint32_t TreePixel::NPoints() {
  return point_count_;
}

uint32_t TreePixel::NPoints(Pixel& pix) {
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
	for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	  total_points += (*iter)->NPoints(pix);
	}
      } else {
	for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	  if (pix.Contains(*(*iter))) total_points++;
	}
      }
    }
  }

  return total_points;
}

double TreePixel::PixelWeight(Pixel& pix) {
  double total_weight = 0.0;

  // First check to see if the input pixel contains the current pixel.
  if (pix.Contains(Resolution(), PixelX(), PixelY())) {
    // If so, then it also contains all of the points in the current pixel.
    total_weight = Weight();
  } else {
    // If not, then either the input pixel doesn't overlap the current one or
    // it's a sub-pixel of the current one.
    if (Contains(pix)) {
      // If we contain the input pixel, then we either iterate over the
      // sub-nodes to this pixel or iterate over the points contained in this
      // pixel.
      if (initialized_subpixels_) {
	for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	  total_weight += (*iter)->PixelWeight(pix);
	}
      } else {
	for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	  if (pix.Contains(*(*iter))) total_weight += (*iter)->Weight();
	}
      }
    }
  }

  return total_weight;
}

double TreePixel::Coverage() {
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
    for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
      total_coverage += 0.25*(*iter)->Coverage();
    }
  }

  return total_coverage;
}

double TreePixel::Coverage(Pixel& pix) {
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
	for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	  if ((*iter)->Contains(pix))
	    total_coverage = (*iter)->Coverage(pix);
	}
      }
    }
  }

  return total_coverage;
}

void TreePixel::Points(WAngularVector& w_ang) {
  if (!w_ang.empty()) w_ang.clear();
  w_ang.reserve(point_count_);

  // If we haven't initialized any sub-nodes, then this is just a matter of
  // creating a copy of all of the points in the current pixel.
  if (!initialized_subpixels_) {
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
      WeightedAngularCoordinate tmp_ang = *(*iter);

      w_ang.push_back(tmp_ang);
    }
  } else {
    // If not, then we need to iterate through our sub-nodes and return an
    // aggregate list.
    for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
      WAngularVector tmp_ang;
      (*iter)->Points(tmp_ang);
      for (WAngularIterator ang_iter=tmp_ang.begin();
	   ang_iter!=tmp_ang.end();++ang_iter) w_ang.push_back(*ang_iter);
    }
  }
}

void TreePixel::Points(WAngularVector& w_ang, Pixel& pix) {
  if (!w_ang.empty()) w_ang.clear();
  w_ang.reserve(point_count_);

  // First, we need to check to verify that the input pixel is either contained
  // in the current pixel or contains it.
  if (Contains(pix) || pix.Contains(Resolution(), PixelX(), PixelY())) {
    // If we haven't initialized any sub-nodes, then this is just a matter of
    // creating a copy of all of the points in the current pixel that are
    // contained in the input pixel.
    if (!initialized_subpixels_) {
      for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	if (pix.Contains(*(*iter))) {
	  WeightedAngularCoordinate tmp_ang((*iter)->Lambda(),
					    (*iter)->Eta(),
					    (*iter)->Weight(),
					    AngularCoordinate::Survey);
	  if ((*iter)->NFields() > 0) {
	    std::vector<std::string> field_names;
	    (*iter)->FieldNames(field_names);
	    for (uint16_t i=0;i<(*iter)->NFields();i++)
	      tmp_ang.SetField(field_names[i], (*iter)->Field(field_names[i]));
	  }
	  w_ang.push_back(tmp_ang);
	}
      }
    } else {
      // If not, then we need to iterate through our sub-nodes and return an
      // aggregate list.
      for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	WAngularVector tmp_ang;
	(*iter)->Points(tmp_ang, pix);
	for (WAngularIterator ang_iter=tmp_ang.begin();
	     ang_iter!=tmp_ang.end();++ang_iter) w_ang.push_back(*ang_iter);
      }
    }
  }
}

uint16_t TreePixel::Nodes() {
  uint16_t n_nodes = 1;

  for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter)
    (*iter)->_AddSubNodes(n_nodes);

  return n_nodes;
}

void TreePixel::_AddSubNodes(uint16_t& n_nodes) {
  n_nodes++;

  for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter)
    (*iter)->_AddSubNodes(n_nodes);
}

void TreePixel::AddToWeight(double weight) {
  SetWeight(Weight() + weight);
}

double TreePixel::FieldTotal(const std::string& field_name) {
  return (field_total_.find(field_name) != field_total_.end() ?
	  field_total_[field_name] : 0.0);
}

double TreePixel::FieldTotal(const std::string& field_name, Pixel& pix) {
  double total_field = 0.0;

  // First check to see if the input pixel contains the current pixel.
  if (pix.Contains(Resolution(), PixelX(), PixelY())) {
    // If so, then it also contains all of the points in the current pixel.
    total_field = FieldTotal(field_name);
  } else {
    // If not, then either the input pixel doesn't overlap the current one or
    // it's a sub-pixel of the current one.
    if (Contains(pix)) {
      // If we contain the input pixel, then we either iterate over the
      // sub-nodes to this pixel or iterate over the points contained in this
      // pixel.
      if (initialized_subpixels_) {
	for (TreePtrIterator iter=subpix_.begin();iter!=subpix_.end();++iter) {
	  total_field += (*iter)->FieldTotal(field_name, pix);
	}
      } else {
	for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter) {
	  if (pix.Contains(*(*iter)))
	    total_field += (*iter)->Field(field_name);
	}
      }
    }
  }

  return total_field;
}

void TreePixel::AddToField(const std::string& field_name, double weight) {
  if (field_total_.find(field_name) != field_total_.end()) {
    field_total_[field_name] += weight;
  } else {
    field_total_[field_name] = weight;
  }
}

uint16_t TreePixel::NField() {
  return field_total_.size();
}

bool TreePixel::HasFields() {
  return (field_total_.size() > 0 ? true : false);
}

void TreePixel::FieldNames(std::vector<std::string>& field_names) {
  field_names.clear();
  for (FieldIterator iter=field_total_.begin();
       iter!=field_total_.end();++iter) field_names.push_back(iter->first);
}

void TreePixel::SetPixelCapacity(uint16_t maximum_points) {
  maximum_points_ = maximum_points;
}

uint16_t TreePixel::PixelCapacity() {
  return maximum_points_;
}

WAngularPtrIterator TreePixel::PointsBegin() {
  return ang_.begin();
}

WAngularPtrIterator TreePixel::PointsEnd() {
  return ang_.end();
}

TreePtrIterator TreePixel::NodesBegin() {
  return subpix_.begin();
}

TreePtrIterator TreePixel::NodesEnd() {
  return subpix_.end();
}

bool TreePixel::HasPoints() {
  return (ang_.empty() ? false : true);
}

bool TreePixel::HasNodes() {
  return (subpix_.empty() ? false : true);
}

void TreePixel::Clear() {
  if (!ang_.empty())
    for (WAngularPtrIterator iter=ang_.begin();iter!=ang_.end();++iter)
      delete *iter;
  ang_.clear();
  if (!subpix_.empty())
    for (uint32_t i=0;i<subpix_.size();i++) {
      subpix_[i]->Clear();
      delete subpix_[i];
    }
  subpix_.clear();
}

double TreePixel::UnitSphereX() {
  return unit_sphere_x_;
}

double TreePixel::UnitSphereY() {
  return unit_sphere_y_;
}

double TreePixel::UnitSphereZ() {
  return unit_sphere_z_;
}

double TreePixel::UnitSphereX_UL() {
  return unit_sphere_x_ul_;
}

double TreePixel::UnitSphereY_UL() {
  return unit_sphere_y_ul_;
}

double TreePixel::UnitSphereZ_UL() {
  return unit_sphere_z_ul_;
}

double TreePixel::UnitSphereX_UR() {
  return unit_sphere_x_ur_;
}

double TreePixel::UnitSphereY_UR() {
  return unit_sphere_y_ur_;
}

double TreePixel::UnitSphereZ_UR() {
  return unit_sphere_z_ur_;
}

double TreePixel::UnitSphereX_LL() {
  return unit_sphere_x_ll_;
}

double TreePixel::UnitSphereY_LL() {
  return unit_sphere_y_ll_;
}

double TreePixel::UnitSphereZ_LL() {
  return unit_sphere_z_ll_;
}

double TreePixel::UnitSphereX_LR() {
  return unit_sphere_x_lr_;
}

double TreePixel::UnitSphereY_LR() {
  return unit_sphere_y_lr_;
}

double TreePixel::UnitSphereZ_LR() {
  return unit_sphere_z_lr_;
}

void TreePixel::WithinAnnulus(AngularBin& theta, PixelVector& pix,
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
      TreePixel tree_pix(x, y, Resolution());
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

TreeNeighbor::TreeNeighbor(AngularCoordinate& reference_ang,
			   uint8_t n_neighbor) {
  reference_ang_ = reference_ang;
  n_neighbors_ = n_neighbor;
  max_distance_ = 100.0;
  n_nodes_visited_ = 0;
}

TreeNeighbor::TreeNeighbor(AngularCoordinate& reference_ang,
			   uint8_t n_neighbor, double max_distance) {
  reference_ang_ = reference_ang;
  n_neighbors_ = n_neighbor;
  max_distance_ = sin(DegToRad*max_distance)*sin(DegToRad*max_distance);
  n_nodes_visited_ = 0;
}

TreeNeighbor::~TreeNeighbor() {
  n_neighbors_ = 0;
  max_distance_ = 100.0;
  n_nodes_visited_ = 0;
}

void TreeNeighbor::NearestNeighbors(WAngularVector& w_ang,
				    bool save_neighbors) {
  if (!w_ang.empty()) w_ang.clear();

  std::vector<DistancePointPair> backup_copy;

  while (!ang_queue_.empty()) {
    DistancePointPair dist_pair = ang_queue_.top();
    ang_queue_.pop();

    WeightedAngularCoordinate tmp_ang(dist_pair.second->UnitSphereX(),
				      dist_pair.second->UnitSphereY(),
				      dist_pair.second->UnitSphereZ(),
				      dist_pair.second->Weight());
    tmp_ang.CopyFields(dist_pair.second);

    w_ang.push_back(tmp_ang);
    backup_copy.push_back(dist_pair);
  }

  if (save_neighbors) {
    for (uint8_t i=0;i<backup_copy.size();i++) {
      ang_queue_.push(backup_copy[i]);
    }
  }
}

uint8_t TreeNeighbor::Neighbors() {
  return ang_queue_.size();
}

uint8_t TreeNeighbor::MaxNeighbors() {
  return n_neighbors_;
}

bool TreeNeighbor::TestPoint(WeightedAngularCoordinate* test_ang) {
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
    DistancePointPair dist_pair(sin2theta, test_ang);
    ang_queue_.push(dist_pair);

    // And reset our maximum distance using the new top of the heap.
    max_distance_ = ang_queue_.top().first;
  }

  return kept_point;
}

double TreeNeighbor::MaxDistance() {
  return max_distance_;
}

double TreeNeighbor::MaxAngularDistance() {
  // Numerical precision can sometimes make the max_distance_ negative.
  return RadToDeg*asin(sqrt(fabs(max_distance_)));
}

uint16_t TreeNeighbor::NodesVisited() {
  return n_nodes_visited_;
}

void TreeNeighbor::AddNode() {
  n_nodes_visited_++;
}

} // end namespace Stomp
