// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the Map class.  Maps are intended to describe an
// arbitrary area on the sky as precisely as possible, given the limits of
// pixel resolution, physical memory, etc.  Internally, this information is
// encoded using pixels at a range of resolutions, depending on how the
// boundaries of the area in question and the pixelization scheme interact.
// However, the goal of the class is to abstract away those details, allowing
// the user to treat Maps as a pure representative of spherical geometry.

#include "stomp_core.h"
#include "stomp_map.h"

namespace Stomp {

SubMap::SubMap(uint32_t superpixnum) {
  superpixnum_ = superpixnum;
  area_ = 0.0;
  size_ = 0;
  min_resolution_ = MaxPixelResolution;
  max_resolution_ = HPixResolution;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  Pixel::PixelBound(HPixResolution, superpixnum, lambda_min_,
		    lambda_max_, eta_min_, eta_max_);
  z_min_ = sin(lambda_min_*DegToRad);
  z_max_ = sin(lambda_max_*DegToRad);
  initialized_ = false;
  unsorted_ = false;

  for (uint16_t resolution=HPixResolution, i=0;
       i<ResolutionLevels;resolution*=2, i++) {
    pixel_count_[resolution] = 0;
  }
}

SubMap::~SubMap() {
  if (!pix_.empty()) pix_.clear();
  superpixnum_ = MaxSuperpixnum;
  initialized_ = false;
}

void SubMap::AddPixel(Pixel& pix) {
  // If our pixels are input in proper order, then we don't need to resolve
  // things down the line.  Provided that every input pixel comes after the
  // last pixel input, then we're assured that the list is sorted.
  if (!pix_.empty())
    if (!Pixel::LocalOrder(pix_[pix_.size()-1], pix)) unsorted_ = true;

  // If we're dealing with a sorted input list, then we can go ahead and
  // collect our summary statistics as we go.  Otherwise, we don't bother
  // since the results will be over-written by the Resolve method.
  if (!unsorted_) {
    area_ += pix.Area();
    if (pix.Resolution() < min_resolution_) min_resolution_ = pix.Resolution();
    if (pix.Resolution() > max_resolution_) max_resolution_ = pix.Resolution();
    if (pix.Weight() < min_weight_) min_weight_ = pix.Weight();
    if (pix.Weight() > max_weight_) max_weight_ = pix.Weight();
    pixel_count_[pix.Resolution()]++;
  }

  pix_.push_back(pix);
  size_ = pix_.size();
  initialized_ = true;
}

void SubMap::Resolve(bool force_resolve) {
  if (pix_.size() != size_) unsorted_ = true;

  if (unsorted_ || force_resolve) {
    Pixel::ResolveSuperPixel(pix_);

    area_ = 0.0;
    min_resolution_ = MaxPixelResolution;
    max_resolution_ = HPixResolution;
    min_weight_ = 1.0e30;
    max_weight_ = -1.0e30;
    for (uint16_t resolution=HPixResolution, i=0;
	 i<ResolutionLevels;resolution*=2, i++) {
      pixel_count_[resolution] = 0;
    }

    for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
      area_ += iter->Area();
      if (iter->Resolution() < min_resolution_)
        min_resolution_ = iter->Resolution();
      if (iter->Resolution() > max_resolution_)
        max_resolution_ = iter->Resolution();
      if (iter->Weight() < min_weight_) min_weight_ = iter->Weight();
      if (iter->Weight() > max_weight_) max_weight_ = iter->Weight();
      pixel_count_[iter->Resolution()]++;
    }
  }

  unsorted_ = false;
  size_ = pix_.size();
  if (pix_.size() > 0) initialized_ = true;
}

void SubMap::SetMinimumWeight(double min_weight) {
  PixelVector pix;
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    if (DoubleGE(iter->Weight(), min_weight)) pix.push_back(*iter);
  }

  Clear();
  for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) AddPixel(*iter);
  Resolve();
}

void SubMap::SetMaximumWeight(double max_weight) {
  PixelVector pix;
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    if (DoubleLE(iter->Weight(), max_weight)) pix.push_back(*iter);
  }

  Clear();
  for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) AddPixel(*iter);
  Resolve();
}

void SubMap::SetMaximumResolution(uint16_t max_resolution,
				  bool average_weights) {
  PixelVector pix;
  pix.reserve(Size());

  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    Pixel tmp_pix = *iter;
    if (average_weights) {
      if (iter->Resolution() > max_resolution) {
	tmp_pix.SetToSuperPix(max_resolution);
	tmp_pix.SetWeight(FindAverageWeight(tmp_pix));
      }
    } else {
      if (iter->Resolution() > max_resolution) {
	tmp_pix.SetToSuperPix(max_resolution);
	tmp_pix.SetWeight(FindUnmaskedFraction(tmp_pix));
      } else {
	tmp_pix.SetWeight(1.0);
      }
    }
    pix.push_back(tmp_pix);
  }

  Clear();
  for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) AddPixel(*iter);
  Resolve();
}

bool SubMap::FindLocation(AngularCoordinate& ang, double& weight) {
  bool keep = false;
  weight = -1.0e-30;

  for (uint16_t resolution=min_resolution_;
       resolution<=max_resolution_;resolution*=2) {
    Pixel tmp_pix(ang,resolution);
    PixelPair iter = equal_range(pix_.begin(), pix_.end(), tmp_pix,
                                 Pixel::SuperPixelBasedOrder);
    if (iter.first != iter.second) {
      keep = true;
      weight = iter.first->Weight();
    }
    if (keep) resolution = max_resolution_*2;
  }

  return keep;
}

double SubMap::FindUnmaskedFraction(Pixel& pix) {
  PixelIterator iter;
  if (pix.Resolution() == max_resolution_) {
    iter = pix_.end();
  } else {
    Pixel tmp_pix(pix.PixelX0()*2, pix.PixelY0()*2,
		  pix.Resolution()*2, 1.0);
    iter = lower_bound(pix_.begin(),pix_.end(),
                       tmp_pix,Pixel::SuperPixelBasedOrder);
  }

  uint16_t resolution = min_resolution_;
  double unmasked_fraction = 0.0;
  bool found_pixel = false;
  while (resolution <= pix.Resolution() && !found_pixel) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution);
    PixelPair super_iter = equal_range(pix_.begin(), iter, tmp_pix,
                                       Pixel::SuperPixelBasedOrder);
    if (super_iter.first != super_iter.second) {
      found_pixel = true;
      unmasked_fraction = 1.0;
    }
    resolution *= 2;
  }

  while (iter != pix_.end() && !found_pixel) {
    if (pix.Contains(*iter)) {
      double pixel_fraction =
          static_cast<double> (pix.Resolution()*pix.Resolution())/
          (iter->Resolution()*iter->Resolution());
      unmasked_fraction += pixel_fraction;
    }
    ++iter;
  }

  return unmasked_fraction;
}

int8_t SubMap::FindUnmaskedStatus(Pixel& pix) {
  PixelIterator iter;
  if (pix.Resolution() == max_resolution_) {
    iter = pix_.end();
  } else {
    Pixel tmp_pix(pix.PixelX0()*2, pix.PixelY0()*2,
		  pix.Resolution()*2, 1.0);
    iter = lower_bound(pix_.begin(), pix_.end(), tmp_pix,
                       Pixel::SuperPixelBasedOrder);
  }

  uint16_t resolution = min_resolution_;
  int8_t unmasked_status = 0;
  while ((resolution <= pix.Resolution()) && (unmasked_status == 0)) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution);
    PixelPair super_iter = equal_range(pix_.begin(),iter,tmp_pix,
                                       Pixel::SuperPixelBasedOrder);
    if (super_iter.first != super_iter.second) unmasked_status = 1;
    resolution *= 2;
  }

  while ((iter != pix_.end()) && (unmasked_status == 0)) {
    if (pix.Contains(*iter)) unmasked_status = -1;
    ++iter;
  }

  return unmasked_status;
}

double SubMap::FindAverageWeight(Pixel& pix) {
  PixelIterator iter;

  if (pix.Resolution() == max_resolution_) {
    iter = pix_.end();
  } else {
    Pixel tmp_pix(pix.Resolution()*2,0,pix.Superpixnum(),1.0);
    iter = lower_bound(pix_.begin(),pix_.end(),tmp_pix,
                       Pixel::SuperPixelBasedOrder);
  }

  double unmasked_fraction = 0.0, weighted_average = 0.0;
  bool found_pixel = false;
  uint16_t resolution = min_resolution_;
  while (resolution <= pix.Resolution() && !found_pixel) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution);
    PixelPair super_iter = equal_range(pix_.begin(), iter, tmp_pix,
                                       Pixel::SuperPixelBasedOrder);
    if (super_iter.first != super_iter.second) {
      found_pixel = true;
      weighted_average = super_iter.first->Weight();
      unmasked_fraction = 1.0;
    }
    resolution *= 2;
  }

  while (iter != pix_.end() && !found_pixel) {
    if (pix.Contains(*iter)) {
      double pixel_fraction =
	static_cast<double> (pix.Resolution()*pix.Resolution())/
	(iter->Resolution()*iter->Resolution());
      unmasked_fraction += pixel_fraction;
      weighted_average += iter->Weight()*pixel_fraction;
    }
    ++iter;
  }

  if (unmasked_fraction > 0.000000001) weighted_average /= unmasked_fraction;

  return weighted_average;
}

void SubMap::FindMatchingPixels(Pixel& pix, PixelVector& match_pix,
				bool use_local_weights) {
  if (!match_pix.empty()) match_pix.clear();

  bool found_pixel = false;
  PixelIterator iter, find_iter;
  if (pix.Resolution() == max_resolution_) {
    iter = pix_.end();
  } else {
    Pixel tmp_pix(pix.PixelX0()*2, pix.PixelY0()*2, pix.Resolution()*2, 1.0);
    iter = lower_bound(pix_.begin(),pix_.end(),tmp_pix,
                       Pixel::SuperPixelBasedOrder);
  }

  uint16_t resolution = min_resolution_;
  while (resolution <= pix.Resolution() && !found_pixel) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution);
    find_iter = lower_bound(pix_.begin(), iter, tmp_pix,
                            Pixel::SuperPixelBasedOrder);
    if (Pixel::PixelMatch(*find_iter,tmp_pix)) {
      found_pixel = true;
      tmp_pix = pix;
      if (use_local_weights) tmp_pix.SetWeight(find_iter->Weight());
      match_pix.push_back(tmp_pix);
    }
    resolution *= 2;
  }

  while (iter != pix_.end() && !found_pixel) {
    if (pix.Contains(*iter)) {
      Pixel tmp_pix = *iter;
      tmp_pix = *iter;
      if (!use_local_weights) tmp_pix.SetWeight(pix.Weight());
      match_pix.push_back(tmp_pix);
    }

    ++iter;
  }
}

double SubMap::AverageWeight() {
  double unmasked_fraction = 0.0, weighted_average = 0.0;
  if (initialized_) {
    for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
      weighted_average += iter->Weight()*iter->Area();
      unmasked_fraction += iter->Area();
    }
    weighted_average /= unmasked_fraction;
  }
  return weighted_average;
}

void SubMap::Soften(PixelVector& output_pix, uint16_t max_resolution,
		    bool average_weights) {
  if (!output_pix.empty()) output_pix.clear();
  output_pix.reserve(pix_.size());

  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    Pixel tmp_pix = *iter;
    if (average_weights) {
      if (iter->Resolution() > max_resolution) {
	tmp_pix.SetToSuperPix(max_resolution);
	tmp_pix.SetWeight(FindAverageWeight(tmp_pix));
      }
    } else {
      if (iter->Resolution() > max_resolution) {
	tmp_pix.SetToSuperPix(max_resolution);
	tmp_pix.SetWeight(FindUnmaskedFraction(tmp_pix));
      } else {
	tmp_pix.SetWeight(1.0);
      }
    }
    output_pix.push_back(tmp_pix);
  }

  Pixel::ResolveSuperPixel(output_pix);
}

bool SubMap::Add(Map& stomp_map, bool drop_single) {
  PixelVector keep_pix;
  PixelVector resolve_pix;

  // Iterate over all of our pixels and check each against the input map.
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    int8_t status = stomp_map.FindUnmaskedStatus(*iter);

    // If the pixel is completely outside the input map, then we can keep
    // the whole pixel, provided that we're not dropping area that's not in
    // both Maps.
    if (status == 0 && !drop_single) {
      keep_pix.push_back(*iter);
    }

    // If the pixel is completely inside of the input map, then we add it to
    // the keep array.  There's a complication here in that we need to account
    // for cases where the current pixel may be covering multiple pixels with
    // different weights in the input map.
    if (status == 1) {
      PixelVector match_pix;
      stomp_map.FindMatchingPixels(*iter, match_pix);
      for (PixelIterator match_iter=match_pix.begin();
	   match_iter!=match_pix.end();++match_iter) {
	match_iter->SetWeight(match_iter->Weight() + iter->Weight());
	keep_pix.push_back(*match_iter);
      }
    }

    // If there's partial overlap, then we need to refine the pixel and check
    // again until we find the parts that don't overlap.
    if (status == -1) {
      PixelVector sub_pix;
      iter->SubPix(2*iter->Resolution(), sub_pix);
      for (PixelIterator sub_iter=sub_pix.begin();
	   sub_iter!=sub_pix.end();++sub_iter) {
	sub_iter->SetWeight(iter->Weight());
	resolve_pix.push_back(*sub_iter);
      }
    }
  }

  // Now we check those pixels that need resolving and iterate until there are
  // no more to check.
  while (resolve_pix.size() > 0) {
    PixelVector tmp_pix;
    tmp_pix.reserve(resolve_pix.size());
    for (PixelIterator iter=resolve_pix.begin();
	 iter!=resolve_pix.end();++iter) tmp_pix.push_back(*iter);

    resolve_pix.clear();

    for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter) {
      int8_t status = stomp_map.FindUnmaskedStatus(*iter);

      if (status == 0 && !drop_single) {
	keep_pix.push_back(*iter);
      }

      if (status == 1) {
	PixelVector match_pix;
	stomp_map.FindMatchingPixels(*iter, match_pix);
	for (PixelIterator match_iter=match_pix.begin();
	     match_iter!=match_pix.end();++match_iter) {
	  match_iter->SetWeight(match_iter->Weight() + iter->Weight());
	  keep_pix.push_back(*match_iter);
	}
      }

      if (status == -1) {
	PixelVector sub_pix;
	iter->SubPix(2*iter->Resolution(), sub_pix);
	for (PixelIterator sub_iter=sub_pix.begin();
	     sub_iter!=sub_pix.end();++sub_iter) {
	  sub_iter->SetWeight(iter->Weight());
	  resolve_pix.push_back(*sub_iter);
	}
      }
    }
  }

  // That covers the area in our current map.  However, if we're not dropping
  // area that is contained in just one map, we need to find those pixels in
  // the input Map that didn't overlap with anything in our current Map and
  // add those to the array of keep pixels.
  if (!drop_single) {
    PixelVector stomp_pix;
    stomp_map.Pixels(stomp_pix, Superpixnum());

    for (PixelIterator iter=stomp_pix.begin();iter!=stomp_pix.end();++iter) {
      // Only keep those pixels that are completely outside of our current Map.
      if (FindUnmaskedStatus(*iter)) keep_pix.push_back(*iter);
    }
  }

  // Now we clear out our current set of pixels and replace them by the ones
  // that weren't in the input map.
  Clear();

  for (PixelIterator iter=keep_pix.begin();iter!=keep_pix.end();++iter) {
    AddPixel(*iter);
  }

  if (unsorted_) Resolve();

  return true;
}

bool SubMap::Multiply(Map& stomp_map, bool drop_single) {
  PixelVector keep_pix;
  PixelVector resolve_pix;

  // Iterate over all of our pixels and check each against the input map.
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    int8_t status = stomp_map.FindUnmaskedStatus(*iter);

    // If the pixel is completely outside the input map, then we can keep
    // the whole pixel, provided that we're not dropping area that's not in
    // both Maps.
    if (status == 0 && !drop_single) {
      keep_pix.push_back(*iter);
    }

    // If the pixel is completely inside of the input map, then we add it to
    // the keep array.  There's a complication here in that we need to account
    // for cases where the current pixel may be covering multiple pixels with
    // different weights in the input map.
    if (status == 1) {
      PixelVector match_pix;
      stomp_map.FindMatchingPixels(*iter, match_pix);
      for (PixelIterator match_iter=match_pix.begin();
	   match_iter!=match_pix.end();++match_iter) {
	match_iter->SetWeight(match_iter->Weight()*iter->Weight());
	keep_pix.push_back(*match_iter);
      }
    }

    // If there's partial overlap, then we need to refine the pixel and check
    // again until we find the parts that don't overlap.
    if (status == -1) {
      PixelVector sub_pix;
      iter->SubPix(2*iter->Resolution(), sub_pix);
      for (PixelIterator sub_iter=sub_pix.begin();
	   sub_iter!=sub_pix.end();++sub_iter) {
	sub_iter->SetWeight(iter->Weight());
	resolve_pix.push_back(*sub_iter);
      }
    }
  }

  // Now we check those pixels that need resolving and iterate until there are
  // no more to check.
  while (resolve_pix.size() > 0) {
    PixelVector tmp_pix;
    tmp_pix.reserve(resolve_pix.size());
    for (PixelIterator iter=resolve_pix.begin();
	 iter!=resolve_pix.end();++iter) tmp_pix.push_back(*iter);

    resolve_pix.clear();

    for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter) {
      int8_t status = stomp_map.FindUnmaskedStatus(*iter);

      if (status == 0 && !drop_single) {
	keep_pix.push_back(*iter);
      }

      if (status == 1) {
	PixelVector match_pix;
	stomp_map.FindMatchingPixels(*iter, match_pix);
	for (PixelIterator match_iter=match_pix.begin();
	     match_iter!=match_pix.end();++match_iter) {
	  match_iter->SetWeight(match_iter->Weight()*iter->Weight());
	  keep_pix.push_back(*match_iter);
	}
      }

      if (status == -1) {
	PixelVector sub_pix;
	iter->SubPix(2*iter->Resolution(), sub_pix);
	for (PixelIterator sub_iter=sub_pix.begin();
	     sub_iter!=sub_pix.end();++sub_iter) {
	  sub_iter->SetWeight(iter->Weight());
	  resolve_pix.push_back(*sub_iter);
	}
      }
    }
  }

  // That covers the area in our current map.  However, if we're not dropping
  // area that is contained in just one map, we need to find those pixels in
  // the input Map that didn't overlap with anything in our current Map and
  // add those to the array of keep pixels.
  if (!drop_single) {
    PixelVector stomp_pix;
    stomp_map.Pixels(stomp_pix, Superpixnum());

    for (PixelIterator iter=stomp_pix.begin();iter!=stomp_pix.end();++iter) {
      // Only keep those pixels that are completely outside of our current Map.
      if (FindUnmaskedStatus(*iter)) keep_pix.push_back(*iter);
    }
  }

  // Now we clear out our current set of pixels and replace them by the ones
  // that weren't in the input map.
  Clear();

  for (PixelIterator iter=keep_pix.begin();iter!=keep_pix.end();++iter) {
    AddPixel(*iter);
  }

  if (unsorted_) Resolve();

  return true;
}

bool SubMap::Exclude(Map& stomp_map) {

  PixelVector keep_pix;
  PixelVector resolve_pix;

  // Iterate over all of our pixels and check each against the input map.
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    int8_t status = stomp_map.FindUnmaskedStatus(*iter);

    // If the pixel is completely outside the input map, then we can keep
    // the whole pixel.
    if (status == 0) {
      keep_pix.push_back(*iter);
    }

    // If there's partial overlap, then we need to refine the pixel and check
    // again until we find the parts that don't overlap.
    if (status == -1) {
      PixelVector sub_pix;
      iter->SubPix(2*iter->Resolution(), sub_pix);
      for (PixelIterator sub_iter=sub_pix.begin();
	   sub_iter!=sub_pix.end();++sub_iter) {
	sub_iter->SetWeight(iter->Weight());
	resolve_pix.push_back(*sub_iter);
      }
    }
  }

  // Now we check those pixels that need resolving and iterate until there are
  // no more to check.
  while (resolve_pix.size() > 0) {
    PixelVector tmp_pix;
    tmp_pix.reserve(resolve_pix.size());
    for (PixelIterator iter=resolve_pix.begin();
	 iter!=resolve_pix.end();++iter) tmp_pix.push_back(*iter);

    resolve_pix.clear();

    for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter) {
      int8_t status = stomp_map.FindUnmaskedStatus(*iter);

      if (status == 0) {
	keep_pix.push_back(*iter);
      }

      if (status == -1) {
	PixelVector sub_pix;
	iter->SubPix(2*iter->Resolution(), sub_pix);
	for (PixelIterator sub_iter=sub_pix.begin();
	     sub_iter!=sub_pix.end();++sub_iter) {
	  sub_iter->SetWeight(iter->Weight());
	  resolve_pix.push_back(*sub_iter);
	}
      }
    }
  }

  // Now we clear out our current set of pixels and replace them by the ones
  // that weren't in the input map.
  Clear();

  for (PixelIterator iter=keep_pix.begin();iter!=keep_pix.end();++iter) {
    AddPixel(*iter);
  }

  if (unsorted_) Resolve();

  return true;
}

void SubMap::ScaleWeight(const double weight_scale) {
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter)
    iter->SetWeight(iter->Weight()*weight_scale);
}

void SubMap::AddConstantWeight(const double add_weight) {
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter)
    iter->SetWeight(iter->Weight()+add_weight);
}

void SubMap::InvertWeight() {
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    if ((iter->Weight() > 1.0e-15) || (iter->Weight() < -1.0e-15)) {
      iter->SetWeight(1.0/iter->Weight());
    } else {
      iter->SetWeight(0.0);
    }
  }
}

void SubMap::Pixels(PixelVector& pix) {
  if (!pix.empty()) pix.clear();
  pix.reserve(Size());
  for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter)
    pix.push_back(*iter);
}

void SubMap::Clear() {
  area_ = 0.0;
  size_ = 0;
  min_resolution_ = MaxPixelResolution;
  max_resolution_ = HPixResolution;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  if (!pix_.empty()) pix_.clear();
  initialized_ = false;
  unsorted_ = false;
}

uint32_t SubMap::Superpixnum() {
  return superpixnum_;
}

PixelIterator SubMap::Begin() {
  return (initialized_ ? pix_.begin() : pix_.end());
}

PixelIterator SubMap::End() {
  return pix_.end();
}

double SubMap::Area() {
  return area_;
}

bool SubMap::Initialized() {
  return initialized_;
}

bool SubMap::Unsorted() {
  return unsorted_;
}

uint16_t SubMap::MinResolution() {
  return min_resolution_;
}

uint16_t SubMap::MaxResolution() {
  return max_resolution_;
}

double SubMap::MinWeight() {
  return min_weight_;
}

double SubMap::MaxWeight() {
  return max_weight_;
}

double SubMap::LambdaMin() {
  return lambda_min_;
}

double SubMap::LambdaMax() {
  return lambda_max_;
}

double SubMap::EtaMin() {
  return eta_min_;
}

double SubMap::EtaMax() {
  return eta_max_;
}

double SubMap::ZMin() {
  return z_min_;
}

double SubMap::ZMax() {
  return z_max_;
}

uint32_t SubMap::Size() {
  return size_;
}

uint32_t SubMap::PixelCount(uint16_t resolution) {
  return (!(resolution % 2) ? pixel_count_[resolution] : 0);
}

Map::Map() {
  area_ = 0.0;
  size_ = 0;
  min_resolution_ = MaxPixelResolution;
  max_resolution_ = HPixResolution;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();

  for (uint16_t resolution=HPixResolution, i=0;
       i<ResolutionLevels;resolution*=2, i++)
    pixel_count_[resolution] = 0;

  sub_map_.reserve(MaxSuperpixnum);

  for (uint32_t k=0;k<MaxSuperpixnum;k++)
    sub_map_.push_back(SubMap(k));

  begin_ = MapIterator(0, sub_map_[0].Begin());
  end_ = begin_;
}

Map::Map(PixelVector& pix, bool force_resolve) {
  area_ = 0.0;
  size_ = 0;
  min_resolution_ = MaxPixelResolution;
  max_resolution_ = HPixResolution;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();
  for (uint16_t resolution=HPixResolution, i=0;
       i<ResolutionLevels;resolution*=2, i++) {
    pixel_count_[resolution] = 0;
  }

  sub_map_.reserve(MaxSuperpixnum);

  for (uint32_t k=0;k<MaxSuperpixnum;k++)
    sub_map_.push_back(SubMap(k));

  for (PixelIterator iter=pix.begin();iter!=pix.end();++iter)
    sub_map_[iter->Superpixnum()].AddPixel(*iter);

  bool found_beginning = false;
  for (SubMapIterator iter=sub_map_.begin();iter!=sub_map_.end();++iter) {
    if (iter->Initialized()) {
      if (iter->Unsorted() || force_resolve) iter->Resolve(force_resolve);

      if (!found_beginning) {
	begin_ = MapIterator(iter->Superpixnum(), iter->Begin());
	found_beginning = true;
      }
      end_ = MapIterator(iter->Superpixnum(), iter->End());

      area_ += iter->Area();
      size_ += iter->Size();
      if (min_resolution_ > iter->MinResolution())
        min_resolution_ = iter->MinResolution();
      if (max_resolution_ < iter->MaxResolution())
        max_resolution_ = iter->MaxResolution();
      if (iter->MinWeight() < min_weight_) min_weight_ = iter->MinWeight();
      if (iter->MaxWeight() > max_weight_) max_weight_ = iter->MaxWeight();
      for (uint16_t resolution=HPixResolution, i=0;
	   i<ResolutionLevels;resolution*=2, i++) {
	pixel_count_[resolution] += iter->PixelCount(resolution);
      }
    }
  }
}

Map::Map(std::string& InputFile, bool hpixel_format, bool weighted_map) {
  Read(InputFile, hpixel_format, weighted_map);
}

Map::~Map() {
  min_resolution_ = max_resolution_ = 0;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  Clear();
}

bool Map::Initialize() {
  area_ = 0.0;
  size_ = 0;
  min_resolution_ = MaxPixelResolution;
  max_resolution_ = HPixResolution;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();
  for (uint16_t resolution=HPixResolution, i=0;
       i<ResolutionLevels;resolution*=2, i++) {
    pixel_count_[resolution] = 0;
  }

  bool found_valid_superpixel = false;
  bool found_beginning = false;
  for (SubMapIterator iter=sub_map_.begin();iter!=sub_map_.end();++iter) {
    if (iter->Initialized()) {
      found_valid_superpixel = true;
      if (iter->Unsorted()) iter->Resolve();

      if (!found_beginning) {
	begin_ = MapIterator(iter->Superpixnum(), iter->Begin());
	found_beginning = true;
      }
      end_ = MapIterator(iter->Superpixnum(), iter->End());

      area_ += iter->Area();
      size_ += iter->Size();
      if (min_resolution_ > iter->MinResolution())
        min_resolution_ = iter->MinResolution();
      if (max_resolution_ < iter->MaxResolution())
        max_resolution_ = iter->MaxResolution();
      if (iter->MinWeight() < min_weight_) min_weight_ = iter->MinWeight();
      if (iter->MaxWeight() > max_weight_) max_weight_ = iter->MaxWeight();
      for (uint16_t resolution=HPixResolution, i=0;
	   i<ResolutionLevels;resolution*=2, i++) {
	pixel_count_[resolution] += iter->PixelCount(resolution);
      }
    }
  }

  return found_valid_superpixel;
}

bool Map::Initialize(PixelVector& pix, bool force_resolve) {
  area_ = 0.0;
  size_ = 0;
  min_resolution_ = MaxPixelResolution;
  max_resolution_ = HPixResolution;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();
  for (uint16_t resolution=HPixResolution, i=0;
       i<ResolutionLevels;resolution*=2, i++) {
    pixel_count_[resolution] = 0;
  }

  if (!sub_map_.empty()) {
    for (SubMapIterator iter=sub_map_.begin();iter!=sub_map_.end();++iter)
      iter->Clear();
    sub_map_.clear();
  }

  sub_map_.reserve(MaxSuperpixnum);

  for (uint32_t k=0;k<MaxSuperpixnum;k++) {
    SubMap tmp_sub_map(k);
    sub_map_.push_back(tmp_sub_map);
  }

  for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) {
    uint32_t k = iter->Superpixnum();

    sub_map_[k].AddPixel(*iter);
  }

  bool found_valid_superpixel = false;
  bool found_beginning = false;
  for (SubMapIterator iter=sub_map_.begin();iter!=sub_map_.end();++iter) {
    if (iter->Initialized()) {
      if (iter->Unsorted() || force_resolve) iter->Resolve(force_resolve);

      found_valid_superpixel = true;

      if (!found_beginning) {
	begin_ = MapIterator(iter->Superpixnum(), iter->Begin());
	found_beginning = true;
      }
      end_ = MapIterator(iter->Superpixnum(), iter->End());

      area_ += iter->Area();
      size_ += iter->Size();
      if (min_resolution_ > iter->MinResolution())
        min_resolution_ = iter->MinResolution();
      if (max_resolution_ < iter->MaxResolution())
        max_resolution_ = iter->MaxResolution();
      if (iter->MinWeight() < min_weight_) min_weight_ = iter->MinWeight();
      if (iter->MaxWeight() > max_weight_) max_weight_ = iter->MaxWeight();
      for (uint16_t resolution=HPixResolution, i=0;
	   i<ResolutionLevels;resolution*=2, i++) {
	pixel_count_[resolution] += iter->PixelCount(resolution);
      }
    }
  }

  return found_valid_superpixel;
}

void Map::Coverage(PixelVector& superpix, uint16_t resolution) {
  if (!superpix.empty()) superpix.clear();

  if (resolution == HPixResolution) {
    // If we're dealing with a coverage map at superpixel resolution (the
    // default behavior), then this is easy.  Just iterate over the submaps
    // and keep those that have been initialized.
    for (uint32_t k=0;k<MaxSuperpixnum;k++) {
      if (sub_map_[k].Initialized()) {
	// We store the unmasked fraction of each superpixel in the weight
	// value in case that's useful.
	Pixel tmp_pix(HPixResolution, k,
		      sub_map_[k].Area()/HPixArea);
	superpix.push_back(tmp_pix);
      }
    }
  } else {
    for (uint32_t k=0;k<MaxSuperpixnum;k++) {
      if (sub_map_[k].Initialized()) {
	Pixel tmp_pix(HPixResolution, k, 1.0);

	PixelVector sub_pix;
	tmp_pix.SubPix(resolution, sub_pix);
	for (PixelIterator iter=sub_pix.begin();iter!=sub_pix.end();++iter) {
	  // For each of the pixels in the superpixel, we check its status
	  // against the current map.  This is faster than finding the unmasked
	  // fraction directly and immediately tells us which pixels we can
	  // eliminate and which of those we do keep require further
	  // calculations to find the unmasked fraction.
	  int8_t unmasked_status = FindUnmaskedStatus(*iter);
	  if (unmasked_status != 0) {
	    if (unmasked_status == 1) {
	      iter->SetWeight(1.0);
	    } else {
	      iter->SetWeight(FindUnmaskedFraction(*iter));
	    }
	    superpix.push_back(*iter);
	  }
	}
      }
    }
    sort(superpix.begin(), superpix.end(), Pixel::SuperPixelBasedOrder);
  }
}

bool Map::Covering(Map& stomp_map, uint32_t maximum_pixels) {
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

    FindAverageWeight(pix);
    stomp_map.Initialize(pix);
  } else {
    // Ok, in this case, we can definitely produce a map that has at most
    // maximum_pixels in it.
    met_pixel_requirement = true;

    // To possibly save us as much work as possible, we check to see if
    // maximum_pixels is larger than our current set.  If so, we just return
    // a copy of the current map.
    if (maximum_pixels > Size()) {
      pix.clear();
      pix.reserve(Size());
      Pixels(pix);
      stomp_map.Initialize(pix);
    } else {
      // Ok, in this case, we have to do actual work.  As a rough rule of
      // thumb, we should expect that the number of pixels in any given
      // resolution level would double if we were to combine all of the pixels
      // at finer resolution into the current level.  So, we proceed from
      // coarse to fine until adding twice a given level would take us over
      // the maximum_pixels limit.  Then we re-sample all of the pixels below
      // that level and check that against our limit, iterating again at a
      // coarser resolution limit if necessary.  This should minimize
      // the number of times that we need to re-create the pixel map.
      uint16_t maximum_resolution = HPixResolution;
      uint32_t reduced_map_pixels = pixel_count_[maximum_resolution];

      while (reduced_map_pixels + 2*pixel_count_[2*maximum_resolution] <
	     maximum_pixels &&
	     maximum_resolution < MaxPixelResolution/2 ) {
	reduced_map_pixels += pixel_count_[maximum_resolution];
	maximum_resolution *= 2;
      }

      reduced_map_pixels = maximum_pixels + 1;
      while ((reduced_map_pixels > maximum_pixels) &&
	     (maximum_resolution >= HPixResolution)) {
	pix.clear();
	pix.reserve(Size());
	Pixels(pix);

	for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) {
	  if (iter->Resolution() > maximum_resolution) {
	    iter->SuperPix(maximum_resolution);
	  }
	}

	// Resolve the array of pixels ignoring their current weight values.
	// The consequence of this is that we may end up under-sampling any
	// scalar field that's been encoded onto the current map and that we
	// won't be producing a map that contains the most resolution possible
	// given the pixel limit.  However, given the usage, these are probably
	// not show-stoppers.
	Pixel::ResolvePixel(pix, true);

	reduced_map_pixels = pix.size();

	if (reduced_map_pixels < maximum_pixels) {
	  FindAverageWeight(pix);
	  stomp_map.Initialize(pix);
	} else {
	  maximum_resolution /= 2;
	}
      }
    }
  }

  return met_pixel_requirement;
}

void Map::Soften(Map& stomp_map, uint16_t maximum_resolution,
		 bool average_weights) {
  if (!stomp_map.Empty()) stomp_map.Clear();

  PixelVector pix;
  for (uint32_t k=0;k<MaxSuperpixnum;k++) {
    if (sub_map_[k].Initialized()) {
      PixelVector tmp_pix;
      sub_map_[k].Soften(tmp_pix, maximum_resolution, average_weights);

      for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter)
	pix.push_back(*iter);
    }
  }

  stomp_map.Initialize(pix);
}

void Map::Soften(uint16_t maximum_resolution, bool average_weights) {
  for (uint32_t k=0;k<MaxSuperpixnum;k++) {
    if (sub_map_[k].Initialized()) {
      sub_map_[k].SetMaximumResolution(maximum_resolution, average_weights);
    }
  }
  Initialize();
}

void Map::SetMinimumWeight(double min_weight) {
  for (uint32_t k=0;k<MaxSuperpixnum;k++) {
    if (sub_map_[k].Initialized()) {
      sub_map_[k].SetMinimumWeight(min_weight);
    }
  }
  Initialize();
}

void Map::SetMaximumWeight(double max_weight) {
  for (uint32_t k=0;k<MaxSuperpixnum;k++) {
    if (sub_map_[k].Initialized()) {
      sub_map_[k].SetMaximumWeight(max_weight);
    }
  }
  Initialize();
}

bool Map::RegionOnlyMap(int16_t region_index, Map& stomp_map) {
  if (!stomp_map.Empty()) stomp_map.Clear();

  bool map_generation_success = false;
  if ((region_index >= 0) && (region_index < NRegion())) {
    // First, we need to create a copy of our current map.
    PixelVector pix;
    Pixels(pix);
    stomp_map.Initialize(pix);
    pix.clear();

    // Now we generate a PixelVector for the input region.
    PixelVector region_pix;
    RegionArea(region_index, region_pix);

    // And finally, we find the intersection of these two maps.  This should
    // just be the part of our current map in the specified region.
    map_generation_success = stomp_map.IntersectMap(region_pix);
    region_pix.clear();
  }

  return map_generation_success;
}

bool Map::RegionExcludedMap(int16_t region_index, Map& stomp_map) {
  if (!stomp_map.Empty()) stomp_map.Clear();

  bool map_generation_success = false;
  if ((region_index >= 0) && (region_index < NRegion())) {
    // First, we need to create a copy of our current map.
    PixelVector pix;
    Pixels(pix);
    stomp_map.Initialize(pix);
    pix.clear();

    // Now we generate a Map for the input region.
    PixelVector region_pix;
    RegionArea(region_index, region_pix);

    // And finally, we exclude the region map from the map copy.  This should
    // just be the part of our current map outside the specified region.
    map_generation_success = stomp_map.ExcludeMap(region_pix);
    region_pix.clear();
  }

  return map_generation_success;
}

bool Map::FindLocation(AngularCoordinate& ang) {
  bool keep = false;
  double weight;

  uint32_t k;
  Pixel::Ang2Pix(HPixResolution, ang, k);

  if (sub_map_[k].Initialized()) keep = sub_map_[k].FindLocation(ang, weight);

  return keep;
}

bool Map::FindLocation(AngularCoordinate& ang, double& weight) {
  bool keep = false;

  uint32_t k;
  Pixel::Ang2Pix(HPixResolution, ang, k);

  if (sub_map_[k].Initialized()) keep = sub_map_[k].FindLocation(ang, weight);

  return keep;
}

double Map::FindLocationWeight(AngularCoordinate& ang) {
  bool keep = false;
  double weight = -1.0e-30;

  uint32_t k;
  Pixel::Ang2Pix(HPixResolution,ang,k);

  if (sub_map_[k].Initialized()) keep = sub_map_[k].FindLocation(ang, weight);

  return weight;
}

bool Map::Contains(Pixel& pix) {
  return (FindUnmaskedStatus(pix) == 1 ? true : false);
}

bool Map::Contains(Map& stomp_map) {
  return (FindUnmaskedStatus(stomp_map) == 1 ? true : false);
}

double Map::FindUnmaskedFraction(Pixel& pix) {
  double unmasked_fraction = 0.0;

  uint32_t k = pix.Superpixnum();

  if (sub_map_[k].Initialized())
    unmasked_fraction = sub_map_[k].FindUnmaskedFraction(pix);

  return unmasked_fraction;
}

void Map::FindUnmaskedFraction(PixelVector& pix,
			       std::vector<double>& unmasked_fraction) {

  if (!unmasked_fraction.empty()) unmasked_fraction.clear();

  unmasked_fraction.reserve(pix.size());

  for (uint32_t i=0;i<pix.size();i++){
    double pixel_unmasked_fraction = 0.0;
    uint32_t k = pix[i].Superpixnum();

    if (sub_map_[k].Initialized())
      pixel_unmasked_fraction = sub_map_[k].FindUnmaskedFraction(pix[i]);

    unmasked_fraction.push_back(pixel_unmasked_fraction);
  }
}

void Map::FindUnmaskedFraction(PixelVector& pix) {
  for (uint32_t i=0;i<pix.size();i++){
    double pixel_unmasked_fraction = 0.0;
    uint32_t k = pix[i].Superpixnum();

    if (sub_map_[k].Initialized())
      pixel_unmasked_fraction = sub_map_[k].FindUnmaskedFraction(pix[i]);
    pix[i].SetWeight(pixel_unmasked_fraction);
  }
}

double Map::FindUnmaskedFraction(Map& stomp_map) {
  double total_unmasked_area = 0.0;
  for (MapIterator iter=stomp_map.Begin();
       iter!=stomp_map.End();stomp_map.Iterate(&iter)) {
    double pixel_unmasked_fraction = 0.0;
    uint32_t k = iter.second->Superpixnum();

    if (sub_map_[k].Initialized())
      pixel_unmasked_fraction =
	sub_map_[k].FindUnmaskedFraction(*(iter.second));
    total_unmasked_area += pixel_unmasked_fraction*iter.second->Area();
  }

  return total_unmasked_area/stomp_map.Area();
}

int8_t Map::FindUnmaskedStatus(Pixel& pix) {
  int8_t unmasked_status = 0;

  uint32_t k = pix.Superpixnum();

  if (sub_map_[k].Initialized())
    unmasked_status = sub_map_[k].FindUnmaskedStatus(pix);

  return unmasked_status;
}

void Map::FindUnmaskedStatus(PixelVector& pix,
			     std::vector<int8_t>& unmasked_status) {

  if (!unmasked_status.empty()) unmasked_status.clear();
  unmasked_status.reserve(pix.size());

  for (uint32_t i=0;i<pix.size();i++){
    int8_t pixel_unmasked_status = 0;
    uint32_t k = pix[i].Superpixnum();

    if (sub_map_[k].Initialized())
      pixel_unmasked_status = sub_map_[k].FindUnmaskedStatus(pix[i]);

    unmasked_status.push_back(pixel_unmasked_status);
  }
}

int8_t Map::FindUnmaskedStatus(Map& stomp_map) {
  // Get the status for the first pixel in the map to seed our initial status.
  MapIterator iter=stomp_map.Begin();
  int8_t map_unmasked_status = 0;
  uint32_t k = iter.second->Superpixnum();

  if (sub_map_[k].Initialized())
    map_unmasked_status = sub_map_[k].FindUnmaskedStatus(*(iter.second));

  stomp_map.Iterate(&iter);

  // Now iterate over the rest of the map to figure out the global status.
  for (;iter!=stomp_map.End();stomp_map.Iterate(&iter)) {
    int8_t unmasked_status = 0;
    k = iter.second->Superpixnum();

    if (sub_map_[k].Initialized())
      unmasked_status = sub_map_[k].FindUnmaskedStatus(*(iter.second));

    if (map_unmasked_status == 1) {
      // If we currently thought that the input Map was completely inside of our
      // Map, but find that this Pixel is either outside the Map or only
      // partially inside the Map, then we switch the global status to partial.
      if (unmasked_status == 0 || unmasked_status == -1)
	map_unmasked_status = -1;
    }

    if (map_unmasked_status == 0) {
      // If we currently thought that the input Map was completely outside the
      // Map, but find that this Pixel is either fully or partially contained
      // in the Map, then we switch the global status to partial.
      if (unmasked_status == 1 || unmasked_status == -1)
	map_unmasked_status = -1;
    }

    if (map_unmasked_status == -1) {
      // If we find that the input Map's unmasked status is partial, then no
      // further testing will change that status.  At this point, we can break
      // our loop and return the global status.
      break;
    }
  }

  return map_unmasked_status;
}

double Map::FindAverageWeight(Pixel& pix) {
  double weighted_average = 0.0;
  uint32_t k = pix.Superpixnum();

  if (sub_map_[k].Initialized())
    weighted_average = sub_map_[k].FindAverageWeight(pix);

  return weighted_average;
}

void Map::FindAverageWeight(PixelVector& pix,
			    std::vector<double>& weighted_average) {
  if (!weighted_average.empty()) weighted_average.clear();

  for (uint32_t i=0;i<pix.size();i++) {
    double pixel_weighted_average = 0.0;
    uint32_t k = pix[i].Superpixnum();

    if (sub_map_[k].Initialized())
      pixel_weighted_average = sub_map_[k].FindAverageWeight(pix[i]);

    weighted_average.push_back(pixel_weighted_average);
  }
}

void Map::FindAverageWeight(PixelVector& pix) {
  for (uint32_t i=0;i<pix.size();i++) {
    double pixel_weighted_average = 0.0;
    uint32_t k = pix[i].Superpixnum();

    if (sub_map_[k].Initialized())
      pixel_weighted_average = sub_map_[k].FindAverageWeight(pix[i]);
    pix[i].SetWeight(pixel_weighted_average);
  }
}

void Map::FindMatchingPixels(Pixel& pix, PixelVector& match_pix,
			     bool use_local_weights) {
  if (!match_pix.empty()) match_pix.clear();

  uint32_t k = pix.Superpixnum();

  if (sub_map_[k].Initialized()) {
    PixelVector tmp_pix;

    sub_map_[k].FindMatchingPixels(pix,tmp_pix,use_local_weights);

    if (!tmp_pix.empty())
      for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter)
        match_pix.push_back(*iter);
  }
}

void Map::FindMatchingPixels(PixelVector& pix, PixelVector& match_pix,
			     bool use_local_weights) {
  if (!match_pix.empty()) match_pix.clear();

  for (uint32_t i=0;i<pix.size();i++) {

    uint32_t k = pix[i].Superpixnum();

    if (sub_map_[k].Initialized()) {
      PixelVector tmp_pix;

      sub_map_[k].FindMatchingPixels(pix[i],tmp_pix,use_local_weights);

      if (!tmp_pix.empty())
        for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter)
          match_pix.push_back(*iter);
    }
  }
}

void Map::GenerateRandomPoints(AngularVector& ang, uint32_t n_point,
			       bool use_weighted_sampling) {
  if (!ang.empty()) ang.clear();
  ang.reserve(n_point);

  double minimum_probability = -0.0001;
  double probability_slope = 0.0;

  if (use_weighted_sampling) {
    if (max_weight_ - min_weight_ < 0.0001) {
      use_weighted_sampling = false;
    } else {
      minimum_probability = 1.0/(max_weight_ - min_weight_ + 1.0);
      probability_slope =
          (1.0 - minimum_probability)/(max_weight_ - min_weight_);
    }
  }

  PixelVector superpix;
  Coverage(superpix);

  MTRand mtrand;
  mtrand.seed();

  for (uint32_t m=0;m<n_point;m++) {
    bool keep = false;
    double lambda, eta, z, weight, probability_limit;
    AngularCoordinate tmp_ang(0.0,0.0);
    uint32_t n,k;

    while (!keep) {
      n = mtrand.randInt(superpix.size()-1);
      k = superpix[n].Superpixnum();

      z = sub_map_[k].ZMin() + mtrand.rand(sub_map_[k].ZMax() -
					   sub_map_[k].ZMin());
      lambda = asin(z)*RadToDeg;
      eta = sub_map_[k].EtaMin() + mtrand.rand(sub_map_[k].EtaMax() -
					       sub_map_[k].EtaMin());
      tmp_ang.SetSurveyCoordinates(lambda,eta);

      keep = sub_map_[k].FindLocation(tmp_ang,weight);

      if (use_weighted_sampling && keep) {
        probability_limit =
            minimum_probability + (weight - min_weight_)*probability_slope;
        if (mtrand.rand(1.0) > probability_limit) keep = false;
      }
    }

    ang.push_back(tmp_ang);
  }
}

void Map::GenerateRandomPoints(WAngularVector& ang, WAngularVector& input_ang) {
  if (!ang.empty()) ang.clear();
  ang.reserve(input_ang.size());

  PixelVector superpix;
  Coverage(superpix);

  MTRand mtrand;
  mtrand.seed();

  WeightedAngularCoordinate tmp_ang;
  for (uint32_t m=0;m<input_ang.size();m++) {
    bool keep = false;
    double lambda, eta, z, map_weight;
    uint32_t n,k;

    while (!keep) {
      n = mtrand.randInt(superpix.size()-1);
      k = superpix[n].Superpixnum();

      z = sub_map_[k].ZMin() + mtrand.rand(sub_map_[k].ZMax() -
					   sub_map_[k].ZMin());
      lambda = asin(z)*RadToDeg;
      eta = sub_map_[k].EtaMin() + mtrand.rand(sub_map_[k].EtaMax() -
					       sub_map_[k].EtaMin());
      tmp_ang.SetSurveyCoordinates(lambda,eta);

      keep = sub_map_[k].FindLocation(tmp_ang, map_weight);
    }
    tmp_ang.SetWeight(input_ang[m].Weight());
    ang.push_back(tmp_ang);
  }
}

void Map::GenerateRandomPoints(WAngularVector& ang,
			       std::vector<double>& weights) {
  if (!ang.empty()) ang.clear();
  ang.reserve(weights.size());

  PixelVector superpix;
  Coverage(superpix);

  MTRand mtrand;
  mtrand.seed();

  WeightedAngularCoordinate tmp_ang;
  for (uint32_t m=0;m<weights.size();m++) {
    bool keep = false;
    double lambda, eta, z, map_weight;
    uint32_t n,k;

    while (!keep) {
      n = mtrand.randInt(superpix.size()-1);
      k = superpix[n].Superpixnum();

      z = sub_map_[k].ZMin() + mtrand.rand(sub_map_[k].ZMax() -
					   sub_map_[k].ZMin());
      lambda = asin(z)*RadToDeg;
      eta = sub_map_[k].EtaMin() + mtrand.rand(sub_map_[k].EtaMax() -
					       sub_map_[k].EtaMin());
      tmp_ang.SetSurveyCoordinates(lambda,eta);

      keep = sub_map_[k].FindLocation(tmp_ang, map_weight);
    }
    tmp_ang.SetWeight(weights[m]);
    ang.push_back(tmp_ang);
  }
}

bool Map::Write(std::string& OutputFile, bool hpixel_format,
		bool weighted_map) {
  std::ofstream output_file(OutputFile.c_str());

  if (output_file.is_open()) {
    for (uint32_t k=0;k<MaxSuperpixnum;k++) {
      if (sub_map_[k].Initialized()) {
        PixelVector pix;

        Pixels(pix,k);

        for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) {
          if (hpixel_format) {
            if (weighted_map) {
              output_file << iter->HPixnum() << " " <<
                  iter->Superpixnum() << " " <<
                  iter->Resolution() << " " <<
                  iter->Weight() << "\n";
            } else {
              output_file << iter->HPixnum() << " " <<
                  iter->Superpixnum() << " " <<
                  iter->Resolution() << "\n";
            }
          } else {
            if (weighted_map) {
              output_file << iter->Pixnum() << " " <<
                  iter->Resolution() << " " <<
                  iter->Weight() << "\n";
            } else {
              output_file << iter->Pixnum() << " " <<
                  iter->Resolution() << "\n";
            }
          }
        }
      }
    }

    output_file.close();

    return true;
  } else {
    return false;
  }
}

bool Map::Read(std::string& InputFile, bool hpixel_format,
	       bool weighted_map) {
  Clear();

  std::ifstream input_file(InputFile.c_str());

  uint32_t hpixnum, superpixnum, pixnum;
  int resolution;  // have to use an int here because of old map formats
  double weight;
  bool found_file = false;

  if (input_file) {
    found_file = true;
    while (!input_file.eof()) {
      if (hpixel_format) {
	if (weighted_map) {
	  input_file >> hpixnum >> superpixnum >> resolution >> weight;
	} else {
	  input_file >> hpixnum >> superpixnum >> resolution;
	  weight = 1.0;
	}
      } else {
	if (weighted_map) {
	  input_file >> pixnum >> resolution >> weight;
	} else {
	  input_file >> pixnum >> resolution;
	  weight = 1.0;
	}
      }

      if (!input_file.eof() && (resolution % 2 == 0) &&
	  (resolution > 0)) {
	if (!hpixel_format)
	  Pixel::Pix2HPix(static_cast<uint16_t>(resolution), pixnum,
			  hpixnum, superpixnum);
	Pixel tmp_pix(static_cast<uint16_t>(resolution), hpixnum,
		      superpixnum, weight);
	sub_map_[superpixnum].AddPixel(tmp_pix);
      }
    }

    input_file.close();

    bool found_beginning = false;
    for (SubMapIterator iter=sub_map_.begin();iter!=sub_map_.end();++iter) {
      if (iter->Initialized()) {
	if (iter->Unsorted()) iter->Resolve();

	if (!found_beginning) {
	  begin_ = MapIterator(iter->Superpixnum(), iter->Begin());
	  found_beginning = true;
	}
	end_ = MapIterator(iter->Superpixnum(), iter->End());

	area_ += iter->Area();
	size_ += iter->Size();
	if (min_resolution_ > iter->MinResolution())
	  min_resolution_ = iter->MinResolution();
	if (max_resolution_ < iter->MaxResolution())
	  max_resolution_ = iter->MaxResolution();
	if (iter->MinWeight() < min_weight_) min_weight_ = iter->MinWeight();
	if (iter->MaxWeight() > max_weight_) max_weight_ = iter->MaxWeight();
	for (uint16_t resolution_iter=HPixResolution, i=0;
	     i<ResolutionLevels;resolution_iter*=2, i++) {
	  pixel_count_[resolution_iter] += iter->PixelCount(resolution_iter);
	}
      }
    }

    if (!found_beginning) found_file = false;
  } else {
    std::cout << InputFile << " does not exist!.  No Map ingested\n";
  }

  return found_file;
}

double Map::AverageWeight() {
  double unmasked_fraction = 0.0, weighted_average = 0.0;

  for (SubMapIterator iter=sub_map_.begin();iter!=sub_map_.end();++iter) {
    if (iter->Initialized()) {
      weighted_average += iter->AverageWeight()*iter->Area();
      unmasked_fraction += iter->Area();
    }
  }

  if (unmasked_fraction > 0.000000001) weighted_average /= unmasked_fraction;

  return weighted_average;
}

bool Map::IngestMap(PixelVector& pix, bool destroy_copy) {
  for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) {
    uint32_t k = iter->Superpixnum();
    sub_map_[k].AddPixel(*iter);
  }

  if (destroy_copy) pix.clear();

  for (uint32_t k=0;k<MaxSuperpixnum;k++)
    if (sub_map_[k].Unsorted()) sub_map_[k].Resolve();

  return Initialize();
}

bool Map::IngestMap(Map& stomp_map, bool destroy_copy) {
  PixelVector super_pix;

  stomp_map.Coverage(super_pix);

  for (PixelIterator iter=super_pix.begin();iter!=super_pix.end();++iter) {
    PixelVector tmp_pix;

    stomp_map.Pixels(tmp_pix,iter->Superpixnum());

    for (uint32_t i=0;i<tmp_pix.size();i++) {
      sub_map_[iter->Superpixnum()].AddPixel(tmp_pix[i]);
    }
  }

  if (destroy_copy) stomp_map.Clear();

  for (uint32_t k=0;k<MaxSuperpixnum;k++)
    if (sub_map_[k].Unsorted()) sub_map_[k].Resolve();

  return Initialize();
}

bool Map::IntersectMap(PixelVector& pix) {
  Map stomp_map;

  stomp_map.Initialize(pix);

  return IntersectMap(stomp_map);
}

bool Map::IntersectMap(Map& stomp_map) {
  bool found_overlapping_area = false;

  uint32_t superpixnum = 0;

  // First, just check to see that we've got some overlapping area between
  // the two maps.
  while (superpixnum < MaxSuperpixnum && !found_overlapping_area) {
    if (sub_map_[superpixnum].Initialized() &&
	stomp_map.ContainsSuperpixel(superpixnum)) {
      if (Area(superpixnum) < stomp_map.Area(superpixnum)) {
	PixelVector tmp_pix;
        PixelVector match_pix;

        sub_map_[superpixnum].Pixels(tmp_pix);

	// FindMatchingPixels using the weights from this map
	stomp_map.FindMatchingPixels(tmp_pix, match_pix, false);

	if (!match_pix.empty()) found_overlapping_area = true;
      } else {
	PixelVector tmp_pix;
        PixelVector match_pix;

        stomp_map.Pixels(tmp_pix, superpixnum);

	// FindMatchingPixels using the weights from this map
	FindMatchingPixels(tmp_pix, match_pix, true);

	if (!match_pix.empty()) found_overlapping_area = true;
      }
    }
    superpixnum++;
  }

  // Provided that we've got some overlap, now do a full calculation for the
  // whole map.
  if (found_overlapping_area) {
    for (superpixnum=0;superpixnum<MaxSuperpixnum;superpixnum++) {
      if (sub_map_[superpixnum].Initialized() &&
	  stomp_map.ContainsSuperpixel(superpixnum)) {
        if (Area(superpixnum) < stomp_map.Area(superpixnum)) {
          PixelVector tmp_pix;
          PixelVector match_pix;

          sub_map_[superpixnum].Pixels(tmp_pix);

          stomp_map.FindMatchingPixels(tmp_pix, match_pix, false);

          sub_map_[superpixnum].Clear();

          if (!match_pix.empty()) {
            found_overlapping_area = true;

            for (PixelIterator match_iter=match_pix.begin();
                 match_iter!=match_pix.end();++match_iter) {
              sub_map_[superpixnum].AddPixel(*match_iter);
            }
            sub_map_[superpixnum].Resolve();

            match_pix.clear();
          }

	  tmp_pix.clear();
        } else {
          PixelVector tmp_pix;
          PixelVector match_pix;

          stomp_map.Pixels(tmp_pix, superpixnum);

          FindMatchingPixels(tmp_pix, match_pix, true);

          sub_map_[superpixnum].Clear();

          if (!match_pix.empty()) {
            found_overlapping_area = true;

            for (PixelIterator match_iter=match_pix.begin();
                 match_iter!=match_pix.end();++match_iter) {
              sub_map_[superpixnum].AddPixel(*match_iter);
            }
            sub_map_[superpixnum].Resolve();

            match_pix.clear();
          }

	  tmp_pix.clear();
        }
      } else {
	// If there are no pixels in the input map for this superpixel, then
	// clear it out.
        if (sub_map_[superpixnum].Initialized())
	  sub_map_[superpixnum].Clear();
      }
    }
    found_overlapping_area = Initialize();
  }

  return found_overlapping_area;
}

bool Map::AddMap(PixelVector& pix, bool drop_single) {
  Map stomp_map;
  stomp_map.Initialize(pix);

  return AddMap(stomp_map,drop_single);
}

bool Map::AddMap(Map& stomp_map, bool drop_single) {
  for (uint32_t k=0;k<MaxSuperpixnum;k++) {
    if (sub_map_[k].Initialized() && stomp_map.ContainsSuperpixel(k)) {
      // Ok, we've got 2 maps in this superpixel, so we have to break
      // both down and calculate the overlap.
      sub_map_[k].Add(stomp_map, drop_single);
    } else {
      // Ok, only one map covers this superpixel, so we can just copy
      // all of the pixels directly into the final map.  If it's only in
      // our current map, then we don't want to do anything, so we skip that
      // case (unless we're dropping non-overlapping area, in which case we
      // clear that superpixel out).

      if (drop_single) {
        if (sub_map_[k].Initialized()) sub_map_[k].Clear();
      } else {
        if (stomp_map.ContainsSuperpixel(k)) {
          PixelVector added_pix;

          stomp_map.Pixels(added_pix,k);

          for (PixelIterator iter=added_pix.begin();
               iter!=added_pix.end();++iter) sub_map_[k].AddPixel(*iter);

          sub_map_[k].Resolve();
        }
      }
    }
  }

  return Initialize();
}

bool Map::MultiplyMap(PixelVector& pix, bool drop_single) {
  Map stomp_map;
  stomp_map.Initialize(pix);

  return MultiplyMap(stomp_map,drop_single);
}

bool Map::MultiplyMap(Map& stomp_map, bool drop_single) {
  for (uint32_t k=0;k<MaxSuperpixnum;k++) {
    if (sub_map_[k].Initialized() && stomp_map.ContainsSuperpixel(k)) {
      // Ok, we've got 2 maps in this superpixel, so we have to break
      // both down and calculate the overlap.
      sub_map_[k].Multiply(stomp_map, drop_single);
    } else {
      // Ok, only one map covers this superpixel, so we can just copy
      // all of the pixels directly into the final map.  If it's only in
      // our current map, then we don't want to do anything, so we skip that
      // case (unless we're dropping non-overlapping area, in which case we
      // clear that superpixel out).

      if (drop_single) {
        if (sub_map_[k].Initialized()) sub_map_[k].Clear();
      } else {
        if (stomp_map.ContainsSuperpixel(k)) {
          PixelVector multi_pix;

          stomp_map.Pixels(multi_pix,k);

          for (PixelIterator iter=multi_pix.begin();
               iter!=multi_pix.end();++iter) sub_map_[k].AddPixel(*iter);

          sub_map_[k].Resolve();
        }
      }
    }
  }

  return Initialize();
}

bool Map::ExcludeMap(PixelVector& pix, bool destroy_copy) {
  Map stomp_map;
  stomp_map.Initialize(pix);

  if (destroy_copy) pix.clear();

  return ExcludeMap(stomp_map,destroy_copy);
}

bool Map::ExcludeMap(Map& stomp_map, bool destroy_copy) {
  PixelVector super_pix;
  stomp_map.Coverage(super_pix);

  for (PixelIterator iter=super_pix.begin();iter!=super_pix.end();++iter) {
    uint32_t superpixnum = iter->Superpixnum();
    if (sub_map_[superpixnum].Initialized()) {
      sub_map_[superpixnum].Exclude(stomp_map);
    }
  }

  if (destroy_copy) stomp_map.Clear();

  return Initialize();
}

bool Map::ImprintMap(PixelVector& pix) {
  Map stomp_map;
  stomp_map.Initialize(pix);

  return ImprintMap(stomp_map);
}

bool Map::ImprintMap(Map& stomp_map) {
  bool found_overlapping_area = false;

  uint32_t k = 0;
  while (k<MaxSuperpixnum && !found_overlapping_area) {
    if (sub_map_[k].Initialized() && stomp_map.ContainsSuperpixel(k)) {
      if (Area(k) < stomp_map.Area(k)) {
	PixelVector tmp_pix;
        PixelVector match_pix;

        sub_map_[k].Pixels(tmp_pix);

	stomp_map.FindMatchingPixels(tmp_pix,match_pix,false);

	if (!match_pix.empty()) found_overlapping_area = true;
      } else {
	PixelVector tmp_pix;
        PixelVector match_pix;

        stomp_map.Pixels(tmp_pix,k);

	FindMatchingPixels(tmp_pix, match_pix, true);

	if (!match_pix.empty()) found_overlapping_area = true;
      }
    }
    k++;
  }

  if (found_overlapping_area) {
    for (k=0;k<MaxSuperpixnum;k++) {
      if (sub_map_[k].Initialized() && stomp_map.ContainsSuperpixel(k)) {
        if (Area(k) < stomp_map.Area(k)) {
          PixelVector tmp_pix;
          PixelVector match_pix;

          sub_map_[k].Pixels(tmp_pix);

          stomp_map.FindMatchingPixels(tmp_pix,match_pix,true);

          sub_map_[k].Clear();

          if (!match_pix.empty()) {
            found_overlapping_area = true;

            for (PixelIterator match_iter=match_pix.begin();
                 match_iter!=match_pix.end();++match_iter) {
              sub_map_[k].AddPixel(*match_iter);
            }
            sub_map_[k].Resolve();

            match_pix.clear();
          }
        } else {
          PixelVector tmp_pix;
          PixelVector match_pix;

          stomp_map.Pixels(tmp_pix,k);

          FindMatchingPixels(tmp_pix,match_pix,false);

          sub_map_[k].Clear();

          if (!match_pix.empty()) {
            found_overlapping_area = true;

            for (PixelIterator match_iter=match_pix.begin();
                 match_iter!=match_pix.end();++match_iter) {
              sub_map_[k].AddPixel(*match_iter);
            }
            sub_map_[k].Resolve();

            match_pix.clear();
          }
        }
      } else {
        if (sub_map_[k].Initialized()) sub_map_[k].Clear();
      }
    }
    found_overlapping_area = Initialize();
  }

  return found_overlapping_area;
}

void Map::ScaleWeight(const double weight_scale) {
  for (uint32_t k=0;k<MaxSuperpixnum;k++)
    if (sub_map_[k].Initialized()) sub_map_[k].ScaleWeight(weight_scale);
}

void Map::AddConstantWeight(const double add_weight) {
  for (uint32_t k=0;k<MaxSuperpixnum;k++)
    if (sub_map_[k].Initialized()) sub_map_[k].AddConstantWeight(add_weight);
}

void Map::InvertWeight() {
  for (uint32_t k=0;k<MaxSuperpixnum;k++)
    if (sub_map_[k].Initialized()) sub_map_[k].InvertWeight();
}

void Map::Pixels(PixelVector& pix, uint32_t superpixnum) {
  if (!pix.empty()) pix.clear();

  if (superpixnum < MaxSuperpixnum) {
    sub_map_[superpixnum].Pixels(pix);
  } else {
    for (uint32_t k=0;k<MaxSuperpixnum;k++) {
      if (sub_map_[k].Initialized()) {
        PixelVector tmp_pix;

        sub_map_[k].Pixels(tmp_pix);

        for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter)
          pix.push_back(*iter);
      }
    }
  }
}

MapIterator Map::Begin() {
  return begin_;
}

MapIterator Map::End() {
  return end_;
}

void Map::Iterate(MapIterator* iter) {
  ++iter->second;
  if (iter->second == sub_map_[iter->first].End() &&
      iter->first != end_.first) {
    bool found_next_iterator = false;
    while (iter->first < MaxSuperpixnum && !found_next_iterator) {
      iter->first++;
      if (sub_map_[iter->first].Initialized()) {
	iter->second = sub_map_[iter->first].Begin();
	found_next_iterator = true;
      }
    }
  }
}

void Map::Clear() {
  area_ = 0.0;
  size_ = 0;
  min_resolution_ = MaxPixelResolution;
  max_resolution_ = HPixResolution;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();

  for (uint16_t resolution=HPixResolution, i=0;
       i<ResolutionLevels;resolution*=2, i++)
    pixel_count_[resolution] = 0;

  if (!sub_map_.empty()) {
    for (SubMapIterator iter=sub_map_.begin();iter!=sub_map_.end();++iter)
      iter->Clear();
    sub_map_.clear();
  }

  sub_map_.reserve(MaxSuperpixnum);

  for (uint32_t k=0;k<MaxSuperpixnum;k++)
    sub_map_.push_back(SubMap(k));

  begin_ = end_;
}

void Map::Clear(uint32_t superpixnum) {
  if (superpixnum < MaxSuperpixnum)
    sub_map_[superpixnum].Clear();
}

bool Map::ContainsSuperpixel(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].Initialized() : false);
}

double Map::Area() {
  return area_;
}

double Map::Area(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].Area() : 0.0);
}

uint16_t Map::MinResolution() {
  return min_resolution_;
}

uint16_t Map::MinResolution(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].MinResolution() : 0);
}

uint16_t Map::MaxResolution() {
  return max_resolution_;
}

uint16_t Map::MaxResolution(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].MaxResolution() : 0);
}

double Map::MinWeight() {
  return min_weight_;
}

double Map::MinWeight(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].MinWeight() : 0.0);
}

double Map::MaxWeight() {
  return max_weight_;
}

double Map::MaxWeight(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].MaxWeight() : 0.0);
}

uint32_t Map::Size() {
  return size_;
}

uint32_t Map::Size(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].Size() : 0);
}

bool Map::Empty() {
  return (size_ == 0 ? true : false);
}

uint32_t Map::PixelCount(uint16_t resolution) {
  return (!(resolution % 2) ? pixel_count_[resolution] : 0);
}


} // end namespace Stomp
