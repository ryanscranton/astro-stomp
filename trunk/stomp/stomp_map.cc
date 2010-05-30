// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)
// vim: set et ts=2 sw=2:
//
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
#include "stomp_geometry.h"

namespace Stomp {

int Map::_INSIDE_MAP=1;
int Map::_FIRST_QUADRANT_OK=2;
int Map::_SECOND_QUADRANT_OK=4;
int Map::_THIRD_QUADRANT_OK=8;
int Map::_FOURTH_QUADRANT_OK=16;


SubMap::SubMap(uint32_t superpixnum) {
  superpixnum_ = superpixnum;
  area_ = 0.0;
  size_ = 0;
  min_level_ = MaxPixelLevel;
  max_level_ = HPixLevel;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  Pixel::PixelBound(HPixResolution, superpixnum, lambda_min_,
		    lambda_max_, eta_min_, eta_max_);
  z_min_ = sin(lambda_min_*DegToRad);
  z_max_ = sin(lambda_max_*DegToRad);
  initialized_ = false;
  unsorted_ = false;

  for (uint32_t resolution=HPixResolution;
       resolution<=MaxPixelResolution;resolution*=2) {
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
    if (pix.Level() < min_level_) min_level_ = pix.Level();
    if (pix.Level() > max_level_) max_level_ = pix.Level();
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
    min_level_ = MaxPixelLevel;
    max_level_ = HPixLevel;
    min_weight_ = 1.0e30;
    max_weight_ = -1.0e30;
    for (uint32_t resolution=HPixResolution;
	 resolution<=MaxPixelResolution;resolution*=2) {
      pixel_count_[resolution] = 0;
    }

    for (PixelIterator iter=pix_.begin();iter!=pix_.end();++iter) {
      area_ += iter->Area();
      if (iter->Level() < min_level_) min_level_ = iter->Level();
      if (iter->Level() > max_level_) max_level_ = iter->Level();
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

void SubMap::SetMaximumResolution(uint32_t max_resolution,
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

  for (uint32_t resolution=MinResolution();
       resolution<=MaxResolution();resolution*=2) {
    Pixel tmp_pix(ang, resolution);
    PixelPair iter = equal_range(pix_.begin(), pix_.end(), tmp_pix,
                                 Pixel::SuperPixelBasedOrder);
    if (iter.first != iter.second) {
      keep = true;
      weight = iter.first->Weight();
    }
    if (keep) break;
  }

  return keep;
}

double SubMap::FindUnmaskedFraction(Pixel& pix) {
  PixelIterator iter;
  if (pix.Level() >= MaxLevel()) {
    iter = pix_.end();
  } else {
    Pixel tmp_pix(pix.PixelX0()*2, pix.PixelY0()*2,
		  pix.Resolution()*2, 1.0);
    iter = lower_bound(pix_.begin(),pix_.end(),
                       tmp_pix,Pixel::SuperPixelBasedOrder);
  }

  uint8_t level = MinLevel();
  double unmasked_fraction = 0.0;
  bool found_pixel = false;
  while (level <= pix.Level() && !found_pixel) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToLevel(level);
    PixelPair super_iter = equal_range(pix_.begin(), iter, tmp_pix,
                                       Pixel::SuperPixelBasedOrder);
    if (super_iter.first != super_iter.second) {
      found_pixel = true;
      unmasked_fraction = 1.0;
    }
    level++;
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
  if (pix.Level() >= MaxLevel()) {
    iter = pix_.end();
  } else {
    Pixel tmp_pix(pix.PixelX0()*2, pix.PixelY0()*2,
		  pix.Resolution()*2, 1.0);
    iter = lower_bound(pix_.begin(), pix_.end(), tmp_pix,
                       Pixel::SuperPixelBasedOrder);
  }

  uint8_t level = MinLevel();
  int8_t unmasked_status = 0;
  while ((level <= pix.Level()) && (unmasked_status == 0)) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToLevel(level);
    PixelPair super_iter = equal_range(pix_.begin(),iter,tmp_pix,
                                       Pixel::SuperPixelBasedOrder);
    if (super_iter.first != super_iter.second) unmasked_status = 1;
    level++;
  }

  while ((iter != pix_.end()) && (unmasked_status == 0)) {
    if (pix.Contains(*iter)) unmasked_status = -1;
    ++iter;
  }

  return unmasked_status;
}

double SubMap::FindAverageWeight(Pixel& pix) {
  PixelIterator iter;

  if (pix.Level() >= MaxLevel()) {
    iter = pix_.end();
  } else {
    Pixel tmp_pix(pix.PixelX0()*2, pix.PixelY0()*2,
		  pix.Resolution()*2, 1.0);
    iter = lower_bound(pix_.begin(),pix_.end(),tmp_pix,
                       Pixel::SuperPixelBasedOrder);
  }

  double unmasked_fraction = 0.0, weighted_average = 0.0;
  bool found_pixel = false;
  uint8_t level = MinLevel();
  while (level <= pix.Level() && !found_pixel) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToLevel(level);
    PixelPair super_iter = equal_range(pix_.begin(), iter, tmp_pix,
                                       Pixel::SuperPixelBasedOrder);
    if (super_iter.first != super_iter.second) {
      found_pixel = true;
      weighted_average = super_iter.first->Weight();
      unmasked_fraction = 1.0;
    }
    level++;
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
  if (pix.Level() >= MaxLevel()) {
    iter = pix_.end();
  } else {
    Pixel tmp_pix(pix.PixelX0()*2, pix.PixelY0()*2, pix.Resolution()*2, 1.0);
    iter = lower_bound(pix_.begin(),pix_.end(),tmp_pix,
                       Pixel::SuperPixelBasedOrder);
  }

  uint8_t level = MinLevel();
  while (level <= pix.Level() && !found_pixel) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToLevel(level);
    find_iter = lower_bound(pix_.begin(), iter, tmp_pix,
                            Pixel::SuperPixelBasedOrder);
    if (Pixel::PixelMatch(*find_iter,tmp_pix)) {
      found_pixel = true;
      tmp_pix = pix;
      if (use_local_weights) tmp_pix.SetWeight(find_iter->Weight());
      match_pix.push_back(tmp_pix);
    }
    level++;
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

void SubMap::Soften(PixelVector& output_pix, uint32_t max_resolution,
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
      stomp_map.FindMatchingPixels(*iter, match_pix, true);
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
	stomp_map.FindMatchingPixels(*iter, match_pix, true);
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
  // area that is contained in just one map, we need to find the pixels in the
  // input Map that weren't contained entirely in our current Map.  The ones
  // that are entirely outside our Map can be added directly.  The ones that are
  // partially inside our Map need to be refined until we find the pieces that
  // are outside our Map (the parts that were inside have already been handled
  // by the procedure above).
  if (!drop_single) {
    resolve_pix.clear();
    PixelVector stomp_pix;
    stomp_map.Pixels(stomp_pix, Superpixnum());

    for (PixelIterator iter=stomp_pix.begin();iter!=stomp_pix.end();++iter) {
      int8_t status = FindUnmaskedStatus(*iter);

      // If the pixel is completely outside our Map, we keep it.
      if (status == 0) keep_pix.push_back(*iter);

      // If it's partially inside, then we need to refine it to find the part
      // not contained in our Map.
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
    stomp_pix.clear();

    // Now we iterate over the partial pixels until we've cleared them out.
    while (resolve_pix.size() > 0) {
      PixelVector tmp_pix;
      tmp_pix.reserve(resolve_pix.size());
      for (PixelIterator iter=resolve_pix.begin();
	   iter!=resolve_pix.end();++iter) tmp_pix.push_back(*iter);

      resolve_pix.clear();

      for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter) {
	int8_t status = FindUnmaskedStatus(*iter);

	if (status == 0) keep_pix.push_back(*iter);

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
      stomp_map.FindMatchingPixels(*iter, match_pix, true);
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
	stomp_map.FindMatchingPixels(*iter, match_pix, true);
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
  // area that is contained in just one map, we need to find the pixels in the
  // input Map that weren't contained entirely in our current Map.  The ones
  // that are entirely outside our Map can be added directly.  The ones that are
  // partially inside our Map need to be refined until we find the pieces that
  // are outside our Map (the parts that were inside have already been handled
  // by the procedure above).
  if (!drop_single) {
    resolve_pix.clear();
    PixelVector stomp_pix;
    stomp_map.Pixels(stomp_pix, Superpixnum());

    for (PixelIterator iter=stomp_pix.begin();iter!=stomp_pix.end();++iter) {
      int8_t status = FindUnmaskedStatus(*iter);

      // If the pixel is completely outside our Map, we keep it.
      if (status == 0) keep_pix.push_back(*iter);

      // If it's partially inside, then we need to refine it to find the part
      // not contained in our Map.
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
    stomp_pix.clear();

    // Now we iterate over the partial pixels until we've cleared them out.
    while (resolve_pix.size() > 0) {
      PixelVector tmp_pix;
      tmp_pix.reserve(resolve_pix.size());
      for (PixelIterator iter=resolve_pix.begin();
	   iter!=resolve_pix.end();++iter) tmp_pix.push_back(*iter);

      resolve_pix.clear();

      for (PixelIterator iter=tmp_pix.begin();iter!=tmp_pix.end();++iter) {
	int8_t status = FindUnmaskedStatus(*iter);

	if (status == 0) keep_pix.push_back(*iter);

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
  min_level_ = MaxPixelLevel;
  max_level_ = HPixLevel;
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

void SubMap::SetUnsorted() {
  unsorted_ = true;
}

uint32_t SubMap::MinResolution() {
  return Pixel::LevelToResolution(min_level_);
}

uint32_t SubMap::MaxResolution() {
  return Pixel::LevelToResolution(max_level_);
}

uint8_t SubMap::MinLevel() {
  return min_level_;
}

uint8_t SubMap::MaxLevel() {
  return max_level_;
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

uint32_t SubMap::PixelCount(uint32_t resolution) {
  return (!(resolution % 2) ? pixel_count_[resolution] : 0);
}

Map::Map() {
  area_ = 0.0;
  size_ = 0;
  min_level_ = MaxPixelLevel;
  max_level_ = HPixLevel;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();

  for (uint32_t resolution=HPixResolution;
       resolution<=MaxPixelResolution;resolution*=2)
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
  min_level_ = MaxPixelLevel;
  max_level_ = HPixLevel;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();
  for (uint32_t resolution=HPixResolution;
       resolution<=MaxPixelResolution;resolution*=2) {
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
      if (min_level_ > iter->MinLevel()) min_level_ = iter->MinLevel();
      if (max_level_ < iter->MaxLevel()) max_level_ = iter->MaxLevel();
      if (iter->MinWeight() < min_weight_) min_weight_ = iter->MinWeight();
      if (iter->MaxWeight() > max_weight_) max_weight_ = iter->MaxWeight();
      for (uint32_t resolution=HPixResolution;
	   resolution<=MaxPixelResolution;resolution*=2) {
	pixel_count_[resolution] += iter->PixelCount(resolution);
      }
    }
  }

}

Map::Map(const std::string& InputFile, bool hpixel_format, bool weighted_map) {
  Read(InputFile, hpixel_format, weighted_map);
}

Map::Map(GeometricBound& bound, double weight, uint32_t max_resolution,
	 bool verbose) {
  if (PixelizeBound(bound, weight, max_resolution) && verbose) {
    std::cout << "Successfully pixelized GeometricBound.\n" <<
      "\tOriginal Area: " << bound.Area() << " sq. degrees; " <<
      "Pixelized Area: " << Area() << " sq. degrees.\n";
  } else {
    if (verbose) std::cout << "Pixelization failed.\n";
  }
}


Map::~Map() {
  min_level_ = max_level_ = 0;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  Clear();
}

bool Map::Initialize() {
  area_ = 0.0;
  size_ = 0;
  min_level_ = MaxPixelLevel;
  max_level_ = HPixLevel;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();
  for (uint32_t resolution=HPixResolution;
       resolution<=MaxPixelResolution;resolution*=2) {
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
      if (min_level_ > iter->MinLevel()) min_level_ = iter->MinLevel();
      if (max_level_ < iter->MaxLevel()) max_level_ = iter->MaxLevel();
      if (iter->MinWeight() < min_weight_) min_weight_ = iter->MinWeight();
      if (iter->MaxWeight() > max_weight_) max_weight_ = iter->MaxWeight();
      for (uint32_t resolution=HPixResolution;
	   resolution<=MaxPixelResolution;resolution*=2) {
	pixel_count_[resolution] += iter->PixelCount(resolution);
      }
    }
  }

  return found_valid_superpixel;
}

bool Map::Initialize(PixelVector& pix, bool force_resolve) {
  area_ = 0.0;
  size_ = 0;
  min_level_ = MaxPixelLevel;
  max_level_ = HPixLevel;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();
  for (uint32_t resolution=HPixResolution;
       resolution<=MaxPixelResolution;resolution*=2) {
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
      if (min_level_ > iter->MinLevel()) min_level_ = iter->MinLevel();
      if (max_level_ < iter->MaxLevel()) max_level_ = iter->MaxLevel();
      if (iter->MinWeight() < min_weight_) min_weight_ = iter->MinWeight();
      if (iter->MaxWeight() > max_weight_) max_weight_ = iter->MaxWeight();
      for (uint32_t resolution=HPixResolution;
	   resolution<=MaxPixelResolution;resolution*=2) {
	pixel_count_[resolution] += iter->PixelCount(resolution);
      }
    }
  }

  return found_valid_superpixel;
}

void Map::AddPixel(Pixel& pix) {
  sub_map_[pix.Superpixnum()].AddPixel(pix);

  area_ += pix.Area();
  size_++;

  if (min_level_ > pix.Level()) min_level_ = pix.Level();
  if (max_level_ < pix.Level()) max_level_ = pix.Level();
  if (pix.Weight() < min_weight_) min_weight_ = pix.Weight();
  if (pix.Weight() > max_weight_) max_weight_ = pix.Weight();
  pixel_count_[pix.Resolution()]++;
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

bool Map::Contains(AngularCoordinate& ang) {
  bool keep = false;
  double weight;

  uint32_t k;
  Pixel::Ang2Pix(HPixResolution, ang, k);

  if (sub_map_[k].Initialized()) keep = sub_map_[k].FindLocation(ang, weight);

  return keep;
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

void Map::Coverage(PixelVector& superpix, uint32_t resolution,
		   bool calculate_fraction) {
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
	      if (calculate_fraction) {
		iter->SetWeight(FindUnmaskedFraction(*iter));
	      } else {
		iter->SetWeight(1.0);
	      }
	    }
	    superpix.push_back(*iter);
	  }
	}
      }
    }
    sort(superpix.begin(), superpix.end(), Pixel::LocalOrder);
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
      uint32_t maximum_resolution = HPixResolution;
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

void Map::Soften(Map& stomp_map, uint32_t maximum_resolution,
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

void Map::Soften(uint32_t maximum_resolution, bool average_weights) {
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

#ifdef WITH_NUMPY  // begin python-only code
// This version returns numerical python arrays in a tuple
// note use_weighted_sampling is optional, default false

// This is the generic version taking a string for the system
PyObject* Map::GenerateRandomPoints(uint32_t n_point, const std::string& system,
				    bool use_weighted_sampling)
  throw (const char*)  {
  Stomp::AngularCoordinate::Sphere sys =
    Stomp::AngularCoordinate::SystemFromString(system);
  return GenerateRandomPoints(n_point,sys,use_weighted_sampling);
}

// This is the generic version taking a Sphere id for the system
PyObject* Map::GenerateRandomPoints(uint32_t n_point,
				    Stomp::AngularCoordinate::Sphere systemid,
				    bool use_weighted_sampling)
  throw (const char*)  {
  double minimum_probability = -0.0001;
  double probability_slope = 0.0;
  std::stringstream err;
  // Make the output numpy arrays
  NumpyVector<double> x1(n_point);
  NumpyVector<double> x2(n_point);

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

      z = sub_map_[k].ZMin() +
	mtrand.rand(sub_map_[k].ZMax() - sub_map_[k].ZMin());
      lambda = asin(z)*RadToDeg;
      eta = sub_map_[k].EtaMin() +
	mtrand.rand(sub_map_[k].EtaMax() - sub_map_[k].EtaMin());
      tmp_ang.SetSurveyCoordinates(lambda,eta);

      keep = sub_map_[k].FindLocation(tmp_ang,weight);

      if (use_weighted_sampling && keep) {
	probability_limit =
	  minimum_probability + (weight - min_weight_)*probability_slope;
	if (mtrand.rand(1.0) > probability_limit) keep = false;
      }
    }

    switch (systemid) {
    case Stomp::AngularCoordinate::Survey:
      x1[m] = lambda;
      x2[m] = eta;
      break;
    case Stomp::AngularCoordinate::Equatorial:
      x1[m] = tmp_ang.RA();
      x2[m] = tmp_ang.DEC();
      break;
    case Stomp::AngularCoordinate::Galactic:
      x1[m] = tmp_ang.GalLon();
      x2[m] = tmp_ang.GalLat();
      break;
    default:
      err << "Bad system id: " << systemid;
      throw err.str().c_str();
    }
  }

  PyObject* output_tuple = PyTuple_New(2);
  PyTuple_SetItem(output_tuple, 0, x1.getref());
  PyTuple_SetItem(output_tuple, 1, x2.getref());

  return output_tuple;
}

// These are wrappers for the more generic function above
PyObject* Map::GenerateRandomEq(uint32_t n_point, bool use_weighted_sampling)
  throw (const char*)  {
  return GenerateRandomPoints(n_point, Stomp::AngularCoordinate::Equatorial,
			      use_weighted_sampling);
}
PyObject* Map::GenerateRandomSurvey(uint32_t n_point,
				    bool use_weighted_sampling)
  throw (const char*)  {
  return GenerateRandomPoints(n_point, Stomp::AngularCoordinate::Survey,
			      use_weighted_sampling);
}
PyObject* Map::GenerateRandomGal(uint32_t n_point, bool use_weighted_sampling)
  throw (const char*)  {
  return GenerateRandomPoints(n_point, Stomp::AngularCoordinate::Galactic,
			      use_weighted_sampling);
}

PyObject* Map::Contains(PyObject* x1obj, PyObject* x2obj,
			const std::string& system, PyObject* radobj)
  throw (const char* ) {
  // convert the string system indicator to a Sphere id
  Stomp::AngularCoordinate::Sphere sys =
    Stomp::AngularCoordinate::SystemFromString(system);

  // Get numpy arrays objects.   No copy made as long as the type and byte
  // order is correct.
  NumpyVector<double> x1(x1obj);
  NumpyVector<double> x2(x2obj);
  if (x1.size() != x2.size()) {
    throw "coordinates must be same size";
  }
  NumpyVector<double> rad;
  npy_intp nrad=0;
  double thisrad=-1;
  if (radobj != NULL) {
    rad.init(radobj);
    nrad = rad.size();
    if (nrad != 1 && nrad != x1.size()) {
      throw "radius must be same length as coordinates or length 1";
    }
    thisrad = rad[0];
    // seed the random number generator.  Distinctive to
    // one second
    std::srand ( std::time(NULL) );
  }

  NumpyVector<npy_int8> maskflags(x1.size());

  Stomp::AngularCoordinate ang;
  for (npy_intp i=0; i<x1.size(); i++) {
    ang.Set(x1[i], x2[i], sys);
    if (Contains(ang)) {
      maskflags[i] |= _INSIDE_MAP;

      // If radii were sent, we will do the quadrant check
      if (nrad > 0) {
        if (nrad > 1) {
          thisrad = rad[i];
        }
        maskflags[i] |= QuadrantsContainedMC(ang,thisrad,sys);
      }
    }
  }
  return maskflags.getref();
}

#endif  // end python-only code

int Map::QuadrantsContainedMC(AngularCoordinate& ang, double radius,
			      Stomp::AngularCoordinate::Sphere coord_system)
  throw (const char*) {
  int maskflags=0;

  if (1) {
    if (QuadrantContainedMC(ang,radius,0)) {
      maskflags |= _FIRST_QUADRANT_OK;
    }
    if (QuadrantContainedMC(ang,radius,1)) {
      maskflags |= _SECOND_QUADRANT_OK;
    }
    if (QuadrantContainedMC(ang,radius,2)) {
      maskflags |= _THIRD_QUADRANT_OK;
    }
    if (QuadrantContainedMC(ang,radius,3)) {
      maskflags |= _FOURTH_QUADRANT_OK;
    }
  } else {
    double amin=0.05;
    double pmax=0.01;
    {
      WedgeBound wb(ang, radius, 0.0, 90.0, coord_system);
      if (Contains(wb, amin, pmax)) 
        maskflags |= _FIRST_QUADRANT_OK;
    }
    {
      WedgeBound wb(ang, radius, 90.0, 180.0, coord_system);
      if (Contains(wb, amin, pmax))
        maskflags |= _SECOND_QUADRANT_OK;
    }
    {
      WedgeBound wb(ang, radius, 180.0, 270.0, coord_system);
      if (Contains(wb, amin, pmax))
        maskflags |= _THIRD_QUADRANT_OK;
    }
    {
      WedgeBound wb(ang, radius, 270.0, 360.0, coord_system);
      if (Contains(wb, amin, pmax))
        maskflags |= _FOURTH_QUADRANT_OK;
    }

  }
  return maskflags;
}


bool Map::QuadrantContainedMC(AngularCoordinate& ang, double radius,
			      int quadrant) throw (const char*) {
  // These tune the test accuracy
  //
  // Minimum size we want to resolve in square degrees This can be pretty
  // big since we are only worried about edges and big holes.  0.05
  // corresponds to half a field
  static double amin=0.05;

  // probability a randomly generated point will *not*
  // fall within our smallest area amin
  double A = Pi*radius*radius/4.0;
  double pmiss = 1. - amin/A;

  // Probability of missing amin sized region?  The number of points used
  // goes as the log of this
  static double pmax=0.01;

  uint32_t nrand=0;
  if (pmiss > 1.e-10) {
    // how many points do we need in order for the
    // probability of missing the hole to be pmax?
    //      We need n such that (1-amin/a)^n = pmax
    double tmp = log10(pmax)/log10(pmiss);
    if (tmp < 20) tmp = 20;
    if (tmp > 20000) tmp = 20000;
    nrand = lround(tmp);
  } else {
    // we reach here often because the search area is very
    // close to or smaller than our minimum resolvable area
    // We don't want nrand to be less than say 20
    nrand = 20;
  }

  // get test the point as clambda,ceta
  double clambda = ang.Lambda();
  double ceta = ang.Eta();

  double rand_clambda=0;
  double rand_ceta=0;

  // use clambda,ceta coords
  Stomp::AngularCoordinate::Sphere sys= Stomp::AngularCoordinate::Survey;
  AngularCoordinate tmp_ang;
  bool quadrant_inside = true;
  for (uint32_t i=0; i<nrand; i++) {
    _GenerateRandLamEtaQuadrant(clambda, ceta, radius, quadrant,
				rand_clambda, rand_ceta);
    tmp_ang.Set(rand_clambda,rand_ceta,sys);
    if (!Contains(tmp_ang)) {
      quadrant_inside = false;
      break;
    }
  }

  return quadrant_inside;
}

void Map::_GenerateRandLamEtaQuadrant(double lambda, double eta,
				      double R, int quadrant,
				      double& rand_lambda, double& rand_eta)
  throw (const char*) {
  // this is crazy slow
  //MTRand mtrand;
  //mtrand.seed();

  std::stringstream err;
  double cospsi;
  double theta,phi,sintheta,costheta,sinphi,cosphi;
  double theta2,costheta2,sintheta2;
  double phi2;
  double Dphi,cosDphi;
  double sinr,cosr;
  double min_theta;

  double rand_r, rand_psi;

  // generate uniformly in R^2
  // random [0,1)
  //rand_r = mtrand.randExc();
  rand_r =
    ( (double)std::rand() / ((double)(RAND_MAX)+(double)(1)) );
  rand_r = sqrt(rand_r)*R*DegToRad;

  // generate theta uniformly from [min_theta,min_theta+90)
  min_theta = quadrant*M_PI/2.;

  rand_psi =
    ( (double)std::rand() / ((double)(RAND_MAX)+(double)(1)) );
  //rand_psi = mtrand.randExc();
  rand_psi = M_PI/2.*rand_psi + min_theta;

  cospsi = cos(rand_psi);


  // [0,180]
  theta = (90.0 - lambda)*DegToRad;
  // [0,360]
  phi   = (eta + 180.0)*DegToRad;

  sintheta = sin(theta);
  costheta = cos(theta);
  sinphi = sin(phi);
  cosphi = cos(phi);

  sinr = sin(rand_r);
  cosr = cos(rand_r);

  costheta2 = costheta*cosr + sintheta*sinr*cospsi;
  if (costheta2 < -1.) costheta2 = -1.;
  if (costheta2 >  1.) costheta2 =  1.;
  theta2 = acos(costheta2);
  sintheta2 = sin(theta2);

  cosDphi = (cosr - costheta*costheta2)/(sintheta*sintheta2);

  if (cosDphi < -1.) cosDphi = -1.;
  if (cosDphi >  1.) cosDphi =  1.;
  Dphi = acos(cosDphi);

  switch(quadrant)
  {
    case 0:
      phi2 = phi + Dphi;
      break;
    case 1:
      phi2 = phi + Dphi;
      break;
    case 2:
      phi2 = phi - Dphi;
      break;
    case 3:
      phi2 = phi - Dphi;
      break;
    default:
      err << "Error: quadrant is undefined: " << quadrant;
      throw err.str().c_str();
  }

  rand_lambda = 90.0 - RadToDeg*theta2;
  rand_eta    = RadToDeg*phi2 - 180.0;
}

//
// amin: Minimum size we want to resolve in square degrees 0.05 corresponds to
// half an SDSS field
//
// pmax: probability a randomly generated point will *not* fall within our
// smallest area amin
bool Map::Contains(GeometricBound& bound, double area_resolution,
		   double precision) {
  double miss_prob = 1.0 - area_resolution/bound.Area();

  uint32_t n_rand = 0;
  if (miss_prob > 1.e-10) {
    // how many points do we need in order for the
    // probability of missing the hole to be pmax?
    //      We need n such that (1-amin/a)^n = pmax
    double points_estimate = log10(precision)/log10(miss_prob);
    if (points_estimate < 20.0) points_estimate = 20.1;
    if (points_estimate > 20000.0)
      points_estimate = 20000.1;  // cap at 20000 random points
    n_rand = static_cast<uint32_t>(points_estimate);
  } else {
    // we reach here often because the search area is very
    // close to or smaller than our minimum resolvable area
    // We don't want nrand to be less than say 20
    n_rand = 20;
  }

  AngularCoordinate tmp_ang;
  bool bound_inside = true;
  for (uint32_t i=0;i<n_rand;i++) {
    bound.GenerateRandomPoint(tmp_ang);
    if (!Contains(tmp_ang)) {
      bound_inside = false;
      break;
    }
  }

  return bound_inside;
}

double Map::FindUnmaskedFraction(GeometricBound& bound, double area_resolution,
				 double precision) {
  double miss_prob = 1.0 - area_resolution/bound.Area();

  uint32_t n_rand = 0;
  if (miss_prob > 1.e-10) {
    // how many points do we need in order for the
    // probability of missing the hole to be pmax?
    //      We need n such that (1-amin/a)^n = pmax
    double points_estimate = log10(precision)/log10(miss_prob);
    if (points_estimate < 20.0) points_estimate = 20.1;
    if (points_estimate > 20000.0)
      points_estimate = 20000.1;  // cap at 20000 random points
    n_rand = static_cast<uint32_t>(points_estimate);
  } else {
    // we reach here often because the search area is very
    // close to or smaller than our minimum resolvable area
    // We don't want nrand to be less than say 20
    n_rand = 20;
  }

  AngularCoordinate tmp_ang;
  uint32_t n_inside = 0;
  for (uint32_t i=0;i<n_rand;i++) {
    bound.GenerateRandomPoint(tmp_ang);
    if (Contains(tmp_ang)) n_inside++;
  }

  return static_cast<double>(n_inside)/static_cast<double>(n_rand);
}

void Map::GenerateRandomPoints(WAngularVector& ang, WAngularVector& input_ang,
			       bool filter_input_points) {
  if (!ang.empty()) ang.clear();
  ang.reserve(input_ang.size());

  PixelVector superpix;
  Coverage(superpix);

  MTRand mtrand;
  mtrand.seed();

  WeightedAngularCoordinate tmp_ang;
  for (uint32_t m=0;m<input_ang.size();m++) {
    if (!filter_input_points ||
	(filter_input_points && Contains(input_ang[m]))) {
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




bool Map::Write(const std::string& OutputFile, bool hpixel_format,
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

bool Map::Read(const std::string& InputFile, bool hpixel_format,
	       bool weighted_map) {
  Clear();

  std::ifstream input_file(InputFile.c_str());

  uint32_t hpixnum, superpixnum, pixnum, x, y;
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
	  Pixel::Pix2HPix(static_cast<uint32_t>(resolution), pixnum,
			  hpixnum, superpixnum);
	Pixel::HPix2XY(static_cast<uint32_t>(resolution), hpixnum, superpixnum,
		       x, y);
	Pixel tmp_pix(x, y, static_cast<uint32_t>(resolution), weight);
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
	if (min_level_ > iter->MinLevel()) min_level_ = iter->MinLevel();
	if (max_level_ < iter->MaxLevel()) max_level_ = iter->MaxLevel();
	if (iter->MinWeight() < min_weight_) min_weight_ = iter->MinWeight();
	if (iter->MaxWeight() > max_weight_) max_weight_ = iter->MaxWeight();
	for (uint32_t resolution_iter=HPixResolution;
	     resolution_iter<=MaxPixelResolution;resolution_iter*=2) {
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

bool Map::PixelizeBound(GeometricBound& bound, double weight,
			uint32_t maximum_resolution) {
  bool pixelized_map = false;

  Clear();

  uint8_t max_resolution_level = MostSignificantBit(maximum_resolution);

  uint8_t starting_resolution_level =
    _FindStartingResolutionLevel(bound.Area());

  if (starting_resolution_level > max_resolution_level)
    starting_resolution_level = max_resolution_level;

  if (starting_resolution_level < HPixLevel)
    starting_resolution_level = HPixLevel;

  // We need to be careful around the poles since the pixels there get
  // very distorted.
  if ((bound.LambdaMin() > 85.0) || (bound.LambdaMax() < -85.0)) {
    starting_resolution_level = MostSignificantBit(512);
    if (max_resolution_level < starting_resolution_level)
      max_resolution_level = starting_resolution_level;
  }

  // std::cout << "Pixelizing from level " <<
  //  static_cast<int>(starting_resolution_level) << " to " <<
  //  static_cast<int>(max_resolution_level) << "...\n";

  uint32_t x_min, x_max, y_min, y_max;
  if (_FindXYBounds(starting_resolution_level, bound,
		    x_min, x_max, y_min, y_max)) {

    PixelVector resolve_pix, previous_pix, kept_pix;
    double pixel_area = 0.0;

    for (uint8_t resolution_level=starting_resolution_level;
         resolution_level<=max_resolution_level;resolution_level++) {
      uint32_t resolution = Pixel::LevelToResolution(resolution_level);

      unsigned n_keep = 0;
      uint32_t nx = Nx0*resolution;
      Pixel tmp_pix;
      tmp_pix.SetResolution(resolution);
      double unit_area = tmp_pix.Area();

      double score;
      AngularCoordinate ang;

      if (resolution_level == starting_resolution_level) {
        resolve_pix.clear();
        previous_pix.clear();

        uint32_t nx_pix;
        if ((x_max < x_min) && (x_min > nx/2)) {
          nx_pix = nx - x_min + x_max + 1;
        } else {
          nx_pix = x_max - x_min + 1;
        }

        for (uint32_t y=y_min;y<=y_max;y++) {
          for (uint32_t m=0,x=x_min;m<nx_pix;m++,x++) {
            if (x==nx) x = 0;
            tmp_pix.SetPixnumFromXY(x,y);

            score = _ScorePixel(bound, tmp_pix);

            if (score < -0.99999) {
              tmp_pix.SetWeight(weight);
              kept_pix.push_back(tmp_pix);
	      pixel_area += unit_area;
              n_keep++;
            } else {
              if (score < -0.00001) {
                tmp_pix.SetWeight(score);
                resolve_pix.push_back(tmp_pix);
              }
            }
            previous_pix.push_back(tmp_pix);
          }
        }
      } else {
        if (resolve_pix.size() == 0) {
          std::cout << "Missed all pixels in initial search; trying again...\n";
          for (PixelIterator iter=previous_pix.begin();
               iter!=previous_pix.end();++iter) {
            PixelVector sub_pix;
            iter->SubPix(resolution,sub_pix);
            for (PixelIterator sub_iter=sub_pix.begin();
                 sub_iter!=sub_pix.end();++sub_iter)
              resolve_pix.push_back(*sub_iter);
          }
        }

        previous_pix.clear();

	previous_pix.reserve(resolve_pix.size());
        for (PixelIterator iter=resolve_pix.begin();
             iter!=resolve_pix.end();++iter) previous_pix.push_back(*iter);

        resolve_pix.clear();

        uint32_t x_min, x_max, y_min, y_max;

        for (PixelIterator iter=previous_pix.begin();
             iter!=previous_pix.end();++iter) {

          iter->SubPix(resolution,x_min,x_max,y_min,y_max);

          for (uint32_t y=y_min;y<=y_max;y++) {
            for (uint32_t x=x_min;x<=x_max;x++) {
              tmp_pix.SetPixnumFromXY(x,y);

              score = _ScorePixel(bound, tmp_pix);

              if (score < -0.99999) {
                tmp_pix.SetWeight(weight);
		kept_pix.push_back(tmp_pix);
		pixel_area += unit_area;
                n_keep++;
              } else {
                if (score < -0.00001) {
                  tmp_pix.SetWeight(score);
                  resolve_pix.push_back(tmp_pix);
                }
              }
            }
          }
        }
      }
    }

    previous_pix.clear();

    if (bound.Area() > pixel_area) {
      sort(resolve_pix.begin(), resolve_pix.end(), Pixel::WeightedOrder);

      uint32_t n=0;
      double ur_weight = resolve_pix[n].Weight();
      double unit_area = resolve_pix[n].Area();
      while ((n < resolve_pix.size()) &&
             ((bound.Area() > pixel_area) ||
              DoubleEQ(resolve_pix[n].Weight(), ur_weight))) {
        ur_weight = resolve_pix[n].Weight();
        resolve_pix[n].SetWeight(weight);
	kept_pix.push_back(resolve_pix[n]);
	pixel_area += unit_area;
        n++;
      }
    }

    Initialize(kept_pix);

    pixelized_map = true;
  }

  return pixelized_map;
}

uint8_t Map::_FindStartingResolutionLevel(double bound_area) {
  if (bound_area < 10.0*Pixel::PixelArea(MaxPixelResolution)) {
    return 0;
  }

  uint32_t starting_resolution = HPixResolution;

  // We want to start things off with the coarsest possible resolution to
  // save time, but we have to be careful that we're not so coarse that we
  // miss parts of the footprint.  This finds the resolution that has pixels
  // about 1/100th the area of the footprint.
  while ((bound_area/Pixel::PixelArea(starting_resolution) <= 100.0) &&
	 (starting_resolution < MaxPixelResolution)) starting_resolution *= 2;

  return Pixel::ResolutionToLevel(starting_resolution);
}

bool Map::_FindXYBounds(const uint8_t resolution_level,
			GeometricBound& bound,
			uint32_t& x_min, uint32_t&x_max,
			uint32_t& y_min, uint32_t&y_max) {

  uint32_t resolution = Pixel::LevelToResolution(resolution_level);
  uint32_t nx = Nx0*resolution, ny = Ny0*resolution;

  Pixel::AreaIndex(resolution,
		   bound.LambdaMin(), bound.LambdaMax(),
		   bound.EtaMin(), bound.EtaMax(),
		   x_min, x_max, y_min, y_max);

  if ((bound.EtaMax() - bound.EtaMin() > 270.0) && bound.ContinuousBounds()) {
    x_min = 0;
    x_max = nx - 1;
  }

  // std::cout << "Starting bounds:\n\tX: " << x_min << " - " << x_max <<
  // ", Y: " << y_min << " - " << y_max << "\n";

  // Checking top border
  bool found_pixel = true;
  bool boundary_failure = false;

  Pixel tmp_pix;
  tmp_pix.SetResolution(resolution);

  uint8_t n_iter = 0;
  uint8_t max_iter = 20;
  uint32_t nx_pix = 0;
  if ((x_max < x_min) && (x_min > nx/2)) {
    nx_pix = nx - x_min + x_max + 1;
  } else {
    nx_pix = x_max - x_min + 1;
  }

  while (found_pixel && n_iter < max_iter) {
    found_pixel = false;
    uint32_t y = y_max;

    for (uint32_t m=0,x=x_min;m<nx_pix;m++,x++) {
      if (x == nx) x = 0;

      tmp_pix.SetPixnumFromXY(x,y);

      // This if statement checks positions within the pixel against the
      // footprint bound.
      if (DoubleLT(_ScorePixel(bound, tmp_pix), 0.0)) {
	// std::cout << m << "," << y << ": " <<
	// tmp_pix.Lambda() << ", " << tmp_pix.Eta() << ": " <<
	// _ScorePixel(bound, tmp_pix) << "\n";
        found_pixel = true;
        m = nx_pix + 1;
      }
    }

    if (found_pixel) {
      // The exception to that case is if we've already reached the maximum
      // y index for the pixels.  In that case, we're just done.
      if (y_max < ny - 1) {
        y_max++;
      } else {
        found_pixel = false;
      }
    }
    n_iter++;
  }
  if (n_iter == max_iter) boundary_failure = true;
  // if (boundary_failure) std::cout << "\t\tBoundary failure on top bound\n";

  // Checking bottom border
  found_pixel = true;
  n_iter = 0;
  while (!boundary_failure && found_pixel && n_iter < max_iter) {
    found_pixel = false;
    uint32_t y = y_min;

    for (uint32_t m=0,x=x_min;m<nx_pix;m++,x++) {
      if (x == nx) x = 0;

      tmp_pix.SetPixnumFromXY(x,y);

      // This if statement checks positions within the pixel against the
      // footprint bound.
      if (_ScorePixel(bound, tmp_pix) < -0.000001) {
        found_pixel = true;
        m = nx_pix + 1;
      }
    }

    if (found_pixel) {
      // The exception to that case is if we've already reached the minimum
      // y index for the pixels.  In that case, we're just done.
      if (y_min > 0) {
        y_min--;
      } else {
        found_pixel = false;
      }
    }
    n_iter++;
  }

  if (n_iter == max_iter) boundary_failure = true;
  // if (boundary_failure) std::cout << "\t\tBoundary failure on lower bound\n";

  // Checking left border
  found_pixel = true;
  n_iter = 0;

  while (!boundary_failure && found_pixel && n_iter < max_iter && nx_pix < nx) {
    found_pixel = false;
    uint32_t x = x_min;

    for (uint32_t y=y_min;y<=y_max;y++) {

      tmp_pix.SetPixnumFromXY(x,y);

      // This if statement checks positions within the pixel against the
      // footprint bound.
      if (_ScorePixel(bound, tmp_pix) < -0.000001) {
        found_pixel = true;
        y = y_max + 1;
      }
    }

    if (found_pixel) {
      if (x_min == 0) {
        x_min = nx - 1;
      } else {
        x_min--;
      }
    }
    n_iter++;
  }
  if (n_iter == max_iter) boundary_failure = true;
  // if (boundary_failure) std::cout << "\t\tBoundary failure on left bound\n";

  // Checking right border
  found_pixel = true;
  n_iter = 0;
  if ((x_max < x_min) && (x_min > nx/2)) {
    nx_pix = nx - x_min + x_max + 1;
  } else {
    nx_pix = x_max - x_min + 1;
  }
  while (!boundary_failure && found_pixel && n_iter < max_iter && nx_pix < nx) {
    found_pixel = false;
    uint32_t x = x_max;

    if ((x_max < x_min) && (x_min > nx/2)) {
      nx_pix = nx - x_min + x_max + 1;
    } else {
      nx_pix = x_max - x_min + 1;
    }

    for (uint32_t y=y_min;y<=y_max;y++) {

      tmp_pix.SetPixnumFromXY(x,y);

      // This if statement checks positions within the pixel against the
      // footprint bound.
      if (_ScorePixel(bound, tmp_pix) < -0.000001) {
        found_pixel = true;
        y = y_max + 1;
      }
    }

    if (found_pixel) {
      if (x_max == nx - 1) {
        x_max = 0;
      } else {
        x_max++;
      }
    }
    n_iter++;
  }

  if (n_iter == max_iter) boundary_failure = true;
  // if (boundary_failure) std::cout << "\t\tBoundary failure on right bound\n";

  // std::cout << "Final bounds:\n\tX: " << x_min << " - " << x_max <<
  // ", Y: " << y_min << " - " << y_max << "\n";

  return !boundary_failure;
}

double Map::_ScorePixel(GeometricBound& bound, Pixel& pix) {

  double inv_nx = 1.0/static_cast<double>(Nx0*pix.Resolution());
  double inv_ny = 1.0/static_cast<double>(Ny0*pix.Resolution());
  double x = static_cast<double>(pix.PixelX());
  double y = static_cast<double>(pix.PixelY());

  double lammid = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.5)*inv_ny);
  double lammin = 90.0 - RadToDeg*acos(1.0-2.0*(y+1.0)*inv_ny);
  double lammax = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.0)*inv_ny);
  double lam_quart = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.75)*inv_ny);
  double lam_three = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.25)*inv_ny);

  double etamid = RadToDeg*(2.0*Pi*(x+0.5))*inv_nx + EtaOffSet;
  if (DoubleGE(etamid, 180.0)) etamid -= 360.0;
  if (DoubleLE(etamid, -180.0)) etamid += 360.0;

  double etamin = RadToDeg*(2.0*Pi*(x+0.0))*inv_nx + EtaOffSet;
  if (DoubleGE(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLE(etamin, -180.0)) etamin += 360.0;

  double etamax = RadToDeg*(2.0*Pi*(x+1.0))*inv_nx + EtaOffSet;
  if (DoubleGE(etamax, 180.0)) etamax -= 360.0;
  if (DoubleLE(etamax, -180.0)) etamax += 360.0;

  double eta_quart = RadToDeg*(2.0*Pi*(x+0.25))*inv_nx + EtaOffSet;
  if (DoubleGE(eta_quart, 180.0)) eta_quart -= 360.0;
  if (DoubleLE(eta_quart, -180.0)) eta_quart += 360.0;

  double eta_three = RadToDeg*(2.0*Pi*(x+0.75))*inv_nx + EtaOffSet;
  if (DoubleGE(eta_three, 180.0)) eta_three -= 360.0;
  if (DoubleLE(eta_three, -180.0)) eta_three += 360.0;

  double score = 0.0;

  AngularCoordinate ang(lammid, etamid, AngularCoordinate::Survey);
  if (bound.CheckPoint(ang)) score -= 4.0;

  ang.SetSurveyCoordinates(lam_quart,etamid);
  if (bound.CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,etamid);
  if (bound.CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lammid,eta_quart);
  if (bound.CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lammid,eta_quart);
  if (bound.CheckPoint(ang)) score -= 3.0;

  ang.SetSurveyCoordinates(lam_quart,eta_quart);
  if (bound.CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,eta_quart);
  if (bound.CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_quart,eta_three);
  if (bound.CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,eta_three);
  if (bound.CheckPoint(ang)) score -= 3.0;

  ang.SetSurveyCoordinates(lammid,etamax);
  if (bound.CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammid,etamin);
  if (bound.CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammax,etamid);
  if (bound.CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammin,etamid);
  if (bound.CheckPoint(ang)) score -= 2.0;

  ang.SetSurveyCoordinates(lammax,etamax);
  if (bound.CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammax,etamin);
  if (bound.CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammin,etamax);
  if (bound.CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammin,etamin);
  if (bound.CheckPoint(ang)) score -= 1.0;

  return score/40.0;
}

bool Map::IngestMap(PixelVector& pix, bool destroy_copy) {
  for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) {
    uint32_t k = iter->Superpixnum();
    sub_map_[k].AddPixel(*iter);
    sub_map_[k].SetUnsorted();
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
      sub_map_[iter->Superpixnum()].SetUnsorted();
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

  return ExcludeMap(stomp_map, destroy_copy);
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
  min_level_ = MaxPixelLevel;
  max_level_ = HPixLevel;
  min_weight_ = 1.0e30;
  max_weight_ = -1.0e30;
  ClearRegions();

  for (uint32_t resolution=HPixResolution;
       resolution<=MaxPixelResolution;resolution*=2)
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

uint32_t Map::MinResolution() {
  return Pixel::LevelToResolution(min_level_);
}

uint32_t Map::MinResolution(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].MinResolution() : 0);
}

uint32_t Map::MaxResolution() {
  return Pixel::LevelToResolution(max_level_);
}

uint32_t Map::MaxResolution(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].MaxResolution() : 0);
}

uint8_t Map::MinLevel() {
  return min_level_;
}

uint8_t Map::MinLevel(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].MinLevel() : 0);
}

uint8_t Map::MaxLevel() {
  return max_level_;
}

uint8_t Map::MaxLevel(uint32_t superpixnum) {
  return (superpixnum < MaxSuperpixnum ?
	  sub_map_[superpixnum].MaxLevel() : 0);
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

uint32_t Map::PixelCount(uint32_t resolution) {
  return (!(resolution % 2) ? pixel_count_[resolution] : 0);
}


} // end namespace Stomp
