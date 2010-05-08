// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the abstract BaseMap class that serves as the
// basis for all of the *Map objects.  BaseMap sets out the basic functionality
// that all of the *Map classes need to describe a given region on the sky and
// do some basic internal maintenance.  Additionally, BaseMap provides a
// common set of methods for dividing that area up into nearly equal-area,
// similarly-shaped regions.  This functionality is the basis for calculating
// jack-knife errors for our various statistical analyses.

#include "stomp_core.h"
#include "stomp_base_map.h"

namespace Stomp {

Section::Section() {
  stripe_min_ = -1;
  stripe_max_ = -1;
}

Section::~Section() {
  stripe_min_ = -1;
  stripe_max_ = -1;
}

void Section::SetMinStripe(uint32_t stripe) {
  stripe_min_ = stripe;
}

void Section::SetMaxStripe(uint32_t stripe) {
  stripe_max_ = stripe;
}

uint32_t Section::MinStripe() {
  return stripe_min_;
}

uint32_t Section::MaxStripe() {
  return stripe_max_;
}

RegionMap::RegionMap() {
  ClearRegions();
};

RegionMap::~RegionMap() {
  ClearRegions();
};

uint16_t RegionMap::InitializeRegions(BaseMap* stomp_map, uint16_t n_region,
				      uint32_t region_resolution) {
  ClearRegions();

  // If we have the default value for the resolution, we need to attempt to
  // find a reasonable value for the resolution based on the area.  We want to
  // shoot for something along the lines of 50 pixels per region to give us a
  // fair chance of getting equal areas without using too many pixels.
  if (region_resolution == 0) {
    double target_area = stomp_map->Area()/(50*n_region);
    region_resolution = Stomp::HPixResolution;
    while ((Pixel::PixelArea(region_resolution) > target_area) &&
	   (region_resolution < 1024)) region_resolution <<= 1;
    // std::cout << "Automatically setting region resolution to " <<
    // region_resolution << " based on requested n_region\n";
  }

  if (region_resolution > 256) {
    std::cout <<
      "WARNING: Attempting to generate region map with resolution " <<
      "above 256!\n";
    std::cout << "This may end badly.\n";
    if (region_resolution > 2048) {
      std::cout <<
	"FAIL: Ok, the resolution is above 2048.  Can't do this.  Try again\n";
      exit(1);
    }
  }

  if (region_resolution > stomp_map->MaxResolution()) {
    std::cout << "WARNING: Re-setting region map resolution to " <<
      stomp_map->MaxResolution() << " to satisfy input map limits.\n";
    region_resolution = stomp_map->MaxResolution();
  }

  region_resolution_ = region_resolution;
  //std::cout << "Generating regionation map at " << region_resolution_ << "\n";

  PixelVector region_pix;
  stomp_map->Coverage(region_pix, region_resolution_);

  // std::cout << "Creating region map with " << region_pix.size() <<
  // " pixels...\n";

  if (static_cast<uint32_t>(n_region) > region_pix.size()) {
    std::cout << "WARNING: Exceeded maximum possible regions.  Setting to " <<
      region_pix.size() << " regions.\n";
    n_region = static_cast<uint32_t>(region_pix.size());
  }

  if (static_cast<uint32_t>(n_region) == region_pix.size()) {
    int16_t i=0;
    for (PixelIterator iter=region_pix.begin();iter!=region_pix.end();++iter) {
      region_map_[iter->Pixnum()] = i;
      i++;
    }
    std::cout <<
      "\tWARNING: Number of regions matches number of regionation pixels.\n";
    std::cout << "\tThis will be dead easy, " <<
      "but won't guarantee an equal area solution...\n";
  } else {
    // std::cout << "\tBreaking up " << stomp_map->Area() <<
    // " square degrees into " << n_region << " equal-area pieces.\n";

    std::vector<uint32_t> tmp_stripe;
    tmp_stripe.reserve(region_pix.size());
    for (PixelIterator iter=region_pix.begin();
	 iter!=region_pix.end();++iter) {
      tmp_stripe.push_back(iter->Stripe(region_resolution_));
    }

    sort(tmp_stripe.begin(), tmp_stripe.end());
    std::vector<uint32_t> stripe;
    stripe.push_back(tmp_stripe[0]);
    for (uint32_t i=1;i<tmp_stripe.size();i++)
      if (tmp_stripe[i] != tmp_stripe[i-1]) stripe.push_back(tmp_stripe[i]);

    tmp_stripe.clear();

    sort(stripe.begin(), stripe.end());

    std::vector<Section> super_section;

    Section tmp_section;
    tmp_section.SetMinStripe(stripe[0]);
    tmp_section.SetMaxStripe(stripe[0]);

    super_section.push_back(tmp_section);

    for (uint32_t i=1,j=0;i<stripe.size();i++) {
      if (stripe[i] == stripe[i-1] + 1) {
        super_section[j].SetMaxStripe(stripe[i]);
      } else {
        tmp_section.SetMinStripe(stripe[i]);
        tmp_section.SetMaxStripe(stripe[i]);
        super_section.push_back(tmp_section);
        j++;
      }
    }

    double region_length = sqrt(stomp_map->Area()/n_region);
    uint8_t region_width =
        static_cast<uint8_t>(region_length*Stomp::Nx0*region_resolution_/360.0);
    if (region_width == 0) region_width = 1;

    std::vector<Section> section;

    int32_t j = -1;
    for (std::vector<Section>::iterator iter=super_section.begin();
         iter!=super_section.end();++iter) {

      for (uint32_t stripe_iter=iter->MinStripe(),section_iter=region_width;
           stripe_iter<=iter->MaxStripe();stripe_iter++) {
        if (section_iter == region_width) {
          tmp_section.SetMinStripe(stripe_iter);
          tmp_section.SetMaxStripe(stripe_iter);
          section.push_back(tmp_section);
          section_iter = 1;
          j++;
        } else {
          section[j].SetMaxStripe(stripe_iter);
          section_iter++;
        }
      }
    }

    double region_area = 0.0, running_area = 0.0;
    double unit_area = Pixel::PixelArea(region_resolution_);
    uint32_t n_pixel = 0;
    int16_t region_iter = 0;
    double mean_area = stomp_map->Area()/region_pix.size();
    double area_break = stomp_map->Area()/n_region;

    // std::cout << "\tAssigning areas...\n\n";
    // std::cout << "\tSample  Pixels  Unmasked Area  Masked Area\n";
    // std::cout << "\t------  ------  -------------  -----------\n";
    for (std::vector<Section>::iterator section_iter=section.begin();
         section_iter!=section.end();++section_iter) {

      for (PixelIterator iter=region_pix.begin();
	   iter!=region_pix.end();++iter) {
        if ((iter->Stripe(region_resolution_) >=
	     section_iter->MinStripe()) &&
            (iter->Stripe(region_resolution_) <=
	     section_iter->MaxStripe())) {
          if ((region_area + 0.75*mean_area < area_break*(region_iter+1)) ||
	      (region_iter == n_region-1)) {
            region_area += iter->Weight()*unit_area;
            region_map_[iter->Pixnum()] = region_iter;
            running_area += iter->Weight()*unit_area;
            n_pixel++;
          } else {
	    // std::cout << "\t" << region_iter << "\t" << n_pixel << "\t" <<
	    // n_pixel*unit_area << "\t\t" << running_area << "\n";
	    region_area_[region_iter] = running_area;

            region_iter++;
            region_area += iter->Weight()*unit_area;
            region_map_[iter->Pixnum()] = region_iter;
            running_area = iter->Weight()*unit_area;
            n_pixel = 1;
          }
        }
      }
    }
    region_area_[region_iter] = running_area;
    // std::cout << "\t" << region_iter << "\t" << n_pixel << "\t" <<
    // n_pixel*unit_area << "\t\t" << running_area << "\n";
  }

  std::vector<uint32_t> region_count_check;

  for (uint16_t i=0;i<n_region;i++) region_count_check.push_back(0);

  for (RegionIterator iter=region_map_.begin();
       iter!=region_map_.end();++iter) {
    if (iter->second < n_region) {
      region_count_check[iter->second]++;
    } else {
      std::cout << "FAIL: Encountered illegal region index: " <<
	iter->second << "\nBailing...\n";
      exit(2);
    }
  }

  n_region_ = static_cast<uint16_t>(region_area_.size());

  return n_region_;
}

bool RegionMap::InitializeRegions(BaseMap* base_map, BaseMap& stomp_map) {
  bool initialized_region_map = true;
  if (!region_map_.empty()) region_map_.clear();

  region_resolution_ = stomp_map.RegionResolution();
  n_region_ = stomp_map.NRegion();

  // Iterate through the current BaseMap to find the region value for each
  // pixel.  If the node is not present in the input map, then
  // we bail and return false.
  PixelVector coverage_pix;
  base_map->Coverage(coverage_pix, stomp_map.RegionResolution());

  for (PixelIterator iter=coverage_pix.begin();
       iter!=coverage_pix.end();++iter) {
    int16_t region = stomp_map.Region(iter->SuperPix(region_resolution_));
    if (region != -1) {
      region_map_[iter->SuperPix(region_resolution_)] = region;
    } else {
      initialized_region_map = false;
      iter = coverage_pix.end();
    }
  }

  if (!initialized_region_map) {
    region_resolution_ = 0;
    n_region_ = -1;
  }

  return initialized_region_map;
}

int16_t RegionMap::FindRegion(AngularCoordinate& ang) {
  Pixel tmp_pix(ang, region_resolution_, 1.0);

  return (region_map_.find(tmp_pix.Pixnum()) != region_map_.end() ?
	  region_map_[tmp_pix.Pixnum()] : -1);
}

void RegionMap::ClearRegions() {
  region_map_.clear();
  n_region_ = 0;
  region_resolution_ = 0;
}

int16_t RegionMap::Region(uint32_t region_idx) {
  return (region_map_.find(region_idx) != region_map_.end() ?
	  region_map_[region_idx] : -1);
}

void RegionMap::RegionArea(int16_t region_index, PixelVector& pix) {
  pix.clear();

  for (RegionIterator iter=Begin();iter!=End();++iter) {
    if (iter->second == region_index) {
      Pixel tmp_pix(Resolution(), iter->first, 1.0);
      pix.push_back(tmp_pix);
    }
  }
}

double RegionMap::RegionArea(int16_t region) {
  return (region_area_.find(region) != region_area_.end() ?
	  region_area_[region] : 0.0);
}

uint16_t RegionMap::NRegion() {
  return n_region_;
}

uint32_t RegionMap::Resolution() {
  return region_resolution_;
}

bool RegionMap::Initialized() {
  return (n_region_ > 0 ? true : false);
}

RegionIterator RegionMap::Begin() {
  return region_map_.begin();
}

RegionIterator RegionMap::End() {
  return region_map_.end();
}

BaseMap::BaseMap() {
  ClearRegions();
}

BaseMap::~BaseMap() {
  ClearRegions();
}

void BaseMap::Coverage(PixelVector& superpix, uint32_t resolution,
		       bool calculate_fraction) {
  superpix.clear();
}

double BaseMap::FindUnmaskedFraction(Pixel& pix) {
  return 0.0;
}

int8_t BaseMap::FindUnmaskedStatus(Pixel& pix) {
  return 0;
}

bool BaseMap::Empty() {
  return true;
}

void BaseMap::Clear() {
  ClearRegions();
}

uint32_t BaseMap::Size() {
  return 0;
}

double BaseMap::Area() {
  return 0.0;
}

uint32_t BaseMap::MinResolution() {
  return Stomp::HPixResolution;
}

uint32_t BaseMap::MaxResolution() {
  return Stomp::MaxPixelResolution;
}

uint8_t BaseMap::MinLevel() {
  return Stomp::HPixLevel;
}

uint8_t BaseMap::MaxLevel() {
  return Stomp::MaxPixelLevel;
}

uint16_t BaseMap::InitializeRegions(uint16_t n_regions,
				    uint32_t region_resolution) {
  return region_map_.InitializeRegions(this, n_regions, region_resolution);
}

bool BaseMap::InitializeRegions(BaseMap& base_map) {
  return region_map_.InitializeRegions(this, base_map);
}

int16_t BaseMap::FindRegion(AngularCoordinate& ang) {
  return region_map_.FindRegion(ang);
}

void BaseMap::ClearRegions() {
  region_map_.ClearRegions();
}

void BaseMap::RegionArea(int16_t region, PixelVector& pix) {
  region_map_.RegionArea(region, pix);
}

int16_t BaseMap::Region(uint32_t region_idx) {
  return region_map_.Region(region_idx);
}

double BaseMap::RegionArea(int16_t region) {
  return region_map_.RegionArea(region);
}

uint16_t BaseMap::NRegion() {
  return region_map_.NRegion();
}

uint32_t BaseMap::RegionResolution() {
  return region_map_.Resolution();
}

bool BaseMap::RegionsInitialized() {
  return region_map_.Initialized();
}

RegionIterator BaseMap::RegionBegin() {
  return region_map_.Begin();
}

RegionIterator BaseMap::RegionEnd() {
  return region_map_.End();
}

} // end namespace Stomp
