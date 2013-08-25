// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the ScalarMap class.  Unlike Maps, the primary
// goal here is to encode a scalar field over some area of the sky.  As such,
// we sacrifice some degree of precision in describing the exact area of the
// field and we use a uniform sampling of the field across the area in
// question.  This makes the class ideal for calculating angular correlation
// functions on the encoded field.

#include "stomp_core.h"
#include "stomp_scalar_map.h"
#include "stomp_map.h"
#include "stomp_angular_correlation.h"

namespace Stomp {

ScalarMap::ScalarMap() {
  area_ = 0.0;
  resolution_ = 0;
  total_points_ = 0;
  mean_intensity_ = 0.0;
  total_intensity_ = 0.0;
  if (!pix_.empty()) pix_.clear();
  ClearRegions();
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;
  map_type_ = ScalarField;
}

ScalarMap::ScalarMap(Map& stomp_map, uint32_t input_resolution,
		     ScalarMapType scalar_map_type,
		     double min_unmasked_fraction,
		     bool use_map_weight_as_intensity) {
  resolution_ = input_resolution;
  unmasked_fraction_minimum_ = min_unmasked_fraction;
  map_type_ = scalar_map_type;

  if (use_map_weight_as_intensity && !(map_type_ == ScalarField)) {
    std::cout <<
      "Stomp::ScalarMap::ScalarMap - " <<
      "WARNING: Converting MapType to ScalarField to sample " <<
      "input Map Weight\n";
    map_type_ = ScalarField;
  };

  PixelVector superpix;
  stomp_map.Coverage(superpix);

  for (PixelIterator iter=superpix.begin();iter!=superpix.end();++iter) {
    PixelVector sub_pix;
    iter->SubPix(resolution_,sub_pix);

    for (PixelIterator sub_iter=sub_pix.begin();
         sub_iter!=sub_pix.end();++sub_iter) {
      double unmasked_fraction = stomp_map.FindUnmaskedFraction(*sub_iter);
      double initial_intensity = 0.0;
      if (unmasked_fraction > unmasked_fraction_minimum_) {
        if (use_map_weight_as_intensity)
          initial_intensity = stomp_map.FindAverageWeight(*sub_iter);
        ScalarPixel tmp_pix(sub_iter->PixelX(), sub_iter->PixelY(),
			    sub_iter->Resolution(), unmasked_fraction,
			    initial_intensity, 0);
        pix_.push_back(tmp_pix);
      }
    }
  }

  pix_.resize(pix_.size());

  sort(pix_.begin(), pix_.end(), Pixel::LocalOrder);
  mean_intensity_ = 0.0;
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;

  area_ = 0.0;
  total_intensity_ = 0.0;
  total_points_ = 0;
  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    area_ += iter->Area()*iter->Weight();
    total_intensity_ += iter->Intensity();
    total_points_ += iter->NPoints();
  }
}

ScalarMap::ScalarMap(ScalarMap& scalar_map,
		     uint32_t input_resolution,
		     double min_unmasked_fraction) {
  if (input_resolution > scalar_map.Resolution()) {
    std::cout << "Stomp::ScalarMap::ScalarMap - " <<
      "Cannot make higher resolution density map " <<
      "by resampling. Exiting.\n";
    exit(1);
  }

  resolution_ = input_resolution;
  unmasked_fraction_minimum_ = min_unmasked_fraction;
  map_type_ = scalar_map.MapType();

  if (scalar_map.Resolution() == resolution_) {
    pix_.reserve(scalar_map.Size());
    for (ScalarIterator iter=scalar_map.Begin();
	 iter!=scalar_map.End();++iter) pix_.push_back(*iter);
  } else {
    PixelVector superpix;
    scalar_map.Coverage(superpix, HPixResolution, false);

    uint32_t x_min, x_max, y_min, y_max;
    for (PixelIterator iter=superpix.begin();iter!=superpix.end();++iter) {
      iter->SubPix(resolution_, x_min, x_max, y_min, y_max);
      for (uint32_t y=y_min;y<=y_max;y++) {
	for (uint32_t x=x_min;x<=x_max;x++) {
          ScalarPixel tmp_pix(x, y, resolution_);
	  scalar_map.Resample(tmp_pix);
	  if (tmp_pix.Weight() > unmasked_fraction_minimum_)
	    pix_.push_back(tmp_pix);
        }
      }
    }
  }

  pix_.resize(pix_.size());

  sort(pix_.begin(), pix_.end(), Pixel::LocalOrder);
  mean_intensity_ = 0.0;
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;

  area_ = 0.0;
  total_intensity_ = 0.0;
  total_points_ = 0;
  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    area_ += iter->Area()*iter->Weight();
    total_intensity_ += iter->Intensity();
    total_points_ += iter->NPoints();
  }
}

ScalarMap::ScalarMap(ScalarVector& pix,
		     ScalarMapType scalar_map_type,
		     double min_unmasked_fraction) {

  resolution_ = pix[0].Resolution();
  unmasked_fraction_minimum_ = min_unmasked_fraction;
  map_type_ = scalar_map_type;

  pix_.reserve(pix.size());

  area_ = 0.0;
  total_intensity_ = 0.0;
  total_points_ = 0;
  for (ScalarIterator iter=pix.begin();iter!=pix.end();++iter) {
    if (iter->Resolution() != resolution_) {
      std::cout << "Stomp::ScalarMap::ScalarMap - " <<
	"Incompatible resolutions in ScalarPixel list.  Exiting.\n";
      exit(2);
    }

    area_ += iter->Area()*iter->Weight();
    total_intensity_ += iter->Intensity();
    total_points_ += iter->NPoints();
    pix_.push_back(*iter);
  }

  sort(pix_.begin(), pix_.end(), Pixel::LocalOrder);
  mean_intensity_ = 0.0;
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;
}

ScalarMap::ScalarMap(Map& stomp_map,
		     AngularCoordinate& center, double theta_max,
		     uint32_t input_resolution,
		     ScalarMapType scalar_map_type,
		     double min_unmasked_fraction,
		     double theta_min) {

  resolution_ = input_resolution;
  unmasked_fraction_minimum_ = min_unmasked_fraction;
  map_type_ = scalar_map_type;

  Pixel pix(center, resolution_);

  PixelVector pixVec;
  pix.WithinAnnulus(theta_min, theta_max, pixVec);

  for (PixelIterator iter=pixVec.begin();iter!=pixVec.end();++iter) {
    double unmasked_fraction = stomp_map.FindUnmaskedFraction(*iter);
    if (unmasked_fraction > unmasked_fraction_minimum_) {
      pix_.push_back(ScalarPixel(iter->PixelX(), iter->PixelY(),
                                 resolution_, unmasked_fraction));
    }
  }

  pix_.resize(pix_.size());

  sort(pix_.begin(), pix_.end(), Pixel::LocalOrder);
  mean_intensity_ = 0.0;
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;

  area_ = 0.0;
  total_intensity_ = 0.0;
  total_points_ = 0;
  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter)
    area_ += iter->Area()*iter->Weight();
}

ScalarMap::~ScalarMap() {
  area_ = 0.0;
  resolution_ = 0;
  mean_intensity_ = 0.0;
  total_intensity_ = 0.0;
  if (!pix_.empty()) pix_.clear();
  ClearRegions();
}

bool ScalarMap::Read(const std::string& InputFile,
		double min_unmasked_fraction) {
  Clear();

  std::ifstream input_file(InputFile.c_str());

  uint32_t hpixnum, superpixnum;
  double unmask, intensity;
  int npoints;
  int resolution;  // have to use an int here because of old map formats
  bool found_file = false;

  uint32_t x, y;

  ScalarVector pix;

  if (input_file) {
  	found_file = true;
    while (!input_file.eof()) {
    	input_file >> hpixnum >> superpixnum >> resolution >> unmask
    						 >> intensity >> npoints;
      if (!input_file.eof() && (resolution % 2 == 0) &&
      		(resolution > 0)) {
      	Pixel::HPix2XY(static_cast<uint32_t>(resolution), hpixnum, superpixnum,
      								 x, y);
      	pix.push_back(
      			ScalarPixel(x, y, static_cast<uint32_t>(resolution),
      			            unmask, intensity, npoints));
      }
    }

    input_file.close();

  } else {
    std::cout << "Stomp::ScalarMap::Read - " << InputFile <<
      " does not exist!.  No Map ingested\n";
  }

  resolution_ = pix[0].Resolution();
  unmasked_fraction_minimum_ = min_unmasked_fraction;
  map_type_ = ScalarField;

  pix_.reserve(pix.size());

  area_ = 0.0;
  total_intensity_ = 0.0;
  total_points_ = 0;
  for (ScalarIterator iter=pix.begin();iter!=pix.end();++iter) {
    if (iter->Resolution() != resolution_) {
      std::cout << "Stomp::ScalarMap::ScalarMap - " <<
          "Incompatible resolutions in input Map file.  Exiting.\n";
      exit(2);
    }

    area_ += iter->Area()*iter->Weight();
    total_intensity_ += iter->Intensity();
    total_points_ += iter->NPoints();
    pix_.push_back(*iter);
  }

  sort(pix_.begin(), pix_.end(), Pixel::LocalOrder);
  mean_intensity_ = 0.0;
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;

  return found_file;
}

void ScalarMap::SetResolution(uint32_t resolution) {
  Clear();
  resolution_ = resolution;
}

void ScalarMap::InitializeFromMap(Map& stomp_map, uint32_t input_resolution,
				  bool use_map_weight_as_intensity) {
  uint32_t current_resolution = resolution_;
  Clear();

  if (input_resolution != 0) {
    SetResolution(input_resolution);
  } else {
    SetResolution(current_resolution);
  }

  if (use_map_weight_as_intensity && !(map_type_ != ScalarField)) {
    std::cout <<
      "Stomp::ScalarMap::InitializeFromMap - " <<
      "WARNING: Converting MapType to ScalarField to sample " <<
      "input Map Weight\n";
    map_type_ = ScalarField;
  };

  PixelVector superpix;
  stomp_map.Coverage(superpix, HPixResolution, false);

  for (PixelIterator iter=superpix.begin();iter!=superpix.end();++iter) {
    PixelVector sub_pix;
    iter->SubPix(resolution_,sub_pix);

    for (PixelIterator sub_iter=sub_pix.begin();
         sub_iter!=sub_pix.end();++sub_iter) {
      double unmasked_fraction = stomp_map.FindUnmaskedFraction(*sub_iter);
      double initial_intensity = 0.0;
      if (unmasked_fraction > unmasked_fraction_minimum_) {
        if (use_map_weight_as_intensity)
          initial_intensity = stomp_map.FindAverageWeight(*sub_iter);
	ScalarPixel tmp_pix(sub_iter->PixelX(), sub_iter->PixelY(),
			    sub_iter->Resolution(), unmasked_fraction,
			    initial_intensity, 0);
	pix_.push_back(tmp_pix);
      }
    }
  }

  pix_.resize(pix_.size());

  sort(pix_.begin(), pix_.end(), Pixel::LocalOrder);
  mean_intensity_ = 0.0;
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;

  area_ = 0.0;
  total_intensity_ = 0.0;
  total_points_ = 0;
  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    area_ += iter->Area()*iter->Weight();
    total_intensity_ += iter->Intensity();
  }
}

void ScalarMap::InitializeFromScalarMap(ScalarMap& scalar_map,
					uint32_t input_resolution) {
  if (input_resolution > scalar_map.Resolution()) {
    std::cout << "Stomp::ScalarMap::InitializeFromScalarMap - " <<
      "Cannot make higher resolution density map " <<
      "by resampling. Exiting.\n";
    exit(1);
  }

  uint32_t current_resolution = resolution_;
  Clear();

  if (input_resolution != 0) {
    SetResolution(input_resolution);
  } else {
    SetResolution(current_resolution);
  }

  map_type_ = scalar_map.MapType();

  if (scalar_map.Resolution() == resolution_) {
    pix_.reserve(scalar_map.Size());
    for (ScalarIterator iter=scalar_map.Begin();
	 iter!=scalar_map.End();++iter) pix_.push_back(*iter);
  } else {
    PixelVector superpix;
    scalar_map.Coverage(superpix, HPixResolution, false);

    uint32_t x_min, x_max, y_min, y_max;

    for (PixelIterator iter=superpix.begin();iter!=superpix.end();++iter) {
      iter->SubPix(resolution_,x_min,x_max,y_min,y_max);
      for (uint32_t y=y_min;y<=y_max;y++) {
	for (uint32_t x=x_min;x<=x_max;x++) {
          ScalarPixel tmp_pix(x, y, resolution_);
	  scalar_map.Resample(tmp_pix);
	  if (tmp_pix.Weight() > unmasked_fraction_minimum_)
	    pix_.push_back(tmp_pix);
        }
      }
    }
  }

  pix_.resize(pix_.size());

  sort(pix_.begin(), pix_.end(), Pixel::LocalOrder);
  mean_intensity_ = 0.0;
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;

  area_ = 0.0;
  total_intensity_ = 0.0;
  total_points_ = 0;
  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
    area_ += iter->Area()*iter->Weight();
    total_intensity_ += iter->Intensity();
    total_points_ += iter->NPoints();
  }
}

void ScalarMap::InitializeFromScalarPixels(ScalarVector& pix,
					   ScalarMapType scalar_map_type) {
  resolution_ = pix[0].Resolution();
  map_type_ = scalar_map_type;

  if (!pix_.empty()) pix_.clear();
  pix_.reserve(pix.size());

  area_ = 0.0;
  total_intensity_ = 0.0;
  total_points_ = 0;
  for (ScalarIterator iter=pix.begin();iter!=pix.end();++iter) {
    if (iter->Resolution() != resolution_) {
      std::cout << "Stomp::ScalarMap::InitializeFromScalarPixels - " <<
	"Incompatible resolutions in ScalarPixel list.  Exiting.\n";
      exit(2);
    }
    area_ += iter->Area()*iter->Weight();
    total_intensity_ += iter->Intensity();
    total_points_ += iter->NPoints();
    pix_.push_back(*iter);
  }

  sort(pix_.begin(), pix_.end(), Pixel::LocalOrder);
  mean_intensity_ = 0.0;
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  use_local_mean_intensity_ = false;
}

bool ScalarMap::AddToMap(AngularCoordinate& ang, double object_weight) {
  ScalarPixel tmp_pix(ang, resolution_, object_weight, 1);
  bool added_point = false;

  ScalarPair iter;
  iter = equal_range(pix_.begin(), pix_.end(), tmp_pix, Pixel::LocalOrder);
  if (iter.first != iter.second) {
    if (map_type_ == ScalarField) {
      iter.first->AddToIntensity(object_weight, 0);
    } else {
      iter.first->AddToIntensity(object_weight, 1);
    }
    total_intensity_ += object_weight;
    total_points_++;
    added_point = true;
  }

  return added_point;
}

bool ScalarMap::AddToMap(WeightedAngularCoordinate& ang) {
  ScalarPixel tmp_pix(ang, resolution_, ang.Weight());
  bool added_point = false;

  ScalarPair iter;
  iter = equal_range(pix_.begin(),pix_.end(),tmp_pix, Pixel::LocalOrder);
  if (iter.first != iter.second) {
    if (map_type_ == ScalarField) {
      iter.first->AddToIntensity(ang.Weight(), 0);
    } else {
      iter.first->AddToIntensity(ang.Weight(), 1);
    }
    total_intensity_ += ang.Weight();
    total_points_++;
    added_point = true;
  }

  return added_point;
}

bool ScalarMap::AddToMap(Pixel& pix) {
  bool added_pixel = false;

  if ((pix.Resolution() <= resolution_) && (map_type_ == ScalarField)) {
    PixelVector pixVec;
    pix.SubPix(resolution_, pixVec);
    for (PixelIterator pix_iter=pixVec.begin();
	 pix_iter!=pixVec.end();++pix_iter) {

      ScalarPixel tmp_pix(pix_iter->PixelX(), pix_iter->PixelY(),
			  pix_iter->Resolution());

      ScalarPair iter;
      iter = equal_range(pix_.begin(), pix_.end(), tmp_pix, Pixel::LocalOrder);
      if (iter.first != iter.second) {
	iter.first->SetIntensity(pix.Weight());
	total_intensity_ += pix.Weight();
      }
    }

    added_pixel = true;
  }

  return added_pixel;
}

void ScalarMap::Coverage(PixelVector& superpix, uint32_t resolution,
			 bool calculate_fraction) {
  if (!superpix.empty()) superpix.clear();

  if (resolution > resolution_) {
    std::cout << "Stomp::ScalarMap::Coverage - " <<
      "WARNING: Requested resolution is higher than " <<
      "the map resolution!\nReseting to map resolution...\n";
    resolution = resolution_;
  }

  std::map<uint32_t, bool> coverage_map;
  for (ScalarIterator iter=pix_.begin(); iter != pix_.end(); ++iter) {
    coverage_map[iter->SuperPix(resolution)] = true;
  }

  for (std::map<uint32_t, bool>::iterator iter=coverage_map.begin();
       iter != coverage_map.end(); ++iter) {
    Pixel tmp_pix(resolution, iter->first, 1.0);

    if (calculate_fraction) {
      tmp_pix.SetWeight(FindUnmaskedFraction(tmp_pix));
    }
    superpix.push_back(tmp_pix);
  }

  sort(superpix.begin(), superpix.end(), Pixel::SuperPixelBasedOrder);
}

bool ScalarMap::Covering(Map& stomp_map, uint32_t maximum_pixels) {
  if (!stomp_map.Empty()) stomp_map.Clear();

  PixelVector pix;
  Coverage(pix, HPixResolution, false);

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
    // equivalent of our current map.
    Coverage(pix, resolution_);
    stomp_map.Initialize(pix);

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

int8_t ScalarMap::FindUnmaskedStatus(Pixel& pix) {
  // By default, there is no overlap between the map and the input pixel.
  int8_t unmasked_status = 0;

  // We have three cases: the input pixel is at lower, equal or higher
  // resolution.  The last two cases can be handled in the same code
  // block since we're looking for matching pixels in our current set.  Note
  // that the higher resolution case output doesn't exactly mean the same
  // thing as it does in the Map class since the pixels in ScalarMaps all have
  // an unmasked fraction which isn't necessarily 1.
  if (pix.Resolution() >= resolution_) {
    Pixel tmp_pix = pix;
    tmp_pix.SetToSuperPix(resolution_);

    ScalarPair iter = equal_range(pix_.begin(), pix_.end(), tmp_pix,
                                  Pixel::LocalOrder);
    if (iter.first != iter.second) unmasked_status = 1;
  }

  // If the input pixel is larger than the ScalarMap pixels, then we scan
  // through the input pixel's sub-pixels, looking for matches.
  if (pix.Resolution() < resolution_) {
    PixelVector pixVec;
    pix.SubPix(resolution_, pixVec);

    for (uint32_t i=0; i < pixVec.size(); i++) {
      ScalarPair iter;
      iter = equal_range(pix_.begin(), pix_.end(), pixVec[i],
                         Pixel::LocalOrder);
      if (iter.first != iter.second) {
        unmasked_status = -1;
        break;
      }
    }
  }

  return unmasked_status;
}

double ScalarMap::FindUnmaskedFraction(Pixel& pix) {
  ScalarPixel scalar_pix(pix.PixelX(), pix.PixelY(), pix.Resolution());

  Resample(scalar_pix);

  return scalar_pix.Weight();
}

void ScalarMap::Resample(ScalarPixel& pix) {
  // Default values if we don't have the input pixel in our ScalarMap.
  double unmasked_fraction = 0.0;
  double total_intensity = 0.0;
  double weighted_intensity = 0.0;
  uint32_t total_points = 0;

  if (pix.Resolution() > resolution_) {
    unmasked_fraction = -1.0;
    total_intensity = -1.0;
    weighted_intensity = -1.0;
    total_points = 0;
  } else {
    if (pix.Resolution() == resolution_) {
      ScalarPixel tmp_pix(pix.PixelX(), pix.PixelY(), pix.Resolution());

      ScalarPair iter = equal_range(pix_.begin(), pix_.end(), tmp_pix,
                                    Pixel::LocalOrder);
      if (iter.first != iter.second) {
        unmasked_fraction = iter.first->Weight();
        total_intensity = iter.first->Intensity();
        weighted_intensity = total_intensity*unmasked_fraction;
        total_points = iter.first->NPoints();
      }
    } else {
      double pixel_fraction =
	  1.0*pix.Resolution()*pix.Resolution()/(resolution_*resolution_);

      PixelVector pixVec;
      pix.SubPix(resolution_, pixVec);

      for (uint32_t i=0; i < pixVec.size(); i++) {
        ScalarPair iter;
        iter = equal_range(pix_.begin(), pix_.end(), pixVec[i],
                           Pixel::LocalOrder);
        if (iter.first != iter.second) {
          unmasked_fraction += pixel_fraction*iter.first->Weight();
          total_intensity += iter.first->Intensity();
          weighted_intensity +=
              iter.first->Intensity()*pixel_fraction*iter.first->Weight();
          total_points += iter.first->NPoints();
        }
      }
    }
  }

  if (unmasked_fraction > 0.0000001) {
    // If we normalize weighted_intensity by the unmasked fraction, we have an
    // area-averaged value of the intensity over the pixel.
    weighted_intensity /= unmasked_fraction;

    // If our pixels were encoding the over-density of the scalar fields, then
    // we need to convert those values back into the raw values.
    if (converted_to_overdensity_) {
      if (map_type_ == DensityField) {
	weighted_intensity =
            weighted_intensity*mean_intensity_ + mean_intensity_;
      } else {
	weighted_intensity += mean_intensity_;
      }
      if (map_type_ == DensityField) {
	// For a DensityField, we want to return the total intensity for all of
	// the area subtended by the input pixel.  If we were dealing with a raw
	// map, then that value is already in total intensity.  Since we've
	// got an over-density map, then we need to convert the weighted
	// average into a total by incorporating the unmasked area.
	total_intensity = weighted_intensity*unmasked_fraction*pix.Area();
      }
      if (map_type_ == SampledField) {
	// Likewise for a SampledField, but the normalization is based on the
	// total number of points in the pixel.
	total_intensity = weighted_intensity*total_points;
      }
    }

    if (map_type_ == ScalarField) {
      // For a ScalarField, we want the average value over the area indicated
      // by the input pixel.
      total_intensity = weighted_intensity;
    }
  } else {
    unmasked_fraction = 0.0;
    total_intensity = 0.0;
    total_points = 0;
  }

  pix.SetWeight(unmasked_fraction);
  pix.SetIntensity(total_intensity);
  pix.SetNPoints(total_points);
}

double ScalarMap::FindIntensity(Pixel& pix) {
  ScalarPixel scalar_pix(pix.PixelX(), pix.PixelY(), pix.Resolution());

  Resample(scalar_pix);

  return scalar_pix.Intensity();
}

double ScalarMap::FindDensity(Pixel& pix) {
  ScalarPixel scalar_pix(pix.PixelX(), pix.PixelY(), pix.Resolution());

  Resample(scalar_pix);

  return scalar_pix.Intensity()/(scalar_pix.Weight()*scalar_pix.Area());
}

double ScalarMap::FindPointDensity(Pixel& pix) {
  ScalarPixel scalar_pix(pix.PixelX(), pix.PixelY(), pix.Resolution());

  Resample(scalar_pix);

  return 1.0*scalar_pix.NPoints()/(scalar_pix.Weight()*scalar_pix.Area());
}

double ScalarMap::FindLocalArea(AngularCoordinate& ang,
				double theta_max, double theta_min) {
  Pixel center_pix(ang,resolution_,0.0);

  PixelVector pixVec;
  if (theta_min < 0.0) {
    center_pix.WithinRadius(theta_max, pixVec);
  } else {
    center_pix.WithinAnnulus(theta_min, theta_max, pixVec);
  }

  double total_area = 0.0;
  for (PixelIterator iter=pixVec.begin();iter!=pixVec.end();++iter) {
    total_area += FindUnmaskedFraction(*iter);
  }

  if (total_area > 0.0000001) {
    return total_area*center_pix.Area();
  } else {
    return 0.0;
  }
}

double ScalarMap::FindLocalAverageIntensity(AngularCoordinate& ang,
				   double theta_max, double theta_min) {
  Pixel center_pix(ang,resolution_,0.0);

  PixelVector pixVec;
  if (theta_min < 0.0) {
    center_pix.WithinRadius(theta_max, pixVec);
  } else {
    center_pix.WithinAnnulus(theta_min, theta_max, pixVec);
  }

  double total_area = 0.0;
  double total_intensity = 0.0;
  for (PixelIterator iter=pixVec.begin();iter!=pixVec.end();++iter) {
  	double unmasked = FindUnmaskedFraction(*iter);
    total_area += unmasked;
    total_intensity += unmasked*FindIntensity(*iter);
  }

  if (total_area > 0.0000001) {
    return total_intensity/(total_area*center_pix.Area());
  } else {
    return 0.0;
  }
}

double ScalarMap::FindLocalIntensity(AngularCoordinate& ang,
				     double theta_max, double theta_min) {
  Pixel center_pix(ang,resolution_,0.0);

  PixelVector pixVec;
  if (theta_min < 0.0) {
    center_pix.WithinRadius(theta_max, pixVec);
  } else {
    center_pix.WithinAnnulus(theta_min, theta_max, pixVec);
  }

  double total_intensity = 0.0;
  for (PixelIterator iter=pixVec.begin();iter!=pixVec.end();++iter) {
    total_intensity += FindIntensity(*iter);
  }

  return total_intensity;
}

double ScalarMap::FindLocalDensity(AngularCoordinate& ang,
				   double theta_max, double theta_min) {
  Pixel center_pix(ang,resolution_,0.0);

  PixelVector pixVec;
  if (theta_min < 0.0) {
    center_pix.WithinRadius(theta_max, pixVec);
  } else {
    center_pix.WithinAnnulus(theta_min, theta_max, pixVec);
  }

  double total_area = 0.0;
  double total_intensity = 0.0;
  for (PixelIterator iter=pixVec.begin();iter!=pixVec.end();++iter) {
    total_area += FindUnmaskedFraction(*iter);
    total_intensity += FindIntensity(*iter);
  }

  if (total_area > 0.0000001) {
    return total_intensity/(total_area*center_pix.Area());
  } else {
    return 0.0;
  }
}

double ScalarMap::FindLocalPointDensity(AngularCoordinate& ang,
					double theta_max, double theta_min) {
  Pixel center_pix(ang,resolution_,0.0);

  PixelVector pixVec;
  if (theta_min < 0.0) {
    center_pix.WithinRadius(theta_max, pixVec);
  } else {
    center_pix.WithinAnnulus(theta_min, theta_max, pixVec);
  }

  double total_area = 0.0;
  double total_point_density = 0.0;
  for (PixelIterator iter=pixVec.begin();iter!=pixVec.end();++iter) {
    double area = FindUnmaskedFraction(*iter);
    if (area > 0.0000001) {
      total_area += area;
      total_point_density += FindPointDensity(*iter)*area*center_pix.Area();
    }
  }

  if (total_area > 0.0000001) {
    return total_point_density/(total_area*center_pix.Area());
  } else {
    return 0.0;
  }
}

void ScalarMap::CalculateMeanIntensity() {
  if (!use_local_mean_intensity_) {
    double sum_pixel = 0.0;
    mean_intensity_ = 0.0;

    switch (map_type_) {
    case ScalarField:
      for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	mean_intensity_ += iter->Intensity()*iter->Weight();
	sum_pixel += iter->Weight();
      }
      break;
    case DensityField:
      for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	mean_intensity_ += iter->Intensity()/(iter->Area()*iter->Weight());
	sum_pixel += 1.0;
      }
    break;
    case SampledField:
      for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	mean_intensity_ += iter->MeanIntensity()*iter->Weight();
	sum_pixel += iter->Weight();
      }
      break;
    }

    mean_intensity_ /= sum_pixel;
    calculated_mean_intensity_ = true;
  } else {
    if (RegionsInitialized()) {
      local_mean_intensity_.clear();
      std::vector<double> sum_pixel_vec;

      double sum_pixel = 0.0;
      mean_intensity_ = 0.0;

      local_mean_intensity_.reserve(NRegion());
      sum_pixel_vec.reserve(NRegion());

      for (uint16_t i=0;i<NRegion();i++) {
	sum_pixel_vec[i] = 0.0;
	local_mean_intensity_[i] = 0.0;
      }

      switch (map_type_) {
      case ScalarField:
	for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	  mean_intensity_ += iter->Intensity()*iter->Weight();
	  sum_pixel += iter->Weight();

	  int16_t region_idx = FindRegion(*iter);
	  if (region_idx != -1 && region_idx < NRegion()) {
	    local_mean_intensity_[region_idx] +=
	      iter->Intensity()*iter->Weight();
	    sum_pixel_vec[region_idx] += iter->Weight();
	  }
	}
	break;
      case DensityField:
	for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	  mean_intensity_ += iter->Intensity()/(iter->Area()*iter->Weight());
	  sum_pixel += 1.0;

	  int16_t region_idx = FindRegion(*iter);
	  if (region_idx != -1 && region_idx < NRegion()) {
	    local_mean_intensity_[region_idx] +=
	      iter->Intensity()/(iter->Area()*iter->Weight());
	    sum_pixel_vec[region_idx] += 1.0;
	  }
	}
	break;
      case SampledField:
	for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	  mean_intensity_ += iter->MeanIntensity()*iter->Weight();
	  sum_pixel += iter->Weight();

	  int16_t region_idx = FindRegion(*iter);
	  if (region_idx != -1 && region_idx < NRegion()) {
	    local_mean_intensity_[region_idx] +=
	      iter->MeanIntensity()*iter->Weight();
	    sum_pixel_vec[region_idx] += iter->Weight();
	  }
	}
	break;
      }

      for (uint16_t i=0;i<NRegion();i++) {
	local_mean_intensity_[i] /= sum_pixel_vec[i];
      }
      mean_intensity_ /= sum_pixel;

      calculated_mean_intensity_ = true;
    }
  }
}

void ScalarMap::ConvertToOverDensity() {
  if (!calculated_mean_intensity_) CalculateMeanIntensity();

  if (!use_local_mean_intensity_) {
    // Only do this conversion if we've got raw intensity values in our map.
    if (!converted_to_overdensity_) {
      if (map_type_ == DensityField) {
	for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter)
	  iter->ConvertToFractionalOverDensity(mean_intensity_);
      } else {
	for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter)
	  iter->ConvertToOverDensity(mean_intensity_);
      }
    }

    converted_to_overdensity_ = true;
  } else {
    if (RegionsInitialized()) {
      // Only do this conversion if we've got raw intensity values in our map.
      if (!converted_to_overdensity_) {
	if (map_type_ == DensityField) {
	  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	    int16_t region_idx = FindRegion(*iter);
	    if (region_idx != -1 && region_idx < NRegion()) {
	      double local_mean_intensity = local_mean_intensity_[region_idx];
	      iter->ConvertToFractionalOverDensity(local_mean_intensity);
	    } else {
	      std::cout <<
		"Stomp::ScalarMap::ConvertToOverDensity - " <<
		"Failed to find region when calculating over-density!\n";
	      exit(2);
	    }
	  }
	} else {
	  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	    int16_t region_idx = FindRegion(*iter);
	    if (region_idx != -1 && region_idx < NRegion()) {
	      double local_mean_intensity = local_mean_intensity_[region_idx];
	      iter->ConvertToOverDensity(local_mean_intensity);
	    } else {
	      std::cout <<
		"Stomp::ScalarMap::ConvertToOverDensity - " <<
		"Failed to find region when calculating over-density!\n";
	      exit(2);
	    }
	  }
	}
      }

      converted_to_overdensity_ = true;
    }
  }
}

void ScalarMap::ConvertFromOverDensity() {
  if (!calculated_mean_intensity_) CalculateMeanIntensity();

  if (!use_local_mean_intensity_) {
    // Only do this conversion if we've got over-density values in our map.
    if (converted_to_overdensity_) {
      if (map_type_ == DensityField) {
	for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter)
	  iter->ConvertFromFractionalOverDensity(mean_intensity_);
      } else {
	for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter)
	  iter->ConvertFromOverDensity(mean_intensity_);
      }
    }

    converted_to_overdensity_ = false;
  } else {
    if (RegionsInitialized()) {
      // Only do this conversion if we've got over-density values in our map.
      if (converted_to_overdensity_) {
	if (map_type_ == DensityField) {
	  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	    int16_t region_idx = FindRegion(*iter);
	    if (region_idx != -1 && region_idx < NRegion()) {
	      double local_mean_intensity = local_mean_intensity_[region_idx];
	      iter->ConvertFromFractionalOverDensity(local_mean_intensity);
	    } else {
	      std::cout <<
		"Stomp::ScalarMap::ConvertFromOverDensity - " <<
		"Failed to find region when calculating over-density!\n";
	      exit(2);
	    }
	  }
	} else {
	  for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	    int16_t region_idx = FindRegion(*iter);
	    if (region_idx != -1 && region_idx < NRegion()) {
	      double local_mean_intensity = local_mean_intensity_[region_idx];
	      iter->ConvertFromOverDensity(local_mean_intensity);
	    } else {
	      std::cout <<
		"Stomp::ScalarMap::ConvertFromOverDensity - " <<
		"Failed to find region when calculating over-density!\n";
	      exit(2);
	    }
	  }
	}
      }

      converted_to_overdensity_ = false;
    }
  }
}

bool ScalarMap::UseLocalMeanIntensity(bool use_local_mean) {
  if (RegionsInitialized()) {
    // Before we potentially change the way that the mean intensity is
    // calculated, we need to make sure that we're working with the raw
    // intensities instead of over-densities.  The latter are a function of
    // the mean intensity, so we need to make sure that we've managed the
    // transition in such a way that we can get back to the raw intensities.
    bool overdensity_map = IsOverDensityMap();

    if (overdensity_map) ConvertFromOverDensity();

    use_local_mean_intensity_ = use_local_mean;
    calculated_mean_intensity_ = false;

    if (overdensity_map) ConvertToOverDensity();

    return true;
  } else {
    return false;
  }
}

bool ScalarMap::UsingLocalMeanIntensity() {
  return use_local_mean_intensity_;
}

bool ScalarMap::ImprintMap(Map& stomp_map, bool use_mean_local_intensity) {
  if (!use_mean_local_intensity) {
    PixelVector pixVec;

    pixVec.reserve(pix_.size());

    for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter)
      pixVec.push_back(Pixel(iter->PixelX(), iter->PixelY(),
			     iter->Resolution(), iter->Intensity()));

    return stomp_map.ImprintMap(pixVec);
  } else {
    if (RegionsInitialized()) {
      if (!UsingLocalMeanIntensity()) UseLocalMeanIntensity(true);

      PixelVector pixVec;

      pixVec.reserve(pix_.size());

      for (ScalarIterator iter=pix_.begin();iter!=pix_.end();++iter) {
	int16_t region_idx = FindRegion(*iter);
	if (region_idx != -1 && region_idx < NRegion()) {
	  double local_mean_intensity = local_mean_intensity_[region_idx];
	  pixVec.push_back(Pixel(iter->PixelX(), iter->PixelY(),
				 iter->Resolution(), local_mean_intensity));
      	}
      }

      return stomp_map.ImprintMap(pixVec);
    } else {
      return false;
    }
  }
}

void ScalarMap::AutoCorrelate(AngularCorrelation& wtheta) {
  ThetaIterator theta_begin = wtheta.Begin(resolution_);
  ThetaIterator theta_end = wtheta.End(resolution_);

  if (theta_begin != theta_end) {
    bool convert_back_to_raw = false;
    if (!converted_to_overdensity_) {
      ConvertToOverDensity();
      convert_back_to_raw = true;
    }

    for (ThetaIterator theta_iter=theta_begin;
	 theta_iter!=theta_end;++theta_iter) AutoCorrelate(theta_iter);

    if (convert_back_to_raw) ConvertFromOverDensity();
  } else {
    std::cout << "Stomp::ScalarMap::AutoCorrelate - " <<
      "No angular bins have resolution " << resolution_ << "...\n";
  }
}

void ScalarMap::AutoCorrelate(ThetaIterator theta_iter) {
  bool convert_back_to_raw = false;
  if (!converted_to_overdensity_) {
    ConvertToOverDensity();
    convert_back_to_raw = true;
  }

  theta_iter->ResetPixelWtheta();

  for (ScalarIterator map_iter=pix_.begin();
       map_iter!=pix_.end();++map_iter) {
    // From a readability standpoint, the simplest thing to do is to
    // use the Pixel.WithinAnnulus method to get a vector of pixel
    // candidates and then check each in turn.  This isn't the fastest thing
    // we could do (probably), but it is easy to read.
    ScalarVector pixVec;
    map_iter->_WithinAnnulus(*theta_iter, pixVec);

    for (ScalarIterator pix_iter=pixVec.begin();
         pix_iter!=pixVec.end();++pix_iter) {
      // We can skip all pixels that preceed our current pixel since the
      // products would be the same.  This keeps us from double counting and
      // saves some time.
      if (Pixel::LocalOrder(*map_iter, *pix_iter)) {
        ScalarPair iter = equal_range(map_iter, pix_.end(), *pix_iter,
                                      Pixel::LocalOrder);
        if (iter.first != iter.second) {
          theta_iter->AddToPixelWtheta(map_iter->Intensity()*map_iter->Weight()*
                                       iter.first->Intensity()*
                                       iter.first->Weight(),
                                       map_iter->Weight()*iter.first->Weight());
        }
      }
    }
  }

  if (convert_back_to_raw) ConvertFromOverDensity();
}

void ScalarMap::AutoCorrelateWithRegions(AngularCorrelation& wtheta) {
  ThetaIterator theta_begin = wtheta.Begin(resolution_);
  ThetaIterator theta_end = wtheta.End(resolution_);

  if (theta_begin != theta_end) {
    bool convert_back_to_raw = false;
    if (!converted_to_overdensity_) {
      ConvertToOverDensity();
      convert_back_to_raw = true;
    }

    for (ThetaIterator theta_iter=theta_begin;
	 theta_iter!=theta_end;++theta_iter)
      AutoCorrelateWithRegions(theta_iter);

    if (convert_back_to_raw) ConvertFromOverDensity();
  } else {
    std::cout << "Stomp::ScalarMap::AutoCorrelateWithRegions - " <<
      "No angular bins have resolution " << resolution_ << "...\n";
  }
}

void ScalarMap::AutoCorrelateWithRegions(ThetaIterator theta_iter) {
  bool convert_back_to_raw = false;
  if (!converted_to_overdensity_) {
    ConvertToOverDensity();
    convert_back_to_raw = true;
  }

  if (theta_iter->NRegion() != NRegion()) {
    theta_iter->ClearRegions();
    theta_iter->InitializeRegions(NRegion());
  }

  theta_iter->ResetPixelWtheta();

  // Same as the method without the regions, we just need to keep track of
  // the region values for the two pixels involved.
  for (ScalarIterator map_iter=pix_.begin();
       map_iter!=pix_.end();++map_iter) {
    uint32_t map_region = Region(map_iter->SuperPix(RegionResolution()));

    ScalarVector pixVec;
    map_iter->_WithinAnnulus(*theta_iter, pixVec);

    for (ScalarIterator pix_iter=pixVec.begin();
         pix_iter!=pixVec.end();++pix_iter) {
      if (Pixel::LocalOrder(*map_iter, *pix_iter)) {
        ScalarPair iter = equal_range(map_iter, pix_.end(), *pix_iter,
                                      Pixel::LocalOrder);
        if (iter.first != iter.second) {
          uint32_t pix_region =
              Region(iter.first->SuperPix(RegionResolution()));
          theta_iter->AddToPixelWtheta(map_iter->Intensity()*map_iter->Weight()*
                                       iter.first->Intensity()*
                                       iter.first->Weight(),
                                       map_iter->Weight()*iter.first->Weight(),
                                       map_region, pix_region);
        }
      }
    }
  }

  if (convert_back_to_raw) ConvertFromOverDensity();
}

void ScalarMap::CrossCorrelate(ScalarMap& scalar_map,
			       AngularCorrelation& wtheta) {
  ThetaIterator theta_begin = wtheta.Begin(resolution_);
  ThetaIterator theta_end = wtheta.End(resolution_);

  if (theta_begin != theta_end) {
    bool convert_back_to_raw = false;
    if (!converted_to_overdensity_) {
      ConvertToOverDensity();
      convert_back_to_raw = true;
    }

    bool convert_input_map_back_to_raw = false;
    if (!scalar_map.IsOverDensityMap()) {
      scalar_map.ConvertToOverDensity();
      convert_input_map_back_to_raw = true;
    }

    for (ThetaIterator theta_iter=theta_begin;
	 theta_iter!=theta_end;++theta_iter)
      CrossCorrelate(scalar_map,theta_iter);

    if (convert_back_to_raw) ConvertFromOverDensity();
    if (convert_input_map_back_to_raw) scalar_map.ConvertFromOverDensity();
  } else {
    std::cout << "Stomp::ScalarMap::CrossCorrelate - " <<
      "No angular bins have resolution " << resolution_ << "...\n";
  }
}

void ScalarMap::CrossCorrelate(ScalarMap& scalar_map,
			       ThetaIterator theta_iter) {
  if (resolution_ != scalar_map.Resolution()) {
    std::cout << "Stomp::ScalarMap::CrossCorrelate - " <<
      "Map resolutions must match!  Exiting...\n";
    exit(1);
  }

  bool convert_back_to_raw = false;
  if (!converted_to_overdensity_) {
    ConvertToOverDensity();
    convert_back_to_raw = true;
  }

  bool convert_input_map_back_to_raw = false;
  if (!scalar_map.IsOverDensityMap()) {
    scalar_map.ConvertToOverDensity();
    convert_input_map_back_to_raw = true;
  }

  theta_iter->ResetPixelWtheta();

  // Same as the auto-correlation case, although this time we're iterating
  // over the pixels in the input ScalarMap.  We also don't exclude pixels
  // if their LocalOrder puts them before our input ScalarMap's pixel since
  // double-counting isn't an issue.
  for (ScalarIterator map_iter=scalar_map.Begin();
       map_iter!=scalar_map.End();++map_iter) {
    ScalarVector pixVec;
    map_iter->_WithinAnnulus(*theta_iter, pixVec);

    for (ScalarIterator pix_iter=pixVec.begin();
         pix_iter!=pixVec.end();++pix_iter) {
      ScalarPair iter = equal_range(pix_.begin(), pix_.end(), *pix_iter,
                                    Pixel::LocalOrder);
      if (iter.first != iter.second) {
        theta_iter->AddToPixelWtheta(map_iter->Intensity()*map_iter->Weight()*
                                     iter.first->Intensity()*
                                     iter.first->Weight(),
                                     map_iter->Weight()*iter.first->Weight());
      }
    }
  }

  if (convert_back_to_raw) ConvertFromOverDensity();
  if (convert_input_map_back_to_raw) scalar_map.ConvertFromOverDensity();
}

void ScalarMap::CrossCorrelateWithRegions(ScalarMap& scalar_map,
					   AngularCorrelation& wtheta) {
  ThetaIterator theta_begin = wtheta.Begin(resolution_);
  ThetaIterator theta_end = wtheta.End(resolution_);

  if (theta_begin != theta_end) {
    bool convert_back_to_raw = false;
    if (!converted_to_overdensity_) {
      ConvertToOverDensity();
      convert_back_to_raw = true;
    }

    bool convert_input_map_back_to_raw = false;
    if (!scalar_map.IsOverDensityMap()) {
      scalar_map.ConvertToOverDensity();
      convert_input_map_back_to_raw = true;
    }

    for (ThetaIterator theta_iter=theta_begin;
	 theta_iter!=theta_end;++theta_iter)
      CrossCorrelateWithRegions(scalar_map,theta_iter);

    if (convert_back_to_raw) ConvertFromOverDensity();
    if (convert_input_map_back_to_raw) scalar_map.ConvertFromOverDensity();
  } else {
    std::cout << "Stomp::ScalarMap::CrossCorrelateWithRegions - " <<
      "No angular bins have resolution " << resolution_ << "...\n";
  }
}

void ScalarMap::CrossCorrelateWithRegions(ScalarMap& scalar_map,
					  ThetaIterator theta_iter) {
  if (resolution_ != scalar_map.Resolution()) {
    std::cout << "Stomp::ScalarMap::CrossCorrelateWithRegions - " <<
      "Map resolutions must match!  Exiting...\n";
    exit(1);
  }

  if (NRegion() != scalar_map.NRegion()) {
    std::cout << "Stomp::ScalarMap::CrossCorrelateWithRegions - " <<
      "Map regionation must match!  Exiting...\n";
    exit(1);
  }

  bool convert_back_to_raw = false;
  if (!converted_to_overdensity_) {
    ConvertToOverDensity();
    convert_back_to_raw = true;
  }

  bool convert_input_map_back_to_raw = false;
  if (!scalar_map.IsOverDensityMap()) {
    scalar_map.ConvertToOverDensity();
    convert_input_map_back_to_raw = true;
  }

  if (theta_iter->NRegion() != NRegion()) {
    theta_iter->ClearRegions();
    theta_iter->InitializeRegions(NRegion());
  }

  theta_iter->ResetPixelWtheta();

  // Same as the method without the regions, we just need to keep track of
  // the region values for the two pixels involved.
  for (ScalarIterator map_iter=scalar_map.Begin();
       map_iter!=scalar_map.End();++map_iter) {
    int16_t map_region = Region(map_iter->SuperPix(RegionResolution()));

    ScalarVector pixVec;
    map_iter->_WithinAnnulus(*theta_iter, pixVec);

    for (ScalarIterator pix_iter=pixVec.begin();
         pix_iter!=pixVec.end();++pix_iter) {
      ScalarPair iter = equal_range(pix_.begin(), pix_.end(), *pix_iter,
                                    Pixel::LocalOrder);
      if (iter.first != iter.second) {
        int16_t pix_region =
            Region(iter.first->SuperPix(RegionResolution()));
        theta_iter->AddToPixelWtheta(map_iter->Intensity()*map_iter->Weight()*
                                     iter.first->Intensity()*
                                     iter.first->Weight(),
                                     map_iter->Weight()*iter.first->Weight(),
                                     map_region, pix_region);
      }
    }
  }

  if (convert_back_to_raw) ConvertFromOverDensity();
  if (convert_input_map_back_to_raw) scalar_map.ConvertFromOverDensity();
}

double ScalarMap::Variance() {
  bool convert_back_to_raw = false;
  if (!converted_to_overdensity_) {
    ConvertToOverDensity();
    convert_back_to_raw = true;
  }

  double variance = 0.0, variance_norm = 0.0;

  for (ScalarIterator map_iter=pix_.begin();
       map_iter!=pix_.end();++map_iter) {
    variance +=
      map_iter->Intensity()*map_iter->Weight()*
      map_iter->Intensity()*map_iter->Weight();
    variance_norm += map_iter->Weight()*map_iter->Weight();
  }

  if (convert_back_to_raw) ConvertFromOverDensity();

  return (DoubleGT(variance_norm, 1.0e-10) ? variance/variance_norm : 0.0);
}

double ScalarMap::Covariance(ScalarMap& scalar_map) {
  if (resolution_ != scalar_map.Resolution()) {
    std::cout << "Stomp::ScalarMap::Covariance - " <<
      "Map resolutions must match!  Exiting...\n";
    exit(1);
  }

  bool convert_back_to_raw = false;
  if (!converted_to_overdensity_) {
    ConvertToOverDensity();
    convert_back_to_raw = true;
  }

  bool convert_input_map_back_to_raw = false;
  if (!scalar_map.IsOverDensityMap()) {
    scalar_map.ConvertToOverDensity();
    convert_input_map_back_to_raw = true;
  }

  ScalarIterator search_begin = pix_.begin();
  double covariance = 0.0, covariance_norm = 0.0;

  for (ScalarIterator map_iter=scalar_map.Begin();
       map_iter!=scalar_map.End();++map_iter) {
    ScalarPair iter = equal_range(search_begin, pix_.end(), *map_iter,
				  Pixel::LocalOrder);
    if (iter.first != iter.second) {
      covariance +=
	iter.first->Intensity()*iter.first->Weight()*
	map_iter->Intensity()*map_iter->Weight();
      covariance_norm += iter.first->Weight()*map_iter->Weight();
    }
    search_begin = iter.second;
  }

  if (convert_back_to_raw) ConvertFromOverDensity();
  if (convert_input_map_back_to_raw) scalar_map.ConvertFromOverDensity();

  return (DoubleGT(covariance_norm, 1.0e-10) ?
	  covariance/covariance_norm : 0.0);
}

void ScalarMap::VarianceWithErrors(double& variance, double& variance_error) {
  bool convert_back_to_raw = false;
  if (!converted_to_overdensity_) {
    ConvertToOverDensity();
    convert_back_to_raw = true;
  }

  std::vector<double> region_variance, region_variance_norm;
  region_variance.reserve(NRegion());
  region_variance_norm.reserve(NRegion());

  for (uint16_t i=0;i<NRegion();i++) {
    region_variance.push_back(0.0);
    region_variance_norm.push_back(0.0);
  }

  for (ScalarIterator map_iter=pix_.begin();
       map_iter!=pix_.end();++map_iter) {
    int16_t pix_region = Region(map_iter->SuperPix(RegionResolution()));
    region_variance[pix_region] +=
      map_iter->Intensity()*map_iter->Weight()*
      map_iter->Intensity()*map_iter->Weight();
    region_variance_norm[pix_region] += map_iter->Weight()*map_iter->Weight();
  }

  variance = 0.0;
  variance_error = 0.0;
  uint16_t n_region = 0;
  for (uint16_t i=0;i<NRegion();i++) {
    if (DoubleGT(region_variance_norm[i], 1.0e-10)) {
      variance += region_variance[i]/region_variance_norm[i];
      n_region++;
    }
  }

  if (n_region > 0) {
    variance /= n_region;
    for (uint16_t i=0;i<NRegion();i++) {
      if (DoubleGT(region_variance_norm[i], 1.0e-10)) {
	variance_error +=
	  (variance - region_variance[i]/region_variance_norm[i])*
	  (variance - region_variance[i]/region_variance_norm[i]);
      }
    }
    variance_error = sqrt(variance_error)/n_region;
  }

  if (convert_back_to_raw) ConvertFromOverDensity();
}

void ScalarMap::CovarianceWithErrors(ScalarMap& scalar_map, double& covariance,
				     double& covariance_error) {
  if (resolution_ != scalar_map.Resolution()) {
    std::cout << "Stomp::ScalarMap::CovarianceWithErrors - " <<
      "Map resolutions must match!  Exiting...\n";
    exit(1);
  }

  if (NRegion() != scalar_map.NRegion()) {
    std::cout << "Stomp::ScalarMap::CovarianceWithErrors - " <<
      "Map regionation must match!  Exiting...\n";
    exit(1);
  }

  bool convert_back_to_raw = false;
  if (!converted_to_overdensity_) {
    ConvertToOverDensity();
    convert_back_to_raw = true;
  }

  bool convert_input_map_back_to_raw = false;
  if (!scalar_map.IsOverDensityMap()) {
    scalar_map.ConvertToOverDensity();
    convert_input_map_back_to_raw = true;
  }

  std::vector<double> region_covariance, region_covariance_norm;
  region_covariance.reserve(NRegion());
  region_covariance_norm.reserve(NRegion());

  for (uint16_t i=0;i<NRegion();i++) {
    region_covariance.push_back(0.0);
    region_covariance_norm.push_back(0.0);
  }

  ScalarIterator search_begin = pix_.begin();
  for (ScalarIterator map_iter=scalar_map.Begin();
       map_iter!=scalar_map.End();++map_iter) {
    ScalarPair iter = equal_range(search_begin, pix_.end(), *map_iter,
				  Pixel::LocalOrder);
    if (iter.first != iter.second) {
      int16_t pix_region = Region(map_iter->SuperPix(RegionResolution()));
      region_covariance[pix_region] +=
	iter.first->Intensity()*iter.first->Weight()*
	map_iter->Intensity()*map_iter->Weight();
      region_covariance_norm[pix_region] +=
	iter.first->Weight()*map_iter->Weight();
    }
    search_begin = iter.second;
  }

  covariance = 0.0;
  covariance_error = 0.0;
  uint16_t n_region = 0;
  for (uint16_t i=0;i<NRegion();i++) {
    if (DoubleGT(region_covariance_norm[i], 1.0e-10)) {
      covariance += region_covariance[i]/region_covariance_norm[i];
      n_region++;
    }
  }

  if (n_region > 0) {
    covariance /= n_region;
    for (uint16_t i=0;i<NRegion();i++) {
      if (DoubleGT(region_covariance_norm[i], 1.0e-10)) {
	covariance_error +=
	  (covariance - region_covariance[i]/region_covariance_norm[i])*
	  (covariance - region_covariance[i]/region_covariance_norm[i]);
      }
    }
    covariance_error = sqrt(covariance_error)/n_region;
  }

  if (convert_back_to_raw) ConvertFromOverDensity();
  if (convert_input_map_back_to_raw) scalar_map.ConvertFromOverDensity();
}

uint32_t ScalarMap::Resolution() {
  return resolution_;
}

double ScalarMap::Intensity() {
  return total_intensity_;
}

int ScalarMap::NPoints() {
  return total_points_;
}

double ScalarMap::Density() {
  return total_intensity_/area_;
}

double ScalarMap::PointDensity() {
  return 1.0*total_points_/area_;
}

ScalarIterator ScalarMap::Begin() {
  return pix_.begin();
}

ScalarIterator ScalarMap::End() {
  return pix_.end();
}

double ScalarMap::Area() {
  return area_;
}

uint32_t ScalarMap::Size() {
  return pix_.size();
}

uint32_t ScalarMap::MinResolution() {
  return resolution_;
}

uint32_t ScalarMap::MaxResolution() {
  return resolution_;
}

uint8_t ScalarMap::MinLevel() {
  return Pixel::ResolutionToLevel(resolution_);
}

uint8_t ScalarMap::MaxLevel() {
  return Pixel::ResolutionToLevel(resolution_);
}

bool ScalarMap::Empty() {
  return (pix_.empty() ? true : false);
}

void ScalarMap::Clear() {
  area_ = 0.0;
  resolution_ = 0;
  total_points_ = 0;
  mean_intensity_ = 0.0;
  total_intensity_ = 0.0;
  if (!pix_.empty()) pix_.clear();
  converted_to_overdensity_ = false;
  calculated_mean_intensity_ = false;
  ClearRegions();
}

double ScalarMap::MeanIntensity() {
  if (!calculated_mean_intensity_) CalculateMeanIntensity();
  return mean_intensity_;
}

bool ScalarMap::IsOverDensityMap() {
  return converted_to_overdensity_;
}

ScalarMap::ScalarMapType ScalarMap::MapType() {
  return map_type_;
}

} // end namespace Stomp
