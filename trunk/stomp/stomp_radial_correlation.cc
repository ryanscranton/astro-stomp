// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the class for calculating radial correlations on
// the sphere.  In general, different methods are more efficient on small vs.
// large angular scales, so this class draws on nearly the entire breadth of
// the STOMP library.

#include "stomp_core.h"
#include "stomp_angular_correlation.h"
#include "stomp_radial_correlation.h"
#include "stomp_map.h"
#include "stomp_scalar_map.h"
#include "stomp_tree_map.h"

namespace Stomp {

RadialCorrelation::RadialCorrelation(double r_min, double r_max,
				     double bins_per_decade) {
  double unit_double = floor(log10(r_min))*bins_per_decade;
  double r = pow(10.0, unit_double/bins_per_decade);

  while (r < r_max) {
    if (DoubleGE(r, r_min) && (r < r_max)) {
      RadialBin radialbin;
      radialbin.SetRadiusMin(r);
      radialbin.SetRadiusMax(pow(10.0,(unit_double+1.0)/bins_per_decade));
      radialbin.SetRadius(pow(10.0,0.5*(log10(radialbin.RadiusMin())+
				       log10(radialbin.RadiusMax()))));
      radialbin_.push_back(radialbin);
    }
    unit_double += 1.0;
    r = pow(10.0,unit_double/bins_per_decade);
  }

  r_min_ = radialbin_[0].RadiusMin();
  r_max_ = radialbin_[radialbin_.size()-1].RadiusMax();

  regionation_resolution_ = 0;

  min_resolution_ = HPixResolution;
  max_resolution_ = HPixResolution;

  manual_resolution_break_ = false;

  //UseOnlyPairs currently hard coded in.
  UseOnlyPairs();
}

RadialCorrelation::RadialCorrelation(uint32_t n_bins,
				     double r_min, double r_max) {
  double dr = (r_max - r_min)/n_bins;

  for (uint32_t i=0;i<n_bins;i++) {
    RadialBin radialbin;
    radialbin.SetRadiusMin(r_min + i*dr);
    radialbin.SetRadiusMin(r_min + (i+1)*dr);
    radialbin.SetRadius(0.5*(radialbin.RadiusMin()+radialbin.RadiusMax()));
    radialbin_.push_back(radialbin);
  }

  r_min_ = radialbin_[0].RadiusMin();
  r_max_ = radialbin_[n_bins-1].RadiusMax();

  regionation_resolution_ = 0;
    
  min_resolution_ = HPixResolution;
  max_resolution_ = HPixResolution;

  manual_resolution_break_ = false;

  //UseOnlyPairs currently hard coded in.
  UseOnlyPairs();
}

void RadialCorrelation::InitializeRegions(int16_t n_regions) {
  n_region_ = n_regions;
  for (RadialIterator iter=Begin();iter!=End();++iter)
    iter->InitializeRegions(n_region_);
}

void RadialCorrelation::ClearRegions() {
  n_region_ = 0;
  for (RadialIterator iter=Begin();iter!=End();++iter)
    iter->ClearRegions();
  regionation_resolution_ = 0;
}

void RadialCorrelation::FindAutoCorrelation(Map& stomp_map,
					    CosmoVector& galaxy,
					    uint8_t random_iterations) {
  //if (!manual_resolution_break_)
  //  AutoMaxResolution(galaxy.size(), stomp_map.Area());

  //if (theta_pixel_begin_ != theta_pixel_end_)
  //  FindPixelAutoCorrelation(stomp_map, galaxy);

  //if (theta_pair_begin_ != theta_pair_end_)
  FindPairAutoCorrelation(stomp_map, galaxy, random_iterations);
}

void RadialCorrelation::FindCrossCorrelation(Map& stomp_map,
					    CosmoVector& galaxy_z,
					    WAngularVector& galaxy_w,
					    uint8_t random_iterations) {
  //if (!manual_resolution_break_)
  //  AutoMaxResolution(galaxy.size(), stomp_map.Area());

  //if (theta_pixel_begin_ != theta_pixel_end_)
  //  FindPixelAutoCorrelation(stomp_map, galaxy);

  //if (theta_pair_begin_ != theta_pair_end_)
  FindPairCrossCorrelation(stomp_map, galaxy_z, galaxy_w, random_iterations);
}

void RadialCorrelation::FindPairAutoCorrelation(Map& stomp_map,
						CosmoVector& galaxy,
						uint8_t random_iterations) {
  int16_t tree_resolution = min_resolution_;
  if (regionation_resolution_ > min_resolution_)
    tree_resolution = regionation_resolution_;

  TreeMap* galaxy_tree = new TreeMap(tree_resolution, 200);

  uint32_t n_kept = 0;
  uint32_t n_fail = 0;
  for (CosmoIterator iter=galaxy.begin();iter!=galaxy.end();++iter) {
    if (stomp_map.Contains(*iter)) {
      n_kept++;
      if (!galaxy_tree->AddPoint(*iter)) {
	std::cout << "Stomp::RadialCorrelation::FindPairAutoCorrelation - " <<
	  "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
	n_fail++;
      }
    }
  }
  std::cout << "Stomp::RadialCorrelation::FindPairAutoCorrelation - " <<
    n_kept - n_fail << "/" << galaxy.size() << " objects added to tree;" <<
    n_fail << " failed adds...\n";


  if (stomp_map.NRegion() > 0) {
    if (!galaxy_tree->InitializeRegions(stomp_map)) {
      std::cout << "Stomp::RadialCorrelation::FindPairAutoCorrelation - " <<
	"Failed to initialize regions on TreeMap  Exiting.\n";
      exit(2);
    }
  }

  // Galaxy-galaxy
  std::cout << "Stomp::RadialCorrelation::FindPairAutoCorrelation - \n";
  std::cout << "\tGalaxy-galaxy pairs...\n";
  for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
    if (stomp_map.NRegion() > 0) {
      galaxy_tree->FindWeightedPairsWithRegions(galaxy, *iter);
    } else {
      for (CosmoIterator gal_iter=galaxy.begin();
	   gal_iter!=galaxy.end();++gal_iter) {
	iter->SetRedshift(gal_iter->Redshift());
	galaxy_tree->FindWeightedPairs(*gal_iter, *iter);
      }
    }
    iter->MoveWeightToGalGal();
  }

  // Done with the galaxy-based tree, so we can delete that memory.
  delete galaxy_tree;

  // Before we start on the random iterations, we'll zero out the data fields
  // for those counts.
  for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
    iter->ResetGalRand();
    iter->ResetRandGal();
    iter->ResetRandRand();
  }

  std::cout << "Stomp::RadialCorrelation::FindPairAutoCorrelation - \n";
  for (uint8_t rand_iter=0;rand_iter<random_iterations;rand_iter++) {
    std::cout << "\tRandom iteration " <<
      static_cast<int>(rand_iter) << "...\n";

    // Generate set of random points based on the input galaxy file and map.
    CosmoVector random_galaxy;
    stomp_map.GenerateRandomPoints(random_galaxy, galaxy);

    // Create the TreeMap from those random points.
    TreeMap* random_tree = new TreeMap(tree_resolution, 200);

    for (CosmoIterator iter=random_galaxy.begin();
	 iter!=random_galaxy.end();++iter) {
      if (!random_tree->AddPoint(*iter)) {
	std::cout << "Stomp::RadialCorrelation::FindPairAutoCorrelation - " <<
	  "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
      }
    }

    if (stomp_map.NRegion() > 0) {
      if (!random_tree->InitializeRegions(stomp_map)) {
	std::cout << "Stomp::RadialCorrelation::FindPairAutoCorrelation - " <<
	  "Failed to initialize regions on TreeMap  Exiting.\n";
	exit(2);
      }
    }

    // Galaxy-Random -- there's a symmetry here, so the results go in GalRand
    // and RandGal.
    for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
      if (stomp_map.NRegion() > 0) {
	random_tree->FindWeightedPairsWithRegions(galaxy, *iter);
      } else {
	for (CosmoIterator gal_iter=galaxy.begin();
	     gal_iter!=galaxy.end();++gal_iter) {
	  iter->SetRedshift(gal_iter->Redshift());
	  random_tree->FindWeightedPairs(*gal_iter, *iter);
	}
      }
      iter->MoveWeightToGalRand(true);
    }

    // Random-Random
    for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
      if (stomp_map.NRegion() > 0) {
	random_tree->FindWeightedPairsWithRegions(random_galaxy, *iter);
      } else {
	for (CosmoIterator rand_iter=random_galaxy.begin();
	     rand_iter!=random_galaxy.end();++rand_iter) {
	  iter->SetRedshift(rand_iter->Redshift());
	  random_tree->FindWeightedPairs(*rand_iter, *iter);
	}
      }
      iter->MoveWeightToRandRand();
    }

    delete random_tree;
  }

  // Finally, we rescale our random pair counts to normalize them to the
  // number of input objects.
  for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
    iter->RescaleGalRand(1.0*random_iterations);
    iter->RescaleRandGal(1.0*random_iterations);
    iter->RescaleRandRand(1.0*random_iterations);
  }
}

void RadialCorrelation::FindPairCrossCorrelation(Map& stomp_map,
						CosmoVector& galaxy_z,
						WAngularVector& galaxy_w,
						uint8_t random_iterations) {
  int16_t tree_resolution = min_resolution_;
  if (regionation_resolution_ > min_resolution_)
    tree_resolution = regionation_resolution_;

  TreeMap* galaxy_tree = new TreeMap(tree_resolution, 200);

  uint32_t n_kept = 0;
  uint32_t n_fail = 0;
  for (WAngularIterator iter=galaxy_w.begin();iter!=galaxy_w.end();++iter) {
    if (stomp_map.Contains(*iter)) {
      n_kept++;
      if (!galaxy_tree->AddPoint(*iter)) {
	std::cout << "Stomp::RadialCorrelation::FindPairCrossCorrelation - " <<
	  "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
	n_fail++;
      }
    }
  }
  std::cout << "Stomp::RadialCorrelation::FindPairCrossCorrelation - " <<
    n_kept - n_fail << "/" << galaxy_w.size() << " objects added to tree;" <<
    n_fail << " failed adds...\n";


  if (stomp_map.NRegion() > 0) {
    if (!galaxy_tree->InitializeRegions(stomp_map)) {
      std::cout << "Stomp::RadialCorrelation::FindPairCrossCorrelation - " <<
	"Failed to initialize regions on TreeMap  Exiting.\n";
      exit(2);
    }
  }

  // Galaxy-galaxy
  std::cout << "Stomp::RadialCorrelation::FindPairCrossCorrelation - \n";
  std::cout << "\tGalaxy-galaxy pairs...\n";
  for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
    if (stomp_map.NRegion() > 0) {
      galaxy_tree->FindWeightedPairsWithRegions(galaxy_z, *iter);
    } else {
      for (CosmoIterator gal_iter=galaxy_z.begin();
	   gal_iter!=galaxy_z.end();++gal_iter) {
	iter->SetRedshift(gal_iter->Redshift());
	galaxy_tree->FindWeightedPairs(*gal_iter, *iter);
      }
    }
    iter->MoveWeightToGalGal();
  }

  // Before we start on the random iterations, we'll zero out the data fields
  // for those counts.
  for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
    iter->ResetGalRand();
    iter->ResetRandGal();
    iter->ResetRandRand();
  }

  std::cout << "Stomp::RadialCorrelation::FindPairCrossCorrelation - \n";
  for (uint8_t rand_iter=0;rand_iter<random_iterations;rand_iter++) {
    std::cout << "\tRandom iteration " <<
      static_cast<int>(rand_iter) << "...\n";

    // Generate set of random points based on the input galaxy file and map.
    CosmoVector random_galaxy_z;
    WAngularVector random_galaxy_w;
    stomp_map.GenerateRandomPoints(random_galaxy_z, galaxy_z);
    stomp_map.GenerateRandomPoints(random_galaxy_w, galaxy_w);

    // Create the TreeMap from those random points.
    TreeMap* random_tree = new TreeMap(tree_resolution, 200);

    for (WAngularIterator iter=random_galaxy_w.begin();
	 iter!=random_galaxy_w.end();++iter) {
      if (!random_tree->AddPoint(*iter)) {
	std::cout << "Stomp::RadialCorrelation::FindPairCrossCorrelation - " <<
	  "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
      }
    }

    if (stomp_map.NRegion() > 0) {
      if (!random_tree->InitializeRegions(stomp_map)) {
	std::cout << "Stomp::RadialCorrelation::FindPairCrossCorrelation - " <<
	  "Failed to initialize regions on TreeMap  Exiting.\n";
	exit(2);
      }
    }

    // Galaxy-Random
    for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
      if (stomp_map.NRegion() > 0) {
	random_tree->FindWeightedPairsWithRegions(galaxy_z, *iter);
      } else {
	for (CosmoIterator gal_iter=galaxy_z.begin();
	     gal_iter!=galaxy_z.end();++gal_iter) {
	  iter->SetRedshift(gal_iter->Redshift());
	  random_tree->FindWeightedPairs(*gal_iter, *iter);
	}
      }
      iter->MoveWeightToGalRand();
    }

    // Random-Galaxy
    for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
      if (stomp_map.NRegion() > 0) {
	galaxy_tree->FindWeightedPairsWithRegions(random_galaxy_z, *iter);
      } else {
	for (CosmoIterator gal_iter=random_galaxy_z.begin();
	     gal_iter!=random_galaxy_z.end();++gal_iter) {
	  iter->SetRedshift(gal_iter->Redshift());
	  random_tree->FindWeightedPairs(*gal_iter, *iter);
	}
      }
      iter->MoveWeightToRandGal();
    }

    // Random-Random
    for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
      if (stomp_map.NRegion() > 0) {
	random_tree->FindWeightedPairsWithRegions(random_galaxy_z, *iter);
      } else {
	for (CosmoIterator rand_iter=random_galaxy_z.begin();
	     rand_iter!=random_galaxy_z.end();++rand_iter) {
	  iter->SetRedshift(rand_iter->Redshift());
	  random_tree->FindWeightedPairs(*rand_iter, *iter);
	}
      }
      iter->MoveWeightToRandRand();
    }

    delete random_tree;
  }

  // Done with the galaxy-based tree, so we can delete that memory.
  delete galaxy_tree;

  // Finally, we rescale our random pair counts to normalize them to the
  // number of input objects.
  for (RadialIterator iter=radial_pair_begin_;iter!=radial_pair_end_;++iter) {
    iter->RescaleGalRand(1.0*random_iterations);
    iter->RescaleRandGal(1.0*random_iterations);
    iter->RescaleRandRand(1.0*random_iterations);
  }
}

void RadialCorrelation::FindAutoCorrelationWithRegions(Map& stomp_map,
						       CosmoVector& gal,
						       uint8_t random_iter,
						       uint16_t n_regions) {

  if (n_regions == 0) n_regions = static_cast<uint16_t>(2*radialbin_.size());
  std::cout << "Stomp::RadialCorrelation::FindAutoCorrelationWithRegions - " <<
    "Regionating with " << n_regions << " regions...\n";
  uint16_t n_true_regions = stomp_map.NRegion();
  if (n_true_regions == 0)
    n_true_regions = stomp_map.InitializeRegions(n_regions);

  if (n_true_regions != n_regions) {
    std::cout << "Stomp::RadialCorrelation::" <<
      "FindAutoCorrelationWithRegions - Splitting into " << n_true_regions <<
      " rather than " << n_regions << "...\n";
    n_regions = n_true_regions;
  }

  regionation_resolution_ = stomp_map.RegionResolution();

  std::cout << "Stomp::RadialCorrelation::FindAutoCorrelationWithRegions - " <<
    "Regionated at " << regionation_resolution_ << "...\n";
  InitializeRegions(n_regions);
  if (regionation_resolution_ > min_resolution_)
    SetMinResolution(regionation_resolution_);

  if (regionation_resolution_ > max_resolution_) {
    std::cout << "Stomp::RadialCorrelation::FindAutoCorrelationWithRegions " <<
      " - regionation resolution (" << regionation_resolution_ <<
      ") exceeds maximum resolution (" << max_resolution_ << ")\n" <<
      "\tReseting to use pair-based estimator only\n";
    UseOnlyPairs();
  }

  FindPairAutoCorrelation(stomp_map, gal, random_iter);
}

void RadialCorrelation::FindCrossCorrelationWithRegions(Map& stomp_map,
						       CosmoVector& galaxy_z,
						       WAngularVector& galaxy_w,
						       uint8_t random_iter,
						       uint16_t n_regions) {

  if (n_regions == 0) n_regions = static_cast<uint16_t>(2*radialbin_.size());
  std::cout << "Stomp::RadialCorrelation::FindCrossCorrelationWithRegions - " <<
    "Regionating with " << n_regions << " regions...\n";
  uint16_t n_true_regions = stomp_map.NRegion();
  if (n_true_regions == 0)
    n_true_regions = stomp_map.InitializeRegions(n_regions);

  if (n_true_regions != n_regions) {
    std::cout << "Stomp::RadialCorrelation::" <<
      "FindAutoCorrelationWithRegions - Splitting into " << n_true_regions <<
      " rather than " << n_regions << "...\n";
    n_regions = n_true_regions;
  }

  regionation_resolution_ = stomp_map.RegionResolution();

  std::cout << "Stomp::RadialCorrelation::FindCrossCorrelationWithRegions - " <<
    "Regionated at " << regionation_resolution_ << "...\n";
  InitializeRegions(n_regions);
  if (regionation_resolution_ > min_resolution_)
    SetMinResolution(regionation_resolution_);

  if (regionation_resolution_ > max_resolution_) {
    std::cout << "Stomp::RadialCorrelation::FindCrossCorrelationWithRegions " <<
      " - regionation resolution (" << regionation_resolution_ <<
      ") exceeds maximum resolution (" << max_resolution_ << ")\n" <<
      "\tReseting to use pair-based estimator only\n";
    UseOnlyPairs();
  }

  FindPairCrossCorrelation(stomp_map, galaxy_z, galaxy_w, random_iter);
}

bool RadialCorrelation::Write(const std::string& output_file_name) {
  bool wrote_file = false;

  std::ofstream output_file(output_file_name.c_str());

  if (output_file.is_open()) {
    wrote_file = true;

    for (RadialIterator iter=Begin();iter!=End();++iter) {
      if (iter->NRegion()) {
	output_file << std::setprecision(6) << iter->Radius() << " " <<
	  iter->MeanWtheta()  << " " << iter->MeanWthetaError() << "\n";
      } else {
	if (iter->Resolution() == 0) {
	  output_file << std::setprecision(6) << iter->Radius() << " " <<
	    iter->Wtheta()  << " " << iter->GalGal() << " " <<
	    iter->GalRand() << " " << iter->RandGal() << " " <<
	    iter->RandRand() << "\n";
	} else {
	  output_file << std::setprecision(6) << iter->Radius() << " " <<
	    iter->Wtheta()  << " " << iter->PixelWtheta() << " " <<
	    iter->PixelWeight() << "\n";
	}
      }
    }

    output_file.close();
  }

  return wrote_file;
}

double RadialCorrelation::Covariance(uint8_t bin_idx_a, uint8_t bin_idx_b) {
  double covariance = 0.0;

  RadialIterator r_a = BinIterator(bin_idx_a);
  RadialIterator r_b = BinIterator(bin_idx_b);

  if ((r_a->NRegion() == r_b->NRegion()) &&
      (r_a->NRegion() > 0)) {
    // We have a valid number of regions and both bins were calculated with
    // the same number of regions, so we calculate the jack-knife covariance.
    uint16_t n_region = r_a->NRegion();
    double mean_wtheta_a = r_a->MeanWtheta();
    double mean_wtheta_b = r_b->MeanWtheta();

    for (uint16_t region_iter=0;region_iter<n_region;region_iter++) {
      covariance +=
	(r_a->Wtheta(region_iter) - mean_wtheta_a)*
	(r_b->Wtheta(region_iter) - mean_wtheta_b);
    }

    covariance *= (n_region - 1.0)*(n_region - 1.0)/(1.0*n_region*n_region);
  } else {
    // If the above doesn't hold, then we're reduced to using Poisson errors.
    // In this case, we only have non-zero elements on the diagonal of the
    // covariance matrix; all others are zero by definition.
    if (bin_idx_a == bin_idx_b) {
      covariance = r_a->WthetaError()*r_a->WthetaError();
    }
  }

  return covariance;
}

bool RadialCorrelation::WriteCovariance(const std::string& output_file_name) {
  bool wrote_file = false;

  std::ofstream output_file(output_file_name.c_str());

  if (output_file.is_open()) {
    wrote_file = true;

    for (uint8_t r_idx_a=0;r_idx_a<radialbin_.size();r_idx_a++) {
      for (uint8_t r_idx_b=0;r_idx_b<radialbin_.size();r_idx_b++) {
	output_file << std::setprecision(6) <<
	  radialbin_[r_idx_a].Radius() << " " <<
	  radialbin_[r_idx_b].Radius() << " " <<
	  Covariance(r_idx_a, r_idx_b) << "\n";
      }
    }

    output_file.close();
  }

  return wrote_file;
}

void RadialCorrelation::UseOnlyPairs() {
  radial_pixel_begin_ = radialbin_.end();
  radial_pixel_end_ = radialbin_.end();

  radial_pair_begin_ = radialbin_.begin();
  radial_pair_end_ = radialbin_.end();
  for (RadialIterator iter=radialbin_.begin();iter!=radialbin_.end();++iter) {
    iter->SetResolution(0);
  }
  manual_resolution_break_ = true;
}

RadialIterator RadialCorrelation::Begin(uint32_t resolution) {
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      return radial_pair_begin_;
    } else {
      return radialbin_.begin();
    }
  } else {
    RadialBin radius;
    radius.SetResolution(resolution);
    RadialPair iter = equal_range(radial_pixel_begin_, radialbin_.end(),
				  radius, RadialBin::ReverseResolutionOrder);
    return iter.first;
  }
}

RadialIterator RadialCorrelation::End(uint32_t resolution) {
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      return radial_pair_end_;
    } else {
      return radialbin_.end();
    }
  } else {
    RadialBin radius;
    radius.SetResolution(resolution);
    RadialPair iter = equal_range(radial_pixel_begin_,radialbin_.end(),
				  radius,RadialBin::ReverseResolutionOrder);
    return iter.second;
  }
}

RadialIterator RadialCorrelation::BinIterator(uint8_t bin_idx) {
  RadialIterator iter = radialbin_.begin();

  for (uint8_t i=0;i<bin_idx;++i) {
    if (iter != radialbin_.end()) ++iter;
  }

  return iter;
}

uint32_t RadialCorrelation::NBins() {
  return radialbin_.size();
}

int16_t RadialCorrelation::NRegion() {
	return n_region_;
}

double RadialCorrelation::RadiusMin(uint32_t resolution) {
  double radius_min = -1.0;
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
    	radius_min = radial_pair_begin_->ThetaMin();
    } else {
      radius_min = r_min_;
    }
  } else {
    RadialBin radius;
    radius.SetResolution(resolution);
    RadialPair iter = equal_range(radial_pixel_begin_,radialbin_.end(),
				 radius,AngularBin::ReverseResolutionOrder);
    if (iter.first != iter.second) {
      radius_min = iter.first->RadiusMin();
    }
  }

  return radius_min;
}

double RadialCorrelation::RadiusMax(uint32_t resolution) {
  double radius_max = -1.0;
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      radius_max = radial_pair_end_->ThetaMin();
    } else {
      radius_max = r_max_;
    }
  } else {
    RadialBin radius;
    radius.SetResolution(resolution);
    RadialPair iter = equal_range(radial_pixel_begin_,radialbin_.end(),
				 radius,AngularBin::ReverseResolutionOrder);
    if (iter.first != iter.second) {
      --iter.second;
      radius_max = iter.second->RadiusMax();
    }
  }

  return radius_max;
}

} // end namespace Stomp
