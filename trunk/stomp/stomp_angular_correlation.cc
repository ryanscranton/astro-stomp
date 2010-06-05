// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the class for calculating angular correlations on
// the sphere.  In general, different methods are more efficient on small vs.
// large angular scales, so this class draws on nearly the entire breadth of
// the STOMP library.

#include "stomp_core.h"
#include "stomp_angular_correlation.h"
#include "stomp_map.h"
#include "stomp_scalar_map.h"
#include "stomp_tree_map.h"

namespace Stomp {

AngularCorrelation::AngularCorrelation(double theta_min, double theta_max,
				       double bins_per_decade,
				       bool assign_resolutions) {
  double unit_double = floor(log10(theta_min))*bins_per_decade;
  double theta = pow(10.0, unit_double/bins_per_decade);

  while (theta < theta_max) {
    if (DoubleGE(theta, theta_min) && (theta < theta_max)) {
      AngularBin thetabin;
      thetabin.SetThetaMin(theta);
      thetabin.SetThetaMax(pow(10.0,(unit_double+1.0)/bins_per_decade));
      thetabin.SetTheta(pow(10.0,0.5*(log10(thetabin.ThetaMin())+
				      log10(thetabin.ThetaMax()))));
      thetabin_.push_back(thetabin);
    }
    unit_double += 1.0;
    theta = pow(10.0,unit_double/bins_per_decade);
  }

  theta_min_ = thetabin_[0].ThetaMin();
  sin2theta_min_ = thetabin_[0].Sin2ThetaMin();
  theta_max_ = thetabin_[thetabin_.size()-1].ThetaMax();
  sin2theta_max_ = thetabin_[thetabin_.size()-1].Sin2ThetaMax();

  if (assign_resolutions) {
    AssignBinResolutions();
    theta_pixel_begin_ = thetabin_.begin();
    theta_pair_begin_ = thetabin_.begin();
    theta_pair_end_ = thetabin_.begin();
  } else {
    min_resolution_ = HPixResolution;
    max_resolution_ = HPixResolution;
    theta_pixel_begin_ = thetabin_.end();
    theta_pair_begin_ = thetabin_.begin();
    theta_pair_end_ = thetabin_.end();
  }

  manual_resolution_break_ = false;
}

AngularCorrelation::AngularCorrelation(uint32_t n_bins,
				       double theta_min, double theta_max,
				       bool assign_resolutions) {
  double dtheta = (theta_max - theta_min)/n_bins;

  for (uint32_t i=0;i<n_bins;i++) {
    AngularBin thetabin;
    thetabin.SetThetaMin(theta_min + i*dtheta);
    thetabin.SetThetaMin(theta_min + (i+1)*dtheta);
    thetabin.SetTheta(0.5*(thetabin.ThetaMin()+thetabin.ThetaMax()));
    thetabin_.push_back(thetabin);
  }

  theta_min_ = thetabin_[0].ThetaMin();
  sin2theta_min_ = thetabin_[0].Sin2ThetaMin();
  theta_max_ = thetabin_[n_bins-1].ThetaMax();
  sin2theta_max_ = thetabin_[n_bins-1].Sin2ThetaMax();

  if (assign_resolutions) {
    AssignBinResolutions();
    theta_pixel_begin_ = thetabin_.begin();
    theta_pair_begin_ = thetabin_.begin();
    theta_pair_end_ = thetabin_.begin();
  } else {
    min_resolution_ = HPixResolution;
    max_resolution_ = HPixResolution;
    theta_pixel_begin_ = thetabin_.end();
    theta_pair_begin_ = thetabin_.begin();
    theta_pair_end_ = thetabin_.end();
  }

  manual_resolution_break_ = false;
}

void AngularCorrelation::AssignBinResolutions(double lammin, double lammax,
					      uint32_t max_resolution) {
  min_resolution_ = MaxPixelResolution;
  max_resolution_ = HPixResolution;

  for (ThetaIterator iter=thetabin_.begin();iter!=thetabin_.end();++iter) {
    iter->CalculateResolution(lammin, lammax, max_resolution);

    if (iter->Resolution() < min_resolution_)
      min_resolution_ = iter->Resolution();
    if (iter->Resolution() > max_resolution_)
      max_resolution_ = iter->Resolution();
  }
}

void AngularCorrelation::SetMaxResolution(uint32_t resolution,
					  bool manual_break) {
  max_resolution_ = resolution;

  // By default every bin is calculated with the pixel-based estimator
  theta_pair_begin_ = thetabin_.begin();
  theta_pair_end_ = thetabin_.begin();
  theta_pixel_begin_ = thetabin_.begin();

  for (ThetaIterator iter=thetabin_.begin();iter!=thetabin_.end();++iter) {
    iter->CalculateResolution();
    if (iter->Resolution() > max_resolution_) {
      iter->SetResolution(0);
      ++theta_pixel_begin_;
      ++theta_pair_end_;
    }
  }

  if (manual_break) manual_resolution_break_ = true;
}

void AngularCorrelation::SetMinResolution(uint32_t resolution) {
  min_resolution_ = resolution;
  for (ThetaIterator iter=theta_pixel_begin_;iter!=thetabin_.end();++iter) {
    if (iter->Resolution() < min_resolution_) {
      iter->SetResolution(min_resolution_);
    }
  }
}

void AngularCorrelation::AutoMaxResolution(uint32_t n_obj, double area) {
  uint32_t max_resolution = 2048;

  if (area > 500.0) {
    // large survey limit
    max_resolution = 512;
    if (n_obj < 500000) max_resolution = 64;
    if ((n_obj > 500000) && (n_obj < 2000000)) max_resolution = 128;
    if ((n_obj > 2000000) && (n_obj < 10000000)) max_resolution = 256;
  } else {
    // small survey limit
    if (n_obj < 500000) max_resolution = 256;
    if ((n_obj > 500000) && (n_obj < 2000000)) max_resolution = 512;
    if ((n_obj > 2000000) && (n_obj < 10000000)) max_resolution = 1024;
  }

  std::cout << "Setting maximum resolution to " <<
    static_cast<int>(max_resolution) << "...\n";

  SetMaxResolution(max_resolution, false);
}

void AngularCorrelation::FindAutoCorrelation(Map& stomp_map,
					     WAngularVector& galaxy,
					     uint8_t random_iterations) {
  if (!manual_resolution_break_)
    AutoMaxResolution(galaxy.size(), stomp_map.Area());

  FindPixelAutoCorrelation(stomp_map, galaxy);
  FindPairAutoCorrelation(stomp_map, galaxy, random_iterations);
}

void AngularCorrelation::FindCrossCorrelation(Map& stomp_map,
					      WAngularVector& galaxy_a,
					      WAngularVector& galaxy_b,
					      uint8_t random_iterations) {
  if (!manual_resolution_break_) {
    uint32_t n_obj =
      static_cast<uint32_t>(sqrt(1.0*galaxy_a.size()*galaxy_b.size()));
    AutoMaxResolution(n_obj, stomp_map.Area());
  }

  FindPixelCrossCorrelation(stomp_map, galaxy_a, galaxy_b);
  FindPairCrossCorrelation(stomp_map, galaxy_a, galaxy_b, random_iterations);
}

void AngularCorrelation::FindAutoCorrelationWithRegions(Map& stomp_map,
							WAngularVector& gal,
							uint8_t random_iter,
							uint16_t n_regions) {
  if (!manual_resolution_break_)
    AutoMaxResolution(gal.size(), stomp_map.Area());

  if (n_regions == 0) n_regions = static_cast<uint16_t>(2*thetabin_.size());
  std::cout << "Regionating with " << n_regions << " regions...\n";
  uint16_t n_true_regions = stomp_map.InitializeRegions(n_regions);
  if (n_true_regions != n_regions) {
    std::cout << "Splitting into " << n_true_regions << " rather than " <<
      n_regions << "...\n";
    n_regions = n_true_regions;
  }

  std::cout << "Regionated at " << stomp_map.RegionResolution() << "...\n";
  for (ThetaIterator iter=Begin();iter!=End();++iter)
    iter->InitializeRegions(n_regions);
  SetMinResolution(stomp_map.RegionResolution());

  std::cout << "Auto-correlating with pixels...\n";
  FindPixelAutoCorrelation(stomp_map, gal);
  std::cout << "Auto-correlating with pairs...\n";
  FindPairAutoCorrelation(stomp_map, gal, random_iter);
}

void AngularCorrelation::FindCrossCorrelationWithRegions(Map& stomp_map,
							 WAngularVector& gal_a,
							 WAngularVector& gal_b,
							 uint8_t random_iter,
							 uint16_t n_regions) {
  if (!manual_resolution_break_) {
    uint32_t n_obj =
      static_cast<uint32_t>(sqrt(1.0*gal_a.size()*gal_b.size()));
    AutoMaxResolution(n_obj, stomp_map.Area());
  }

  if (n_regions == 0) n_regions = static_cast<uint16_t>(2*thetabin_.size());
  uint16_t n_true_regions = stomp_map.InitializeRegions(n_regions);
  if (n_true_regions != n_regions) {
    std::cout << "Splitting into " << n_true_regions << " rather than " <<
      n_regions << "...\n";
    n_regions = n_true_regions;
  }

  std::cout << "Regionated at " << stomp_map.RegionResolution() << "...\n";
  for (ThetaIterator iter=Begin();iter!=End();++iter)
    iter->InitializeRegions(n_regions);
  SetMinResolution(stomp_map.RegionResolution());

  FindPixelCrossCorrelation(stomp_map, gal_a, gal_b);
  FindPairCrossCorrelation(stomp_map, gal_a, gal_b, random_iter);
}

void AngularCorrelation::FindPixelAutoCorrelation(Map& stomp_map,
						  WAngularVector& galaxy) {

  std::cout << "Initializing ScalarMap at " << max_resolution_ << "...\n";
  ScalarMap* scalar_map = new ScalarMap(stomp_map, max_resolution_,
					ScalarMap::DensityField);
  if (stomp_map.NRegion() > 0) {
    std::cout << "Intializing regions...\n";
    scalar_map->InitializeRegions(stomp_map);
  }

  std::cout << "Adding points to ScalarMap...\n";
  uint32_t n_filtered = 0;
  uint32_t n_kept = 0;
  for (WAngularIterator iter=galaxy.begin();iter!=galaxy.end();++iter) {
    if (stomp_map.Contains(*iter)) {
      n_filtered++;
      if (scalar_map->AddToMap(*iter)) n_kept++;
    }
  }

  if (n_filtered != galaxy.size())
    std::cout << "WARNING: " << galaxy.size() - n_filtered <<
      "/" << galaxy.size() << " objects not within input Map.\n";

  if (n_filtered != n_kept)
    std::cout << "WARNING: Failed to place " << n_filtered - n_kept <<
      "/" << n_filtered << " filtered objects into ScalarMap.\n";

  FindPixelAutoCorrelation(*scalar_map);

  delete scalar_map;
}

void AngularCorrelation::FindPixelAutoCorrelation(ScalarMap& stomp_map) {
  for (ThetaIterator iter=Begin(stomp_map.Resolution());
       iter!=End(stomp_map.Resolution());++iter) {
    if (stomp_map.NRegion() > 0) {
      std::cout << "\tAuto-correlating with regions at " <<
	stomp_map.Resolution() << "...\n";
      stomp_map.AutoCorrelateWithRegions(iter);
    } else {
      stomp_map.AutoCorrelate(iter);
    }
  }

  for (uint32_t resolution=stomp_map.Resolution()/2;
       resolution>=min_resolution_;resolution/=2) {
    ScalarMap* sub_scalar_map =
      new ScalarMap(stomp_map,resolution);
    if (stomp_map.NRegion() > 0) sub_scalar_map->InitializeRegions(stomp_map);
    for (ThetaIterator iter=Begin(resolution);iter!=End(resolution);++iter) {
      if (stomp_map.NRegion() > 0) {
	std::cout << "\tAuto-correlating with regions at " <<
	  sub_scalar_map->Resolution() << "...\n";
	sub_scalar_map->AutoCorrelateWithRegions(iter);
      } else {
	std::cout << "\tAuto-correlating at " <<
	  sub_scalar_map->Resolution() << "...\n";
	sub_scalar_map->AutoCorrelate(iter);
      }
    }
    delete sub_scalar_map;
  }
}

void AngularCorrelation::FindPixelCrossCorrelation(Map& stomp_map,
						   WAngularVector& galaxy_a,
						   WAngularVector& galaxy_b) {
  std::cout << "Initialing ScalarMaps at " << max_resolution_ << "...\n";
  ScalarMap* scalar_map_a = new ScalarMap(stomp_map, max_resolution_,
					  ScalarMap::DensityField);
  ScalarMap* scalar_map_b = new ScalarMap(stomp_map, max_resolution_,
					  ScalarMap::DensityField);
  if (stomp_map.NRegion() > 0) {
    std::cout << "Intializing regions...\n";
    scalar_map_a->InitializeRegions(stomp_map);
    scalar_map_b->InitializeRegions(stomp_map);
  }

  uint32_t n_filtered = 0;
  uint32_t n_kept = 0;
  for (WAngularIterator iter=galaxy_a.begin();iter!=galaxy_a.end();++iter) {
    if (stomp_map.Contains(*iter)) {
      n_filtered++;
      if (scalar_map_a->AddToMap(*iter)) n_kept++;
    }
  }

  if (n_filtered != galaxy_a.size())
    std::cout << "WARNING: " << galaxy_a.size() - n_filtered <<
      "/" << galaxy_a.size() << " objects not within input Map.\n";
  if (n_filtered != n_kept)
    std::cout << "WARNING: Failed to place " << n_filtered - n_kept <<
      "/" << n_filtered << " filtered objects into ScalarMap.\n";

  n_filtered = 0;
  n_kept = 0;
  for (WAngularIterator iter=galaxy_b.begin();iter!=galaxy_b.end();++iter) {
    if (stomp_map.Contains(*iter)) {
      n_filtered++;
      if (scalar_map_b->AddToMap(*iter)) n_kept++;
    }
  }

  if (n_filtered != galaxy_b.size())
    std::cout << "WARNING: " << galaxy_b.size() - n_filtered <<
      "/" << galaxy_b.size() << " objects not within input Map.\n";
  if (n_filtered != n_kept)
    std::cout << "WARNING: Failed to place " << n_filtered - n_kept <<
      "/" << n_filtered << " filtered objects into ScalarMap.\n";

  FindPixelCrossCorrelation(*scalar_map_a, *scalar_map_b);

  delete scalar_map_a;
  delete scalar_map_b;
}

void AngularCorrelation::FindPixelCrossCorrelation(ScalarMap& map_a,
						   ScalarMap& map_b) {
  if (map_a.Resolution() != map_b.Resolution()) {
    std::cout << "Incompatible density map resolutions.  Exiting!\n";
    exit(1);
  }

  for (ThetaIterator iter=Begin(map_a.Resolution());
       iter!=End(map_a.Resolution());++iter) {
    if (map_a.NRegion() > 0) {
      map_a.CrossCorrelateWithRegions(map_b, iter);
    } else {
      map_a.CrossCorrelate(map_b, iter);
    }
  }

  for (uint32_t resolution=map_a.Resolution()/2;
       resolution>=min_resolution_;resolution/=2) {
    ScalarMap* sub_map_a = new ScalarMap(map_a, resolution);
    ScalarMap* sub_map_b = new ScalarMap(map_b, resolution);

    if (map_a.NRegion() > 0) {
      sub_map_a->InitializeRegions(map_a);
      sub_map_b->InitializeRegions(map_a);
    }

    for (ThetaIterator iter=Begin(resolution);iter!=End(resolution);++iter) {
      if (map_a.NRegion() > 0) {
	std::cout << "\tCross-correlating with regions at " <<
	  sub_map_a->Resolution() << "...\n";
	sub_map_a->CrossCorrelateWithRegions(*sub_map_b, iter);
      } else {
	std::cout << "\tCross-correlating at " <<
	  sub_map_a->Resolution() << "...\n";
	sub_map_a->CrossCorrelate(*sub_map_b, iter);
      }
    }
    delete sub_map_a;
    delete sub_map_b;
  }
}

void AngularCorrelation::FindPairAutoCorrelation(Map& stomp_map,
						 WAngularVector& galaxy,
						 uint8_t random_iterations) {
  TreeMap* galaxy_tree = new TreeMap(min_resolution_, 200);

  uint32_t n_kept = 0;
  uint32_t n_fail = 0;
  for (WAngularIterator iter=galaxy.begin();iter!=galaxy.end();++iter) {
    if (stomp_map.Contains(*iter)) {
      n_kept++;
      if (!galaxy_tree->AddPoint(*iter)) {
	std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
	n_fail++;
      }
    }
  }
  std::cout << n_kept - n_fail << "/" << galaxy.size() <<
    " objects added to tree;" << n_fail << " failed adds...\n";


  if (stomp_map.NRegion() > 0) {
    if (!galaxy_tree->InitializeRegions(stomp_map)) {
      std::cout << "Failed to initialize regions on TreeMap  Exiting.\n";
      exit(2);
    }
  }

  // Galaxy-galaxy
  std::cout << "\tGalaxy-galaxy pairs...\n";
  for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
    if (stomp_map.NRegion() > 0) {
      galaxy_tree->FindWeightedPairsWithRegions(galaxy, *iter);
    } else {
      galaxy_tree->FindWeightedPairs(galaxy, *iter);
    }
    iter->MoveWeightToGalGal();
  }

  // Done with the galaxy-based tree, so we can delete that memory.
  delete galaxy_tree;

  // Before we start on the random iterations, we'll zero out the data fields
  // for those counts.
  for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
    iter->ResetGalRand();
    iter->ResetRandGal();
    iter->ResetRandRand();
  }

  for (uint8_t rand_iter=0;rand_iter<random_iterations;rand_iter++) {
    std::cout << "\tRandom iteration " <<
      static_cast<int>(rand_iter) << "...\n";

    // Generate set of random points based on the input galaxy file and map.
    WAngularVector random_galaxy;
    stomp_map.GenerateRandomPoints(random_galaxy, galaxy, true);

    // Create the TreeMap from those random points.
    TreeMap* random_tree = new TreeMap(min_resolution_, 200);

    for (WAngularIterator iter=random_galaxy.begin();
	 iter!=random_galaxy.end();++iter) {
      if (!random_tree->AddPoint(*iter)) {
	std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
      }
    }

    if (stomp_map.NRegion() > 0) {
      if (!random_tree->InitializeRegions(stomp_map)) {
	std::cout << "Failed to initialize regions on TreeMap  Exiting.\n";
	exit(2);
      }
    }

    // Galaxy-Random -- there's a symmetry here, so the results go in GalRand
    // and RandGal.
    for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
      if (stomp_map.NRegion() > 0) {
	random_tree->FindWeightedPairsWithRegions(galaxy, *iter);
      } else {
	random_tree->FindWeightedPairs(galaxy, *iter);
      }
      iter->MoveWeightToGalRand(true);
    }

    // Random-Random
    for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
      if (stomp_map.NRegion() > 0) {
	random_tree->FindWeightedPairsWithRegions(random_galaxy, *iter);
      } else {
	random_tree->FindWeightedPairs(random_galaxy, *iter);
      }
      iter->MoveWeightToRandRand();
    }

    delete random_tree;
  }

  // Finally, we rescale our random pair counts to normalize them to the
  // number of input objects.
  for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
    iter->RescaleGalRand(1.0*random_iterations);
    iter->RescaleRandGal(1.0*random_iterations);
    iter->RescaleRandRand(1.0*random_iterations);
  }
}

void AngularCorrelation::FindPairCrossCorrelation(Map& stomp_map,
						  WAngularVector& galaxy_a,
						  WAngularVector& galaxy_b,
						  uint8_t random_iterations) {
  TreeMap* galaxy_tree_a = new TreeMap(min_resolution_, 200);

  uint32_t n_kept = 0;
  uint32_t n_fail = 0;
  for (WAngularIterator iter=galaxy_a.begin();iter!=galaxy_a.end();++iter) {
    if (stomp_map.Contains(*iter)) {
      n_kept++;
      if (!galaxy_tree_a->AddPoint(*iter)) {
	std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
	n_fail++;
      }
    }
  }
  std::cout << n_kept - n_fail << "/" << galaxy_a.size() <<
    " objects added to tree;" << n_fail << " failed adds...\n";

  if (stomp_map.NRegion() > 0) {
    if (!galaxy_tree_a->InitializeRegions(stomp_map)) {
      std::cout << "Failed to initialize regions on TreeMap  Exiting.\n";
      exit(2);
    }
  }

  // Galaxy-galaxy
  std::cout << "\tGalaxy-galaxy pairs...\n";
  for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
    if (stomp_map.NRegion() > 0) {
      galaxy_tree_a->FindWeightedPairsWithRegions(galaxy_b, *iter);
    } else {
      galaxy_tree_a->FindWeightedPairs(galaxy_b, *iter);
    }
    // If the number of random iterations is 0, then we're doing a
    // WeightedCrossCorrelation instead of a cross-correlation between 2
    // population densities.  In that case, we want the ratio between the
    // WeightedPairs and Pairs, so we keep the values in the Weight and
    // Counter fields.
    if (random_iterations > 0) iter->MoveWeightToGalGal();
  }

  // Before we start on the random iterations, we'll zero out the data fields
  // for those counts.
  for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
    iter->ResetGalRand();
    iter->ResetRandGal();
    iter->ResetRandRand();
  }

  for (uint8_t rand_iter=0;rand_iter<random_iterations;rand_iter++) {
    std::cout << "\tRandom iteration " <<
      static_cast<int>(rand_iter) << "...\n";
    WAngularVector random_galaxy_a;
    stomp_map.GenerateRandomPoints(random_galaxy_a, galaxy_a, true);

    WAngularVector random_galaxy_b;
    stomp_map.GenerateRandomPoints(random_galaxy_b, galaxy_b, true);

    // Galaxy-Random
    for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
      if (stomp_map.NRegion() > 0) {
	galaxy_tree_a->FindWeightedPairsWithRegions(random_galaxy_b, *iter);
      } else {
	galaxy_tree_a->FindWeightedPairs(random_galaxy_b, *iter);
      }
      iter->MoveWeightToGalRand();
    }

    TreeMap* random_tree_a = new TreeMap(min_resolution_, 200);

    for (WAngularIterator iter=random_galaxy_a.begin();
	 iter!=random_galaxy_a.end();++iter) {
      if (!random_tree_a->AddPoint(*iter)) {
	std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
      }
    }

    if (stomp_map.NRegion() > 0) {
      if (!random_tree_a->InitializeRegions(stomp_map)) {
	std::cout << "Failed to initialize regions on TreeMap  Exiting.\n";
	exit(2);
      }
    }

    // Random-Galaxy
    for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
      if (stomp_map.NRegion() > 0) {
	random_tree_a->FindWeightedPairsWithRegions(galaxy_b, *iter);
      } else {
	random_tree_a->FindWeightedPairs(galaxy_b, *iter);
      }
      iter->MoveWeightToRandGal();
    }

    // Random-Random
    for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
      if (stomp_map.NRegion() > 0) {
	random_tree_a->FindWeightedPairsWithRegions(random_galaxy_b, *iter);
      } else {
	random_tree_a->FindWeightedPairs(random_galaxy_b, *iter);
      }
      iter->MoveWeightToRandRand();
    }

    delete random_tree_a;
  }

  delete galaxy_tree_a;

  // Finally, we rescale our random pair counts to normalize them to the
  // number of input objects.
  for (ThetaIterator iter=Begin(0);iter!=End(0);++iter) {
    iter->RescaleGalRand(1.0*random_iterations);
    iter->RescaleRandGal(1.0*random_iterations);
    iter->RescaleRandRand(1.0*random_iterations);
  }
}

bool AngularCorrelation::Write(const std::string& output_file_name) {
  bool wrote_file = false;

  std::ofstream output_file(output_file_name.c_str());

  if (output_file.is_open()) {
    wrote_file = true;

    for (ThetaIterator iter=Begin();iter!=End();++iter) {
      if (iter->NRegion()) {
	output_file << std::setprecision(6) << iter->Theta() << " " <<
	  iter->MeanWtheta()  << " " << iter->MeanWthetaError() << "\n";
      } else {
	if (iter->Resolution() == 0) {
	  output_file << std::setprecision(6) << iter->Theta() << " " <<
	    iter->Wtheta()  << " " << iter->GalGal() << " " <<
	    iter->GalRand() << " " << iter->RandGal() << " " <<
	    iter->RandRand() << "\n";
	} else {
	  output_file << std::setprecision(6) << iter->Theta() << " " <<
	    iter->Wtheta()  << " " << iter->PixelWtheta() << " " <<
	    iter->PixelWeight() << "\n";
	}
      }
    }

    output_file.close();
  }

  return wrote_file;
}

double AngularCorrelation::ThetaMin(uint32_t resolution) {
  double theta_min = -1.0;
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      theta_min = theta_pair_begin_->ThetaMin();
    } else {
      theta_min = theta_min_;
    }
  } else {
    AngularBin theta;
    theta.SetResolution(resolution);
    ThetaPair iter = equal_range(theta_pixel_begin_,thetabin_.end(),
				 theta,AngularBin::ReverseResolutionOrder);
    if (iter.first != iter.second) {
      theta_min = iter.first->ThetaMin();
    }
  }

  return theta_min;
}

double AngularCorrelation::ThetaMax(uint32_t resolution) {
  double theta_max = -1.0;
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      theta_max = theta_pair_end_->ThetaMin();
    } else {
      theta_max = theta_max_;
    }
  } else {
    AngularBin theta;
    theta.SetResolution(resolution);
    ThetaPair iter = equal_range(theta_pixel_begin_,thetabin_.end(),
				 theta,AngularBin::ReverseResolutionOrder);
    if (iter.first != iter.second) {
      --iter.second;
      theta_max = iter.second->ThetaMax();
    }
  }

  return theta_max;
}

double AngularCorrelation::Sin2ThetaMin(uint32_t resolution) {
  double sin2theta_min = -1.0;
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      sin2theta_min = theta_pair_begin_->Sin2ThetaMin();
    } else {
      sin2theta_min = sin2theta_min_;
    }
  } else {
    AngularBin theta;
    theta.SetResolution(resolution);
    ThetaPair iter = equal_range(theta_pixel_begin_,thetabin_.end(),
				 theta,AngularBin::ReverseResolutionOrder);
    if (iter.first != iter.second) {
      sin2theta_min = iter.first->Sin2ThetaMin();
    }
  }

  return sin2theta_min;
}

double AngularCorrelation::Sin2ThetaMax(uint32_t resolution) {
  double sin2theta_max = -1.0;
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      sin2theta_max = theta_pair_end_->Sin2ThetaMin();
    } else {
      sin2theta_max = sin2theta_max_;
    }
  } else {
    AngularBin theta;
    theta.SetResolution(resolution);
    ThetaPair iter = equal_range(theta_pixel_begin_,thetabin_.end(),
				 theta,AngularBin::ReverseResolutionOrder);
    if (iter.first != iter.second) {
      --iter.second;
      sin2theta_max = iter.second->Sin2ThetaMax();
    }
  }

  return sin2theta_max;
}

ThetaIterator AngularCorrelation::Begin(uint32_t resolution) {
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      return theta_pair_begin_;
    } else {
      return thetabin_.begin();
    }
  } else {
    AngularBin theta;
    theta.SetResolution(resolution);
    ThetaPair iter = equal_range(theta_pixel_begin_, thetabin_.end(),
				 theta, AngularBin::ReverseResolutionOrder);
    return iter.first;
  }
}

ThetaIterator AngularCorrelation::End(uint32_t resolution) {
  if ((resolution < HPixResolution) ||
      (resolution > MaxPixelResolution) ||
      (resolution % 2 != 0)) {
    if (resolution == 0) {
      return theta_pair_end_;
    } else {
      return thetabin_.end();
    }
  } else {
    AngularBin theta;
    theta.SetResolution(resolution);
    ThetaPair iter = equal_range(theta_pixel_begin_,thetabin_.end(),
				 theta,AngularBin::ReverseResolutionOrder);
    return iter.second;
  }
}

ThetaIterator AngularCorrelation::Find(ThetaIterator begin,
				       ThetaIterator end,
				       double sin2theta) {
  ThetaIterator top = --end;
  ThetaIterator bottom = begin;
  ThetaIterator iter;

  if ((sin2theta < bottom->Sin2ThetaMin()) ||
      (sin2theta > top->Sin2ThetaMax())) {
    iter = ++end;
  } else {
    ++top;
    --bottom;
    while (top-bottom > 1) {
      iter = bottom + (top - bottom)/2;
      if (sin2theta < iter->Sin2ThetaMin()) {
        top = iter;
      } else {
        bottom = iter;
      }
    }
    iter = bottom;
  }

  return iter;
}

ThetaIterator AngularCorrelation::BinIterator(uint8_t bin_idx) {
  ThetaIterator iter = thetabin_.begin();

  for (uint8_t i=0;i<bin_idx;++i) {
    if (iter != thetabin_.end()) ++iter;
  }

  return iter;
}

uint32_t AngularCorrelation::NBins() {
  return thetabin_.size();
}

uint32_t AngularCorrelation::MinResolution() {
  return min_resolution_;
}

uint32_t AngularCorrelation::MaxResolution() {
  return max_resolution_;
}

} // end namespace Stomp
