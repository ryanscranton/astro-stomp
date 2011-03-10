#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>


// Define our command-line flags.
DEFINE_string(gal_file, "",
              "ASCII file containing galaxies with ra,dec,weight,mag");
DEFINE_string(gal_map_output, "",
	      "Outputs galaxy density map");
DEFINE_string(sysmap_file, "",
              "ASCII file containing uniform resolution sysmap");
DEFINE_string(sysmap_cut_file, "",
	      "ASCII file containing sysmap file to cut on");
DEFINE_string(output_tag, "test",
              "Tag for output file: MeanWthetaSys_OUTPUT_TAG_combined");
DEFINE_string(input_bounds, "",
	      "ASCII file containing geometric bound information");
DEFINE_double(seeing_min, 0.5, "Minimum allowable seeing value");
DEFINE_double(seeing_max, 2.0, "Maximum allowable seeing value");
DEFINE_double(extinction_min, 0.0, "Minimum allowable extinction value");
DEFINE_double(extinction_max, 1.0, "Maximum allowable extinction value");
DEFINE_double(sky_min, 0.0, "Minimum allowable sky brightness value");
DEFINE_double(sky_max, 24.0, "Maximum allowable sky brightness value");
DEFINE_double(odds_min, 0.0, "Minimum allowable sky brightness value");
DEFINE_double(odds_max, 1.01, "Maximum allowable sky brightness value");
DEFINE_bool(galaxy_galaxy, false,
	    "Auto-correlate galaxy density.");
DEFINE_bool(galaxy_seeing, false,
	    "Cross-correlate galaxy density and seeing variations.");
DEFINE_bool(galaxy_extinction, false,
	    "Cross-correlate galaxy density and extinction variations.");
DEFINE_bool(galaxy_sky, false,
	    "Cross-correlate galaxy density and sky brightness variations.");
DEFINE_bool(galaxy_odds, false,
	    "Cross-correlate galaxy density and average odds variations.");
DEFINE_double(theta_min, 0.09, "Minimum angular scale (in degrees)");
DEFINE_double(theta_max, 2.0, "Maximum angular scale (in degrees)");
DEFINE_double(galmap_resolution, 2048, "Resolution of Galaxy Map");
DEFINE_int32(n_bins_per_decade, 5, "Number of angular bins per decade.");

typedef std::map<const uint32_t, int> KeepDict;
typedef KeepDict::iterator KeepDictIterator;
Stomp::ScalarVector galaxy_pix;
Stomp::RegionBoundVector region_vector;

enum SystematicsField {
  Seeing,
  Extinction,
  Sky,
  Odds
};

int main(int argc, char **argv) {
  void FindAllowedArea(const char* sysmap_name, KeepDict& keep_map);
  Stomp::ScalarMap* LoadGalaxyData(KeepDict& keep_map,
				   int region_resolution);
  void AutoCorrelateGalaxies(Stomp::ScalarMap* galaxy_map,
			     Stomp::AngularCorrelation& wtheta);
  void CrossCorrelateSystematicsField(Stomp::ScalarMap* galaxy_map,
				      SystematicsField sys_field,
				      Stomp::AngularCorrelation& wtheta,
				      KeepDict& keep_map);

  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --gal_file=<SysMap ASCII>";
  usage += " --sysmap_file=<SysMap ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_gal_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_sysmap_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  std::cout << "Systematics bounds:\n";
  std::cout <<
    "\t" << FLAGS_seeing_min << " < seeing < " << FLAGS_seeing_max << "\n";
  std::cout <<
    "\t" << FLAGS_extinction_min << " < extinction < " <<
    FLAGS_extinction_max << "\n";
  std::cout <<
    "\t" << FLAGS_sky_min << " < sky brightness < " << FLAGS_sky_max << "\n";
  std::cout <<
    "\t" << FLAGS_odds_min << " < odds brightness < " << FLAGS_odds_max << "\n";

  // Now we read in our sysmap.  We only want to keep those pixels that
  // meet our selection criteria, so we filter on those as we read through
  // the file.
  
  if (!FLAGS_input_bounds.empty()) {
    int type;
    double radius;
    double ra[4], dec[4];
    unsigned long idx;
    double weight;
    std::ifstream region_file(FLAGS_input_bounds.c_str());
    
    while (!region_file.eof()) {
      region_file >> idx >> type >> radius >> weight >> ra[0] >> dec[0] >> ra[1] >> dec[1] >> ra[2] >> dec[2] >> ra[3] >> dec[3];
      
      if (!region_file.eof()) {
	
	Stomp::AngularVector ang_vector;
	for (int i=0;i<4;i++)
	  ang_vector.push_back(Stomp::AngularCoordinate(ra[i],dec[i],Stomp::AngularCoordinate::Equatorial));
	
	Stomp::PolygonBound* poly_bound = new Stomp::PolygonBound(ang_vector);
	Stomp::RegionBound region_bound(poly_bound);
	
	region_vector.push_back(region_bound);
      }
    }
  }
  KeepDict keep_map;
  if (FLAGS_sysmap_cut_file.empty()) {
    FindAllowedArea(FLAGS_sysmap_file.c_str(), keep_map);
  } else {
    FindAllowedArea(FLAGS_sysmap_cut_file.c_str(), keep_map);
  }

  // We need some constraints from our angular binning, so we'll make a
  // single AngularCorrelation object at the outset.  If we do the galaxy
  // auto-correlation, then we'll use the object for that calculation.
  Stomp::AngularCorrelation wtheta(FLAGS_theta_min, FLAGS_theta_max,
				   FLAGS_n_bins_per_decade);
 

  // Now, read in the galaxy map.  We only want to keep those pixels that match
  // the ones that passed through our systematics cuts.
  Stomp::ScalarMap* galaxy_map =
    LoadGalaxyData(keep_map, wtheta.MinResolution());

  // Set up the file name suffix that specifies our systematics cuts.
  std::ostringstream output_suffix;
  output_suffix << "_see" << FLAGS_seeing_min << "-" << FLAGS_seeing_max;
  output_suffix << "_ext" << FLAGS_extinction_min << "-" <<
    FLAGS_extinction_max;
  output_suffix << "_sky" << FLAGS_sky_min << "-" << FLAGS_sky_max;
  output_suffix << "_odds" << FLAGS_odds_min << "-" << FLAGS_odds_max;
  output_suffix << "_combined";

  // First, handle the galaxy auto-correlation case
  if (FLAGS_galaxy_galaxy) {
    AutoCorrelateGalaxies(galaxy_map, wtheta);

    // And write out the results...
    std::string output_file_name =
      "MeanWthetaSys_gal-gal_" + FLAGS_output_tag + output_suffix.str();
    std::cout << "Writing results to " << output_file_name << "...\n";

    wtheta.Write(output_file_name);
  }

  // Galaxy-Seeing
  if (FLAGS_galaxy_seeing) {
    Stomp::AngularCorrelation xcorr_seeing(FLAGS_theta_min, FLAGS_theta_max,
					   FLAGS_n_bins_per_decade);

    CrossCorrelateSystematicsField(galaxy_map, Seeing,
				   xcorr_seeing, keep_map);

    std::string output_file_name =
      "MeanWthetaSys_gal-see_" + FLAGS_output_tag + output_suffix.str();
    std::cout << "Writing results to " << output_file_name << "...\n";

    xcorr_seeing.Write(output_file_name);
  }

  // Galaxy-Extinction
  if (FLAGS_galaxy_extinction) {
    Stomp::AngularCorrelation xcorr_extinction(FLAGS_theta_min, FLAGS_theta_max,
					       FLAGS_n_bins_per_decade);

    CrossCorrelateSystematicsField(galaxy_map, Extinction,
				   xcorr_extinction, keep_map);

    std::string output_file_name =
      "MeanWthetaSys_gal-ext_" + FLAGS_output_tag + output_suffix.str();
    std::cout << "Writing results to " << output_file_name << "...\n";

    xcorr_extinction.Write(output_file_name);
  }

  // Galaxy-Sky Brightness
  if (FLAGS_galaxy_sky) {
    Stomp::AngularCorrelation xcorr_sky(FLAGS_theta_min, FLAGS_theta_max,
					   FLAGS_n_bins_per_decade);

    CrossCorrelateSystematicsField(galaxy_map, Sky,
				   xcorr_sky, keep_map);

    std::string output_file_name =
      "MeanWthetaSys_gal-sky_" + FLAGS_output_tag + output_suffix.str();
    std::cout << "Writing results to " << output_file_name << "...\n";

    xcorr_sky.Write(output_file_name);
  }

  // Galaxy-Odds
  if (FLAGS_galaxy_odds) {
    Stomp::AngularCorrelation xcorr_odds(FLAGS_theta_min, FLAGS_theta_max,
					   FLAGS_n_bins_per_decade);

    CrossCorrelateSystematicsField(galaxy_map, Odds,
				   xcorr_odds, keep_map);

    std::string output_file_name =
      "MeanWthetaSys_gal-odds_" + FLAGS_output_tag + output_suffix.str();
    std::cout << "Writing results to " << output_file_name << "...\n";

    xcorr_odds.Write(output_file_name);
  }

  std::cout << "Done.\n";

  return 0;
}

void FindAllowedArea(const char* sysmap_name, KeepDict& keep_map) {
  std::ifstream sysmap_file(sysmap_name);
  double unmasked, seeing, extinction, sky, odds;
  uint32_t x, y, hpixnum, superpixnum, resolution, n_objects;
  double total_area = 0.0;
  

  std::cout << "Reading systematics map from " << FLAGS_sysmap_file << "...\n";
  unsigned long n_pixel = 0, n_keep = 0;
  unsigned long n_seeing = 0, n_extinction = 0, n_sky = 0, n_odds = 0, 
    n_blank = 0;
  while (!sysmap_file.eof()) {
    sysmap_file >> hpixnum >> superpixnum >> resolution >> unmasked >> seeing >>
      extinction >> sky >> odds >> n_objects;

    if (!sysmap_file.eof()) {
      n_pixel++;
      if ((seeing >= FLAGS_seeing_min) &&
	  (seeing <= FLAGS_seeing_max) &&
	  (extinction >= FLAGS_extinction_min) &&
	  (extinction <= FLAGS_extinction_max) &&
	  (sky >= FLAGS_sky_min) &&
	  (sky <= FLAGS_sky_max) &&
	  (odds >= FLAGS_odds_min) &&
	  (odds <= FLAGS_odds_max)
	  && (n_objects > 0)) {
	bool keep_pixel = false;
	Stomp::Pixel::HPix2XY(resolution, hpixnum, superpixnum, x, y);
	Stomp::Pixel tmp_pix(x, y, resolution, 1.0);
	for (uint32_t bound_iter=0;bound_iter<region_vector.size();bound_iter++) {
	  if (region_vector[bound_iter].CheckPixel(tmp_pix)) {
	    if (keep_pixel) {
	      std::cout << "Syspixel found in more than one bound!\n";
	      keep_pixel = false;
	      break;
	    }
	    keep_pixel = true;  
	  }
	}
	//Stomp::ScalarPixel tmp_pix(x, y, resolution, unmasked, 0, 0);
	if (keep_pixel) {
	  total_area += unmasked*tmp_pix.Area();
	  keep_map[tmp_pix.Pixnum()] = 1;
	  n_keep++;
	}
	//galaxy_pix.push_back(tmp_pix);
      } else {
	if ((seeing < FLAGS_seeing_min) || (seeing > FLAGS_seeing_max))
	  n_seeing++;
	if ((extinction < FLAGS_extinction_min) ||
	    (extinction > FLAGS_extinction_max)) n_extinction++;
	if ((sky >= FLAGS_sky_min) || (sky <= FLAGS_sky_max)) n_sky++;
	if ((odds >= FLAGS_odds_min) || (odds <= FLAGS_odds_max))
	if (n_objects == 0) n_blank++;
      }
    }
  }
  sysmap_file.close();
  std::cout << "\tKept " << n_keep << "/" << n_pixel <<
    " systematics pixels; " << total_area << " sq. degrees...\n";
  std::cout << "\t\t" << n_seeing << " failed seeing cut.\n";
  std::cout << "\t\t" << n_extinction << " failed extinction cut.\n";
  std::cout << "\t\t" << n_sky << " failed sky_brightness cut.\n";
  std::cout << "\t\t" << n_odds << " failed odds cut.\n";
  std::cout << "\t\t" << n_blank << " were empty.\n";
}

Stomp::ScalarMap* LoadGalaxyData(KeepDict& keep_map, int region_resolution) {
  std::ifstream sysmap_file(FLAGS_sysmap_file.c_str());
  double unmasked, seeing, extinction, sky, odds;
  uint32_t x, y, hpixnum, superpixnum, resolution, n_objects;
  double total_area = 0.0;
  Stomp::ScalarVector galaxy_pix;
  KeepDictIterator keep_iter;

  while (!sysmap_file.eof()) {
    sysmap_file >> hpixnum >> superpixnum >> resolution >> unmasked >> seeing >>
      extinction >> sky >> odds >> n_objects;

    if (!sysmap_file.eof()) { 
      Stomp::Pixel::HPix2XY(resolution, hpixnum, superpixnum, x, y);
      //Stomp::Pixel tmp_pix(x, y, resolution, 1.0);
      Stomp::ScalarPixel tmp_pix(x, y, resolution, unmasked, 0, 0);
      keep_iter = keep_map.find(tmp_pix.Pixnum());
      if (keep_iter != keep_map.end()) {
	total_area += unmasked*tmp_pix.Area();
	galaxy_pix.push_back(tmp_pix);
      }
    }
  }

  std::ifstream gal_file(FLAGS_gal_file.c_str());
  double ra, dec, weight, mag;
  double n_galaxies = 0;

  std::cout << "Reading galaxies from " << FLAGS_gal_file << "...\n";
  Stomp::ScalarMap* galaxy_map =
    new Stomp::ScalarMap(galaxy_pix, Stomp::ScalarMap::DensityField);
  while (!gal_file.eof()) {
    gal_file >> ra >> dec >> weight >> mag;
    
    Stomp::AngularCoordinate ang(ra, dec, Stomp::AngularCoordinate::Equatorial);
    galaxy_map->AddToMap(ang,1);
    //Stomp::Pixel pix(ang, FLAGS_galmap_resolution, 1);
    //keep_iter = keep_map.find(pix.Pixnum());
    //if (keep_iter != keep_map.end()) {
    //  galaxy_map->AddToMap(pix);
    n_galaxies++;
      //}
  }
  gal_file.close();
  if (!FLAGS_gal_map_output.empty()) {
    std::ofstream output_file(FLAGS_gal_map_output.c_str());
    for (Stomp::ScalarIterator iter=galaxy_pix.begin();
	 iter!=galaxy_pix.end();++iter) {
      output_file << 
	iter->HPixnum() << " " <<
	iter->Superpixnum() <<  " " <<
	iter->Resolution() << " " <<
	galaxy_map->FindUnmaskedFraction(*iter) << " " <<
	galaxy_map->FindIntensity(*iter) << " " <<
	galaxy_map->FindDensity(*iter)/(60.0*60.0) << " " <<
	galaxy_map->FindPointDensity(*iter)/(60.0*60.0) << "\n";
    }
    output_file.close();
  }
  galaxy_pix.clear();

  std::cout << "\tFound " << n_galaxies << " galaxies; " <<
    total_area << " sq. degrees...\n";

  // Convert raw galaxy pixels into a ScalarMap.

  return galaxy_map;
}

void AutoCorrelateGalaxies(Stomp::ScalarMap* galaxy_map,
			   Stomp::AngularCorrelation& wtheta) {
  // Since we're using jack knife samples, we need to initialize that
  // functionality in the AngularCorrelation object.

  // And set up our AngularCorrelation resolution limits around the limits of
  // the galaxy map.
  //wtheta.SetMinResolution(galaxy_map->RegionResolution());
  wtheta.SetMaxResolution(galaxy_map->Resolution());

  // Now we do the requested auto-correlation.
  int max_resolution = galaxy_map->Resolution();
  int min_resolution = wtheta.MinResolution();

  std::cout << "Starting galaxy auto-correlation...\n";
  for (Stomp::ThetaIterator iter=wtheta.Begin(max_resolution);
       iter!=wtheta.End(max_resolution);++iter)
    galaxy_map->AutoCorrelate(iter);

  for (int resolution=max_resolution/2;
       resolution>=min_resolution;resolution/=2) {
    Stomp::ScalarMap* galaxy_sub_map =
      new Stomp::ScalarMap(*galaxy_map, resolution);
    for (Stomp::ThetaIterator iter=wtheta.Begin(resolution);
	 iter!=wtheta.End(resolution);++iter) {
      galaxy_sub_map->AutoCorrelate(iter);
    }
    delete galaxy_sub_map;
  }
}

void CrossCorrelateSystematicsField(Stomp::ScalarMap* galaxy_map,
				    SystematicsField sys_field,
				    Stomp::AngularCorrelation& wtheta,
				    KeepDict& keep_map) {

  //wtheta.SetMinResolution(galaxy_map->RegionResolution());
  wtheta.SetMaxResolution(galaxy_map->Resolution());

  //Scan through the systematics file again and load up the appropriate values.
  Stomp::ScalarVector sys_pix;
  std::ifstream sysmap_file(FLAGS_sysmap_file.c_str());
  double unmasked, seeing, extinction, sky, odds;
  uint32_t x, y, hpixnum, superpixnum, resolution, n_objects;
  KeepDictIterator keep_iter;

  switch (sys_field) {
  case Seeing:
    std::cout << "Reading in seeing values from " <<
      FLAGS_sysmap_file << "...\n";
    break;
  case Extinction:
    std::cout << "Reading in extinction values from " <<
      FLAGS_sysmap_file << "...\n";
    break;
  case Sky:
    std::cout << "Reading in sky brightness values from " <<
      FLAGS_sysmap_file << "...\n";
    break;
  case Odds:
    std::cout << "Reading in sky brightness values from " <<
      FLAGS_sysmap_file << "...\n";
    break;
  }

  // Parse the systematics file and pull out the requested field
  while (!sysmap_file.eof()) {
    sysmap_file >> hpixnum >> superpixnum >> resolution >> unmasked >>
      seeing >> extinction >> sky >> odds >> n_objects;

    if (!sysmap_file.eof()) { 
      Stomp::Pixel::HPix2XY(resolution, hpixnum, superpixnum, x, y);
      //Stomp::Pixel tmp_pix(x, y, resolution, 1.0);
      Stomp::ScalarPixel pix(x, y, resolution, unmasked, 0, 0);
      keep_iter = keep_map.find(pix.Pixnum());
      if (keep_iter != keep_map.end()) {
	switch (sys_field) {
	case Seeing:
	  pix.SetIntensity(seeing);
	  break;
	case Extinction:
	  pix.SetIntensity(extinction);
	  break;
	case Sky:
	  pix.SetIntensity(sky);
	  break;
	case Odds:
	  pix.SetIntensity(odds);
	  break;
	}
	sys_pix.push_back(pix);
      }
    }
  }
  sysmap_file.close();

  // Set up our systematics ScalarMap
  Stomp::ScalarMap* sys_map =
    new Stomp::ScalarMap(sys_pix, Stomp::ScalarMap::ScalarField);
  sys_pix.clear();

  // Now we do the requested cross-correlation.
  int max_resolution = galaxy_map->Resolution();
  int min_resolution = wtheta.MinResolution();

  switch (sys_field) {
  case Seeing:
    std::cout << "Cross-correlating galaxies with seeing...\n";
    break;
  case Extinction:
    std::cout << "Cross-correlating galaxies with extinction...\n";
    break;
  case Sky:
    std::cout << "Cross-correlating galaxies with sky brightness...\n";
    break;
  case Odds:
    std::cout << "Cross-correlating galaxies with odds brightness...\n";
    break;
  }
  for (Stomp::ThetaIterator iter=wtheta.Begin(max_resolution);
       iter!=wtheta.End(max_resolution);++iter) {
    galaxy_map->CrossCorrelate(*sys_map, iter);
  }

  std::cout << "Completed on galaxy map\n";

  for (uint32_t resolution=max_resolution/2;
       resolution>=min_resolution;resolution/=2) {
    std::cout << "Current Resolution: " << resolution << "\n";
    Stomp::ScalarMap* galaxy_sub_map =
      new Stomp::ScalarMap(*galaxy_map, resolution);
    
    Stomp::ScalarMap* sys_sub_map = 
      new Stomp::ScalarMap(*sys_map, resolution);
    
    for (Stomp::ThetaIterator iter=wtheta.Begin(resolution);
	 iter!=wtheta.End(resolution);++iter) {
      galaxy_sub_map->CrossCorrelate(*sys_sub_map, iter);
    }
    delete sys_sub_map;
    delete galaxy_sub_map;
  }
}

