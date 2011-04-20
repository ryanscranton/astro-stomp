#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

enum SystematicsField {
  Seeing,
  Extinction,
  Sky,
  Odds
};

namespace Stomp {
  class SysPixel : public Pixel {
  public:
    SysPixel() {
      SetResolution(0);

      SetPixnumFromXY(0, 0);
      SetWeight(0.0);

      sum_seeing_ = 0.0;
      sum_extinction_ = 0.0;
      sum_sky_ = 0.0;
      sum_odds_ = 0.0;
      n_obj_ = 0;
    };
    SysPixel(uint32_t resolution, uint32_t pixnum, double unmasked_fraction) {
      SetResolution(resolution);

      uint32_t tmp_y = pixnum/(Stomp::Nx0*Resolution());
      uint32_t tmp_x = pixnum - Stomp::Nx0*Resolution()*tmp_y;

      SetPixnumFromXY(tmp_x, tmp_y);
      SetWeight(unmasked_fraction);

      sum_seeing_ = 0.0;
      sum_extinction_ = 0.0;
      sum_sky_ = 0.0;
      sum_odds_ = 0.0;
      n_obj_ = 0;
    };
    SysPixel(uint32_t resolution, uint32_t pixnum, double unmasked_fraction,
	     double mean_seeing, double mean_extinction, double mean_sky,
	     double mean_odds, int n_objects) {
      SetResolution(resolution);

      uint32_t tmp_y = pixnum/(Stomp::Nx0*Resolution());
      uint32_t tmp_x = pixnum - Stomp::Nx0*Resolution()*tmp_y;

      SetPixnumFromXY(tmp_x, tmp_y);
      SetWeight(unmasked_fraction);

      sum_seeing_ = mean_seeing*n_objects;
      sum_extinction_ = mean_extinction*n_objects;
      sum_sky_ = mean_sky*n_objects;
      sum_odds_ = mean_odds*n_objects;
      n_obj_ = n_objects;
    };
    ~SysPixel() {
      SetResolution(0);

      SetPixnumFromXY(0, 0);
      SetWeight(0.0);

      sum_seeing_ = 0.0;
      sum_extinction_ = 0.0;
      sum_sky_ = 0.0;
      sum_odds_ = 0.0;
      n_obj_ = 0;
    }
    inline void AddSysParameters(double seeing, double extinction, 
				 double sky, double odds) {
      sum_seeing_ += seeing;
      sum_extinction_ += extinction;
      sum_sky_ += sky;
      sum_odds_ += odds;
      n_obj_++;
    };
    inline void SetSysParameters(double seeing, double extinction, 
				 double sky, double odds, int n_obj) {
      sum_seeing_ = n_obj*seeing;
      sum_extinction_ = n_obj*extinction;
      sum_sky_ = n_obj*sky;
      sum_odds_ = n_obj*odds;
      n_obj_ = n_obj;
    };
    inline double MeanSeeing() {
      return (n_obj_ > 0 ? sum_seeing_/n_obj_ : -1.0);
    };
    inline double MeanExtinction() {
      return (n_obj_ > 0 ? sum_extinction_/n_obj_ : -1.0);
    };
    inline double MeanSky() {
      return (n_obj_ > 0 ? sum_sky_/n_obj_ : -1.0);
    };
    inline double MeanOdds() {
      return (n_obj_ > 0 ? sum_odds_/n_obj_ : -1.0);
    };
    inline double UnmaskedFraction() {
      return Weight();
    };
    inline int NObjects() {
      return n_obj_;
    };

  private:
    double sum_seeing_, sum_extinction_, sum_sky_, sum_odds_;
    int n_obj_;
  };
} // end namespace Stomp

typedef std::map<const uint32_t, Stomp::SysPixel> SysDict;
typedef SysDict::iterator SysDictIterator;

void LoadMapData(Stomp::Map *&stomp_map, SysDict &sysmap,
		 Stomp::PixelVector &map_pix);
void LoadRegions(Stomp::RegionBoundVector &region_vector);
void LoadGalaxyData(Stomp::WAngularVector &galaxy, 
		    Stomp::Map *stomp_map);
void GalaxyAutoCorrelation(Stomp::WAngularVector &galaxy, 
			   Stomp::Map *stomp_map,
			   std::string const &ouput_suffix);
void SystematicsCorrelation(Stomp::WAngularVector &galaxy,
			    Stomp::WAngularVector &random,
			    SystematicsField const &sys_field,
			    Stomp::Map *stomp_map,
			    SysDict &sysmap,
			    Stomp::PixelVector &map_pix,
			    std::string const &ouput_suffix);

// Define our command-line flags.
DEFINE_string(map_file, "",
              "Name of the ASCII file containing the StompMap geometry");
DEFINE_string(sysmap_file, "",
              "Name of sysmap containing the scalar values to sample.");
DEFINE_string(galaxy_file, "",
              "Name of the ASCII file containing the first galaxies");
DEFINE_string(region_file, "",
              "Name of the ASCII file containing region definitions");
DEFINE_bool(galaxy_radec, true, "Galaxy coordinates are in RA-DEC");
DEFINE_bool(use_only_pairs, false, "Use only the stomp pair based estimator");
DEFINE_string(output_tag, "test",
              "Tag for output file: Wtheta_OUTPUT_TAG");
DEFINE_string(output_map, "",
              "Name of resulting Stomp map from excluding systematics.");
DEFINE_double(theta_min, 0.001, "Minimum angular scale (in degrees)");
DEFINE_double(theta_max, 10.0, "Maximum angular scale (in degrees)");
DEFINE_double(mag_min, 18.0, "Minimum acceptable galaxy magnitude");
DEFINE_double(mag_max, 28.0, "Maximum acceptable galaxy magnitude");
DEFINE_double(prob_min, 0.2, "Minimum acceptable galaxy likelihood");
DEFINE_double(prob_max, 1.00001, "Maximum acceptable galaxy likelihood");
DEFINE_double(see_min, 0.5, "Minimum allowable seeing value");
DEFINE_double(see_max, 2.0, "Maximum allowable seeing value");
DEFINE_double(ext_min, 0.0, "Minimum allowable extinction value");
DEFINE_double(ext_max, 1.0, "Maximum allowable extinction value");
DEFINE_double(sky_min, 0.0, "Minimum allowable sky brightness value");
DEFINE_double(sky_max, 24.0, "Maximum allowable sky brightness value");
DEFINE_double(odds_min, 0.0, "Minimum allowable odds brightness value");
DEFINE_double(odds_max, 1.01, "Maximum allowable odds brightness value");
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
DEFINE_int32(n_bins_per_decade, 5, "Number of angular bins per decade.");
DEFINE_int32(n_random, 1,
	     "Integer number of random points per galaxy to use.");
DEFINE_int32(n_jackknife, 0,
	     "Number of jack-knife samples to use. Defaults to 2*angular bins");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Map file is missing weight column.");
DEFINE_bool(coordinates_only, false, "Galaxy files only contain coordinates.");
DEFINE_int32(maximum_resolution, -1,
	     "Maximum resolution to use for pixel-based estimator");

int main(int argc, char **argv) {

  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --map_file=<StompMap ASCII>";
  usage += " --sample_map=<Weighted StompMap to sample>";
  usage += " --galaxy_file=<list of galaxy input files>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_map_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_sysmap_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_galaxy_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  std::ostringstream output_suffix;
  output_suffix << "_see" << FLAGS_see_min << "-" << FLAGS_see_max;
  output_suffix << "_ext" << FLAGS_ext_min << "-" << FLAGS_ext_max;
  output_suffix << "_sky" << FLAGS_sky_min << "-" << FLAGS_sky_max;
  output_suffix << "_odds" << FLAGS_odds_min << "-" << FLAGS_odds_max;
  output_suffix << "_combined";

  // Read In the stomp map and systematics maps and intersect the area.
  Stomp::Map* stomp_map;
  SysDict sysmap;
  Stomp::PixelVector map_pix;
  LoadMapData(stomp_map, sysmap, map_pix);
  if (!FLAGS_output_map.empty()) {
    stomp_map->Write(FLAGS_output_map);
  }

  if (!FLAGS_region_file.empty()) {
    Stomp::RegionBoundVector region_vector;
    LoadRegions(region_vector);
    std::cout << "Regionating with " << region_vector.size() 
	      << " regions...\n";
    FLAGS_n_jackknife = 
      stomp_map->InitializeRegions(region_vector, region_vector.size());
  }

  //Create random positions to sample from systematics map.
  std::cout << "Requested " << 10*map_pix.size() << " random points...\n";
  Stomp::WAngularVector random_points;
  //while (random_points.size() < 10*map_pix.size()) {
  while (random_points.size() < 10*map_pix.size()) {
    Stomp::WeightedAngularCoordinate tmp_ang(0.0, 0.0, 0.0);
    stomp_map->GenerateSingleRandomPoint(tmp_ang, false, false);
    random_points.push_back(tmp_ang);
  }
  std::cout << "Random Points Returned " << random_points.size() 
	    << " random points...\n";

  //Read in galaxy data using the intersected map
  Stomp::WAngularVector galaxy;
  LoadGalaxyData(galaxy, stomp_map);


  //Run auto-correlation on galaxies.
  if (FLAGS_galaxy_galaxy)
    GalaxyAutoCorrelation(galaxy, stomp_map,output_suffix.str());

  if (FLAGS_galaxy_seeing)
    SystematicsCorrelation(galaxy, random_points, Seeing, stomp_map,
			   sysmap, map_pix, output_suffix.str());

  if (FLAGS_galaxy_extinction)
    SystematicsCorrelation(galaxy, random_points, Extinction, stomp_map,
			   sysmap, map_pix, output_suffix.str());

  if (FLAGS_galaxy_sky)
    SystematicsCorrelation(galaxy, random_points, Sky, stomp_map,
			   sysmap, map_pix, output_suffix.str());

  if (FLAGS_galaxy_odds)
    SystematicsCorrelation(galaxy, random_points, Odds, stomp_map,
			   sysmap, map_pix, output_suffix.str());

  delete stomp_map;

  return 0;
}

void LoadMapData(Stomp::Map *&stomp_map, SysDict &sysmap,
		 Stomp::PixelVector & map_pix) {

  if (FLAGS_single_index) {
    if (FLAGS_no_weight) {
      stomp_map = new Stomp::Map(FLAGS_map_file, false, false);
    } else {
      stomp_map = new Stomp::Map(FLAGS_map_file, false);
    }
  } else {
    if (FLAGS_no_weight) {
      stomp_map = new Stomp::Map(FLAGS_map_file, true, false);
    } else {
      stomp_map = new Stomp::Map(FLAGS_map_file);
    }
  }
  std::cout << "Read map from " << FLAGS_map_file << "; initial area: " <<
    stomp_map->Area() << " sq. deg...\n";

  std::ifstream systematics_file(FLAGS_sysmap_file.c_str());
  std::cout << "Reading systematics map from " << FLAGS_sysmap_file << "...\n";

  double total_area = 0.0;
  unsigned long n_pixel = 0, n_keep = 0;
  unsigned long n_seeing = 0, n_extinction = 0, n_sky = 0, n_odds = 0, 
    n_blank = 0;
  uint32_t x, y, hpixnum, superpixnum, resolution, n_objects;
  double unmasked, seeing, extinction, sky, odds;

  while (!systematics_file.eof()) {
    systematics_file >> hpixnum >> superpixnum >> resolution 
		     >> unmasked >> seeing >> extinction 
		     >> sky >> odds >> n_objects;

    if (!systematics_file.eof()) {
      n_pixel++;
      if ((seeing >= FLAGS_see_min) &&
	  (seeing <= FLAGS_see_max) &&
	  (extinction >= FLAGS_ext_min) &&
	  (extinction <= FLAGS_ext_max) &&
	  (sky >= FLAGS_sky_min) &&
	  (sky <= FLAGS_sky_max) &&
	  (odds >= FLAGS_odds_min) &&
	  (odds <= FLAGS_odds_max) && 
	  (n_objects > 0)) {
	Stomp::Pixel::HPix2XY(resolution, hpixnum, superpixnum, x, y);
	Stomp::Pixel tmp_pix(x, y, resolution, 1.0);
	total_area += unmasked*tmp_pix.Area();
	map_pix.push_back(tmp_pix);
	Stomp::SysPixel sys_pix(resolution, tmp_pix.Pixnum(), unmasked,
				seeing, extinction, sky,
				odds, n_objects);
	sysmap[tmp_pix.Pixnum()] =  sys_pix;
	n_keep++;
      }
      else {
	if ((seeing < FLAGS_see_min) || (seeing > FLAGS_see_max))
	  n_seeing++;
	if ((extinction < FLAGS_ext_min) ||
	    (extinction > FLAGS_ext_max)) n_extinction++;
	if ((sky >= FLAGS_sky_min) || (sky <= FLAGS_sky_max)) n_sky++;
	if ((odds >= FLAGS_odds_min) || (odds <= FLAGS_odds_max))
	if (n_objects == 0) n_blank++;
      }
    }
  }
  systematics_file.close();
  std::cout << "\tKept " << n_keep << "/" << n_pixel <<
    " systematics pixels; " << total_area << " sq. degrees...\n";
  std::cout << "\t\t" << n_seeing << " failed seeing cut.\n";
  std::cout << "\t\t" << n_extinction << " failed extinction cut.\n";
  std::cout << "\t\t" << n_sky << " failed sky_brightness cut.\n";
  std::cout << "\t\t" << n_odds << " failed odds cut.\n";
  std::cout << "\t\t" << n_blank << " were empty.\n";

  std::cout << "Excluding Pixels from map...\n";
  stomp_map->IntersectMap(map_pix);

  std::cout << "Final Map Area: " << stomp_map->Area() << " sq. deg.\n";
}

void LoadRegions(Stomp::RegionBoundVector &region_vector) {

  int type;
  double radius;
  double ra[4], dec[4];
  unsigned long idx;
  double weight;
  std::ifstream region_file(FLAGS_region_file.c_str());
  
  while (!region_file.eof()) {
    region_file >> idx >> type >> radius >> weight 
		>> ra[0] >> dec[0] >> ra[1] >> dec[1] 
		>> ra[2] >> dec[2] >> ra[3] >> dec[3];
    
    if (!region_file.eof()) {
      
      Stomp::AngularVector ang_vector;
      for (int i=0;i<4;i++)
	ang_vector.push_back(Stomp::AngularCoordinate
			     (ra[i],dec[i],					
			      Stomp::AngularCoordinate::Equatorial));
      
      Stomp::PolygonBound* poly_bound = new Stomp::PolygonBound(ang_vector);
      Stomp::RegionBound region_bound(poly_bound);
      
      region_vector.push_back(region_bound);
    }
  }
}

void LoadGalaxyData(Stomp::WAngularVector &galaxy, 
		    Stomp::Map *stomp_map) {

  std::cout << "Parsing " << FLAGS_galaxy_file << " file...\n";

  Stomp::AngularCoordinate::Sphere galaxy_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_galaxy_radec) galaxy_sphere = Stomp::AngularCoordinate::Equatorial;
  std::ifstream galaxy_file(FLAGS_galaxy_file.c_str());
  double theta, phi, prob, mag;
  long int n_galaxy = 0;
  
  prob = 1.0;
  mag = 0.5*(FLAGS_mag_max + FLAGS_mag_min);
  
  while (!galaxy_file.eof()) {
    if (FLAGS_coordinates_only) {
      galaxy_file >> theta >> phi;
    } else {
      galaxy_file >> theta >> phi >> prob >> mag;
    }
    
    if (!galaxy_file.eof()) {
      Stomp::WeightedAngularCoordinate tmp_ang(theta, phi,
					       prob, galaxy_sphere);
      
      if ((mag >= FLAGS_mag_min) && (mag <= FLAGS_mag_max) &&
	  (stomp_map->FindLocation(tmp_ang, prob))) galaxy.push_back(tmp_ang);
      n_galaxy++;
    }
  }
  galaxy_file.close();
  std::cout << "Read " << n_galaxy << " galaxies; kept " <<
    galaxy.size() << "\n";
  n_galaxy = galaxy.size();
}

void GalaxyAutoCorrelation(Stomp::WAngularVector &galaxy, 
			   Stomp::Map *stomp_map,
			   std::string const &output_suffix) {
  
  Stomp::AngularCorrelation wtheta(FLAGS_theta_min, FLAGS_theta_max,
				   FLAGS_n_bins_per_decade);

  if (FLAGS_maximum_resolution > 0) {
    std::cout << "Setting maximum resolution to " <<
      static_cast<uint16_t>(FLAGS_maximum_resolution) << "...\n";
    wtheta.SetMaxResolution(static_cast<uint16_t>(FLAGS_maximum_resolution));
  }
  if (FLAGS_use_only_pairs) {
    wtheta.UseOnlyPairs();
  }
  wtheta.FindAutoCorrelationWithRegions(*stomp_map, galaxy,
					static_cast<uint8_t>(FLAGS_n_random),
					static_cast<uint16_t>(FLAGS_n_jackknife)
					);
  std::string wtheta_file_name = "Wtheta_gal-gal_" + FLAGS_output_tag +
    output_suffix;
  wtheta.Write(wtheta_file_name);
}

void SystematicsCorrelation(Stomp::WAngularVector &galaxy,
			    Stomp::WAngularVector &random_points,
			    SystematicsField const &sys_field,
			    Stomp::Map *stomp_map,
			    SysDict &sysmap,
			    Stomp::PixelVector & map_pix,
			    std::string const &output_suffix) {

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
    std::cout << "Reading in odds values from " <<
      FLAGS_sysmap_file << "...\n";
    break;
  }

  for (Stomp::PixelIterator iter=map_pix.begin();
       iter!=map_pix.end();++iter) {
    SysDictIterator sys_iter = sysmap.find(iter->Pixnum());
    switch (sys_field) {
    case Seeing:
      iter->SetWeight(sys_iter->second.MeanSeeing());
      break;
    case Extinction:
      iter->SetWeight(sys_iter->second.MeanExtinction());
      break;
    case Sky:
      iter->SetWeight(sys_iter->second.MeanSky());
      break;
    case Odds:
      iter->SetWeight(sys_iter->second.MeanOdds());
      break;
    }
  }

  std::cout << "Setting random point weights...\n";
  double weight;
  Stomp::Map sys_sample_map(map_pix);
  for (Stomp::WAngularIterator iter = random_points.begin(); 
       iter != random_points.end(); iter++) {
    weight = sys_sample_map.FindLocationWeight(*iter);
    iter->SetWeight(weight);
  }

  Stomp::AngularCorrelation wtheta(FLAGS_theta_min, FLAGS_theta_max,
				   FLAGS_n_bins_per_decade);

  if (FLAGS_maximum_resolution > 0) {
    std::cout << "Setting maximum resolution to " <<
      static_cast<uint16_t>(FLAGS_maximum_resolution) << "...\n";
    wtheta.SetMaxResolution(static_cast<uint16_t>(FLAGS_maximum_resolution));
  }
  if (FLAGS_use_only_pairs) {
    wtheta.UseOnlyPairs();
  }
  wtheta.FindCrossCorrelationWithRegions(*stomp_map, galaxy, random_points,
					static_cast<uint8_t>(FLAGS_n_random),
					static_cast<uint16_t>(FLAGS_n_jackknife)
					);

  std::string sys_type;
  switch (sys_field) {
  case Seeing:
    sys_type.append("seeing");
    break;
  case Extinction:
    sys_type.append("ext");
    break;
  case Sky:
    sys_type.append("sky");
    break;
  case Odds:
    sys_type.append("odds");
    break;
  }
  std::string wtheta_file_name = "Wtheta_gal-" + sys_type + "_" + 
    FLAGS_output_tag + output_suffix;
  wtheta.Write(wtheta_file_name);
}
