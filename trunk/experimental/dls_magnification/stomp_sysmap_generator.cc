#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

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

// Define our command-line flags.
DEFINE_string(input_file, "",
              "Input systematics file");
DEFINE_bool(galaxy_radec, false, "Galaxy coordinates are in RA-DEC");
DEFINE_string(map_file, "",
              "Name of the ASCII file containing the input map.");
DEFINE_string(sysmap_file, "",
              "Name of the ASCII file to write systematics map to.");
DEFINE_int32(resolution, 256,
	     "Resolution to use for creating the systematics map.");
DEFINE_bool(single_index, false, "Use older single-index file input format.");
DEFINE_bool(no_weight, false, "Input map file is missing weight column.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --input_file=<input systematics file>";
  usage += " --star_file=<input star systematics file>";
  usage += " --map_file=<StompMap ASCII>";
  usage += " --sysmap_file=<SysMap ASCII>";
  usage += "\n Optional Arguments:";
  usage += "--strmap_file=<SysMap ASCII>\n:\tSave Resultant Star Map\n";
  usage += "--resolution=int_32\n:\tSet Resolution of output map\n";
  usage += "--star_res=int_32\n:\tSet Resolution of Star Map\n";   
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_input_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

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

  // First, we read our STOMP map into a map object.  There are a couple
  // permutations based on the various map formats that are out there: with
  // or without a weight column or in the single index or double index format.
  std::cout << "Reading in initial map...\n";
  Stomp::Map* stomp_map;
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
    stomp_map->Area() << " sq. deg.\n";

  // To prime our map, we find the coverage map of the input geometry at the
  // systematics map resolution.
  Stomp::PixelVector coverage_map;
  std::cout<<"Making "<<FLAGS_resolution<<" Coverage Map\n";
  stomp_map->Coverage(coverage_map, FLAGS_resolution);
  

  //stomp_map->Clear();
  //delete stomp_map;

  // Now we create a map object from pixel index to SysPixel...
  SysDict sysmap;

  //... and initialize it with our coverage pixels.
  for (Stomp::PixelIterator iter=coverage_map.begin();
       iter!=coverage_map.end();++iter) {
    Stomp::SysPixel pix(iter->Resolution(), iter->Pixnum(), iter->Weight());
    sysmap[iter->Pixnum()] = pix;
  }
  
  // Now we start reading in our systematics files.
  int32_t n_obj = 0, n_keep = 0;
  std::cout << "Parsing " << FLAGS_input_file.c_str() << "...\n";
  std::ifstream systematics_file(FLAGS_input_file.c_str());
  SysDictIterator sys_iter;
  double theta, phi;
  double seeing, extinction, sky, odds;
  //double inv_sky_2;

  Stomp::AngularCoordinate::Sphere galaxy_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_galaxy_radec) galaxy_sphere = Stomp::AngularCoordinate::Equatorial;
  
  while (!systematics_file.eof()) {
    systematics_file >> theta >> phi >> seeing >> extinction >> sky >> odds;
    
    if (!systematics_file.eof()) {
      n_obj++;
      Stomp::AngularCoordinate ang(theta, phi, galaxy_sphere);
      Stomp::Pixel pix(ang, FLAGS_resolution, 1.0);
      //inv_sky_2 =  stomp_map->FindAverageWeight(pix);
      sys_iter = sysmap.find(pix.Pixnum());
      if (sys_iter != sysmap.end()) {
	sys_iter->second.AddSysParameters(seeing, extinction, 
					  sky, odds);
	n_keep++;
      }
    }
  }
  systematics_file.close();
  std::cout << "Kept " << n_keep << "/" << n_obj << " systematics objects...\n";
 
  std::cout << "FINDING Average Weights...";
  stomp_map->FindAverageWeight(coverage_map);
  
  std::cout << "Writing systematics map to " << FLAGS_sysmap_file << "\n";

  std::ofstream output_file(FLAGS_sysmap_file.c_str());
  // Now we iterate over our coverage pixels and write out the results
  for (Stomp::PixelIterator iter=coverage_map.begin();
       iter!=coverage_map.end();++iter) {
    SysDictIterator sys_iter = sysmap.find(iter->Pixnum());
    if (sys_iter != sysmap.end()) {
      output_file << sys_iter->second.HPixnum() << " " <<
	sys_iter->second.Superpixnum() << " " <<
	sys_iter->second.Resolution() << " " <<
	sys_iter->second.UnmaskedFraction() << " " <<
	sys_iter->second.MeanSeeing() << " " <<
	sys_iter->second.MeanExtinction() << " " <<
	sqrt(1.0/iter->Weight()) << " " <<
        //sys_iter->second.MeanSky() << " " <<
	sys_iter->second.MeanOdds() << " " <<
	sys_iter->second.NObjects() << "\n";
    }
  }

  output_file.close();

  return 0;
}
