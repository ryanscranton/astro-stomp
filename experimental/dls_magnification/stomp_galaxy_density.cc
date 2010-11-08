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

      galaxies_ = 0.0;
      n_obj_ = 0;
    };
    SysPixel(uint32_t resolution, uint32_t pixnum, double unmasked_fraction) {
      SetResolution(resolution);

      uint32_t tmp_y = pixnum/(Stomp::Nx0*Resolution());
      uint32_t tmp_x = pixnum - Stomp::Nx0*Resolution()*tmp_y;

      SetPixnumFromXY(tmp_x, tmp_y);
      SetWeight(unmasked_fraction);

      galaxies_ = 0.0;
      n_obj_ = 0;
    };
    SysPixel(uint32_t resolution, uint32_t pixnum, double unmasked_fraction,
	     int galaxies, int n_objects) {
      SetResolution(resolution);

      uint32_t tmp_y = pixnum/(Stomp::Nx0*Resolution());
      uint32_t tmp_x = pixnum - Stomp::Nx0*Resolution()*tmp_y;

      SetPixnumFromXY(tmp_x, tmp_y);
      SetWeight(unmasked_fraction);

      galaxies_ = galaxies;
      n_obj_ = n_objects;
    };
    ~SysPixel() {
      SetResolution(0);

      SetPixnumFromXY(0, 0);
      SetWeight(0.0);
 
      galaxies_ = 0.0;
      n_obj_ = 0;
    }
    inline void AddSysParameters(int galaxies) {
      galaxies_ += galaxies;
      n_obj_++;
    };
    inline double Galaxies() {
      return galaxies_;
    };
    inline double UnmaskedFraction() {
      return Weight();
    };
    inline int NObjects() {
      return n_obj_;
    };

  private:
    int galaxies_, n_obj_;
  };
} // end namespace Stomp

typedef std::map<const uint32_t, Stomp::SysPixel> SysDict;
typedef SysDict::iterator SysDictIterator;

// Define our command-line flags.
DEFINE_string(input_files, "",
              "CSV list of Galaxy Files");
DEFINE_string(map_file, "",
              "Name of the ASCII file containing the input map.");
DEFINE_string(galaxy_file, "",
              "Name of the ASCII file to write galaxies map to.");
DEFINE_int32(resolution, 256,
	     "Resolution to use for creating the galaxies map.");
DEFINE_bool(single_index, false, "Use older single-index file input format.");
DEFINE_bool(no_weight, false, "Input map file is missing weight column.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --input_files=<CSV list of input galaxies files>";
  usage += " --map_file=<StompMap ASCII>";
  usage += " --galaxy_file=<SysMap ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_input_files.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_map_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_galaxy_file.empty()) {
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
  // galaxies map resolution.
  Stomp::PixelVector coverage_map;
  stomp_map->Coverage(coverage_map, FLAGS_resolution);

  stomp_map->Clear();
  delete stomp_map;

  // Now we create a map object from pixel index to SysPixel...
  SysDict galmap;

  //... and initialize it with our coverage pixels.
  for (Stomp::PixelIterator iter=coverage_map.begin();
       iter!=coverage_map.end();++iter) {
    Stomp::SysPixel pix(iter->Resolution(), iter->Pixnum(), iter->Weight());
    galmap[iter->Pixnum()] = pix;
  }

  // Now we start reading in our galaxies files.
  std::vector<std::string> input_files;

  Stomp::Tokenize(FLAGS_input_files, input_files, ",");

  std::cout << "Parsing " << input_files.size() << " files...\n";

  uint32_t n_obj = 0, n_keep = 0;
  for (std::vector<std::string>::iterator file_iter=input_files.begin();
       file_iter!=input_files.end();++file_iter) {

    std::cout << "\tParsing " << file_iter->c_str() << "...\n";
    std::ifstream galaxies_file(file_iter->c_str());
    double lambda, eta;
    uint32_t pixnum;
    SysDictIterator gal_iter;

    while (!galaxies_file.eof()) {
      galaxies_file >> lambda >> eta;

      if (!galaxies_file.eof()) {
	n_obj++;
	Stomp::AngularCoordinate ang(lambda, eta,
				     Stomp::AngularCoordinate::Equatorial);
	Stomp::Pixel pix(ang, FLAGS_resolution, 1.0);
	gal_iter = galmap.find(pix.Pixnum());
	if (gal_iter != galmap.end()) {
	  gal_iter->second.AddSysParameters(1);
	  n_keep++;
	}
      }
    }
    galaxies_file.close();
  }
  std::cout << "Kept " << n_keep << "/" << n_obj << " galaxies...\n";

  std::cout << "Writing density map to " << FLAGS_galaxy_file << "\n";

  std::ofstream output_file(FLAGS_galaxy_file.c_str());
  // Now we iterate over our coverage pixels and write out the results
  for (Stomp::PixelIterator iter=coverage_map.begin();
       iter!=coverage_map.end();++iter) {
    SysDictIterator gal_iter = galmap.find(iter->Pixnum());
    if (gal_iter != galmap.end()) {
      output_file << 
	gal_iter->second.HPixnum() << " " <<
	gal_iter->second.Superpixnum() << " " <<
	gal_iter->second.Resolution() << " " <<
	gal_iter->second.UnmaskedFraction() << " " <<
	gal_iter->second.Galaxies() << " " <<
	gal_iter->second.NObjects() << "\n";
    }
  }

  output_file.close();

  return 0;
}
