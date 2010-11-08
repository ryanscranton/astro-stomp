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
      sum_galaxies_ = 0.0;
      n_obj_ = 0;
    };
    SysPixel(uint16_t resolution, uint32_t hpixnum, uint16_t superpixnum,
	     double unmasked_fraction, double mean_seeing,
	     double mean_extinction, double mean_sky, uint16_t n_objects) {
      SetResolution(resolution);

      uint32_t tmp_x, tmp_y;
      Pixel::HPix2XY(resolution, hpixnum, superpixnum, tmp_x, tmp_y);
      SetPixnumFromXY(tmp_x, tmp_y);

      SetWeight(unmasked_fraction);

      sum_seeing_ = mean_seeing*n_objects;
      sum_extinction_ = mean_extinction*n_objects;
      sum_sky_ = mean_sky*n_objects;
      sum_galaxies_ = 0.0;

      n_obj_ = n_objects;
    };
    ~SysPixel() {
      SetResolution(0);

      SetPixnumFromXY(0, 0);
      SetWeight(0.0);

      sum_seeing_ = 0.0;
      sum_extinction_ = 0.0;
      sum_sky_ = 0.0;
      sum_galaxies_ = 0.0;

      n_obj_ = 0;
    }
    inline void AddGalaxies(double galaxies) {
      sum_galaxies_ += galaxies;
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
    inline double TotalGalaxies() {
      return (n_obj_ > 0 ? sum_galaxies_ : 0.0);
    };
    inline double UnmaskedFraction() {
      return Weight();
    };
    inline double UnmaskedArea() {
      return Weight()*Area();
    };
    inline uint16_t NObjects() {
      return n_obj_;
    };

  private:
    double sum_seeing_, sum_extinction_, sum_sky_, sum_galaxies_;
    uint16_t n_obj_;
  };
} // end namespace Stomp

typedef std::map<const unsigned long, Stomp::SysPixel> SysDict;
typedef SysDict::iterator SysDictIterator;

// Define our command-line flags.
DEFINE_string(output_file, "",
              "Name of ASCII output file");
DEFINE_string(gal_file, "",
              "Name of the ASCII file containing the galaxies.");
DEFINE_string(sysmap_file, "",
              "Name of the ASCII file containing the systematics map.");
DEFINE_int32(resolution, 256,
	     "resolution of galaxy map");
DEFINE_double(seeing_min, 0.5, "Minimum allowable seeing value");
DEFINE_double(seeing_max, 2.0, "Maximum allowable seeing value");
DEFINE_double(extinction_min, 0.0, "Minimum allowable extinction value");
DEFINE_double(extinction_max, 1.0, "Maximum allowable extinction value");
DEFINE_double(sky_min, 0.0, "Minimum allowable sky brightness value");
DEFINE_double(sky_max, 24.0, "Maximum allowable sky brightness value");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --output_file=<Output file suffix>";
  usage += " --gal_file=<Catalog ASCII>";
  usage += " --sysmap_file=<SysMap ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_output_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

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

  // First, we create a map object from pixel index to SysPixel...
  SysDict sysmap;

  // Now we read in our sysmap.  We only want to keep those pixels that
  // meet our selection criteria, so we filter on those as we read through
  // the file.
  std::cout << "Reading in systematics map from " <<
    FLAGS_sysmap_file << "...\n";
  std::ifstream sysmap_file(FLAGS_sysmap_file.c_str());
  double unmasked, seeing, extinction, sky;
  uint32_t hpixnum, superpixnum;
  uint16_t resolution, n_objects;

  uint32_t n_pixel = 0, n_keep = 0;
  while (!sysmap_file.eof()) {
    sysmap_file >> hpixnum >> superpixnum >> resolution >> unmasked >>
      seeing >> extinction >> sky >> n_objects;

    if (!sysmap_file.eof()) {
      n_pixel++;
      if ((seeing >= FLAGS_seeing_min) &&
	  (seeing <= FLAGS_seeing_max) &&
	  (extinction >= FLAGS_extinction_min) &&
	  (extinction <= FLAGS_extinction_max) &&
	  (sky >= FLAGS_sky_min) &&
	  (sky <= FLAGS_sky_max) &&
	  (n_objects > 0)) {
	Stomp::SysPixel pix(resolution, hpixnum, superpixnum, unmasked,
			    seeing, extinction, sky, n_objects);
	sysmap[pix.Pixnum()] = pix;
	n_keep++;
      }
    }
  }
  sysmap_file.close();
  std::cout << "Kept " << n_keep << "/" << n_pixel <<
    " systematics pixels...\n";

  // Now we scan through our galaxy map and add the galaxy total to each
  // corresponding systematics pixel.
  std::ifstream gal_file(FLAGS_gal_file.c_str());
  long double ra, dec, mag;
  std::cout << "Writing output catalog to " << FLAGS_output_file << "\n";
  std::ofstream out_file(FLAGS_output_file.c_str());

  uint32_t matched_galaxies = 0;
  uint32_t n_galaxies = 0;
  double id;
  std::cout << "Reading galaxies map from " << FLAGS_gal_file << "...\n";
  while (!gal_file.eof()) {
    gal_file >> id >> ra >> dec >> mag;

    n_galaxies++;
    Stomp::AngularCoordinate ang(ra,dec,
				 Stomp::AngularCoordinate::Equatorial);
    Stomp::Pixel tmp_pix(ang, FLAGS_resolution, 1.0);
    SysDictIterator sys_iter = sysmap.find(tmp_pix.Pixnum());     
    if (sys_iter != sysmap.end()) {
      out_file << id << " " << ra << " "  << dec << " " << mag << " " <<
	sys_iter->second.MeanSeeing() << " " <<
	sys_iter->second.MeanExtinction() << " " <<
	sys_iter->second.MeanSky() << "\n";
      matched_galaxies++;
    } else {
      out_file << id << " " << ra << " " << dec << " " << mag << " " <<
	-99 << " " << -99 << " " << -99 << "\n";
    }
  }
  gal_file.close();
  out_file.close();
  std::cout << "Done.  Matched " << matched_galaxies << "/" << n_galaxies <<
    " galaxies...\n";

  return 0;
}
