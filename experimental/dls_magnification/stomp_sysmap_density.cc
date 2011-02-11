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
      sum_galaxies_ = 0.0;
      n_obj_ = 0;
    };
    SysPixel(uint16_t resolution, uint32_t hpixnum, uint16_t superpixnum,
	     double unmasked_fraction, double mean_seeing,
	     double mean_extinction, double mean_sky, double mean_odds, uint16_t n_objects) {
      SetResolution(resolution);

      uint32_t tmp_x, tmp_y;
      Pixel::HPix2XY(resolution, hpixnum, superpixnum, tmp_x, tmp_y);
      SetPixnumFromXY(tmp_x, tmp_y);

      SetWeight(unmasked_fraction);

      sum_seeing_ = mean_seeing*n_objects;
      sum_extinction_ = mean_extinction*n_objects;
      sum_sky_ = mean_sky*n_objects;
      sum_odds_ = mean_odds*n_objects;
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
      sum_odds_ = 0.0;
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
    inline double MeanOdds() {
      return (n_obj_ > 0 ? sum_odds_/n_obj_ : -1.0);
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
    double sum_seeing_, sum_extinction_, sum_sky_, sum_odds_, sum_galaxies_;
    uint16_t n_obj_;
  };
} // end namespace Stomp

class AreaBin {
public:
  AreaBin() {
    total_area_ = 0.0;
    total_galaxies_ = 0.0;
    bin_min_ = -1.0;
    bin_max_ = -1.0;
  };
  AreaBin(double bin_min, double bin_max) {
    total_area_ = 0.0;
    total_galaxies_ = 0.0;
    bin_min_ = bin_min;
    bin_max_ = bin_max;
  };
  ~AreaBin() {
    total_area_ = 0.0;
    total_galaxies_ = 0.0;
    bin_min_ = -1.0;
    bin_max_ = -1.0;
  };
  void AddToBin(double galaxies, double area) {
    total_galaxies_ += galaxies;
    total_area_ += area;
  };
  double BinMinimum() {
    return bin_min_;
  };
  double BinMaximum() {
    return bin_max_;
  };
  double BinCenter() {
    return 0.5*(bin_min_ + bin_max_);
  };
  bool WithinBin(double bin_value) {
    return (Stomp::DoubleLE(bin_value, bin_max_) &&
	    Stomp::DoubleGE(bin_value, bin_min_) ? true : false);
  };
  double BinGalaxies() {
    return total_galaxies_;
  };
  double BinArea() {
    return total_area_;
  };
  double BinDensity() {
    return (total_area_ > 0.001 ? total_galaxies_/total_area_ : 0.0);
  };
  double BinDensityError() {
    return (total_area_ > 0.001 ? sqrt(total_galaxies_)/total_area_ : 0.0);
  };

private:
  double bin_min_, bin_max_, total_galaxies_, total_area_;
};

typedef std::vector<AreaBin> BinVector;
typedef BinVector::iterator BinIterator;

class SysDensity {
public:
  SysDensity() {
    sys_min_ = -1.0;
    sys_max_ = -1.0;
    n_bins_ = 0;
  };
  SysDensity(double sys_min, double sys_max, int n_bins) {
    sys_min_ = sys_min;
    sys_max_ = sys_max;
    n_bins_ = n_bins;

    bins_.reserve(n_bins);
    double dbin = (sys_max_ - sys_min_)/n_bins_;
    for (int i=0;i<n_bins_;i++)
      bins_.push_back(AreaBin(sys_min_ + dbin*i, sys_min_ + dbin*(i+1)));
  };
  ~SysDensity() {
    sys_min_ = -1.0;
    sys_max_ = -1.0;
    n_bins_ = 0;

    bins_.clear();
  };
  bool AddToBin(double sys, double galaxies, double area) {
    bool found_bin = false;

    // Do a simple binary search on the bins to find the appropriate one.
    BinIterator top = bins_.end();
    BinIterator bottom = bins_.begin();
    BinIterator iter;
    --top;

    if ((sys < bottom->BinMinimum()) ||
	(sys > top->BinMaximum())) {
      iter = bins_.end();
    } else {
      ++top;
      --bottom;
      while (top-bottom > 1) {
	iter = bottom + (top - bottom)/2;
	if (sys < iter->BinMinimum()) {
	  top = iter;
	} else {
	  bottom = iter;
	}
      }
      iter = bottom;
    }

    // Provided that we've found the bin where the input belongs, add the
    // galaxies and area to that bin.
    if (iter!=bins_.end() && iter->WithinBin(sys)) {
      iter->AddToBin(galaxies, area);
      found_bin = true;
    }

    return found_bin;
  };
  BinIterator Begin() {
    return bins_.begin();
  };
  BinIterator End() {
    return bins_.end();
  };

private:
  double sys_min_, sys_max_;
  int n_bins_;
  BinVector bins_;
};

typedef std::map<const unsigned long, Stomp::SysPixel> SysDict;
typedef SysDict::iterator SysDictIterator;

// Define our command-line flags.
DEFINE_string(output_tag, "",
              "Output files will be e.g. SeeingDensity_OUTPUT_TAG.");
DEFINE_string(gal_file, "",
              "Name of the ASCII file containing the galaxy map.");
DEFINE_string(sysmap_file, "",
              "Name of the ASCII file containing the systematics map.");
DEFINE_bool(galaxy_radec, false, "Galaxy coordinates are in RA-DEC");
DEFINE_int32(resolution, -1,
	     "resolution of galaxy map");
DEFINE_int32(n_bins, 20,
	     "Number of bins to use when calculating densities.");
DEFINE_double(mag_min,18,
	      "Minimum allowable Magnitude");
DEFINE_double(mag_max,28,
	      "Maximum allowable Magnitude");
DEFINE_double(seeing_min, 0.5, "Minimum allowable seeing value");
DEFINE_double(seeing_max, 2.0, "Maximum allowable seeing value");
DEFINE_double(extinction_min, 0.0, "Minimum allowable extinction value");
DEFINE_double(extinction_max, 1.0, "Maximum allowable extinction value");
DEFINE_double(sky_min, 0.0, "Minimum allowable sky brightness value");
DEFINE_double(sky_max, 24.0, "Maximum allowable sky brightness value");
DEFINE_double(odds_min, 0.0, "Minimum allowable BPZ odds value");
DEFINE_double(odds_max, 1.0, "Maximum allowable BPZ odds value");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --output_tag=<Output file suffix>";
  usage += " --gal_file=<GalMap ASCII>";
  usage += " --sysmap_file=<SysMap ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_output_tag.empty()) {
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
  double unmasked, seeing, extinction, sky, odds;
  uint32_t hpixnum, superpixnum;
  uint16_t resolution, n_objects;

  uint32_t n_pixel = 0, n_keep = 0;
  while (!sysmap_file.eof()) {
    sysmap_file >> hpixnum >> superpixnum >> resolution >> unmasked >>
      seeing >> extinction >> sky >> odds >> n_objects;

    if (!sysmap_file.eof()) {
      n_pixel++;
      if ((seeing >= FLAGS_seeing_min) &&
	  (seeing <= FLAGS_seeing_max) &&
	  (extinction >= FLAGS_extinction_min) &&
	  (extinction <= FLAGS_extinction_max) &&
	  (sky >= FLAGS_sky_min) &&
	  (sky <= FLAGS_sky_max) &&
	  (odds >= FLAGS_odds_min) &&
	  (odds <= FLAGS_odds_max) &&
	  (n_objects > 0)) {
	Stomp::SysPixel pix(resolution, hpixnum, superpixnum, unmasked,
			    seeing, extinction, sky, odds, n_objects);
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
  double theta, phi, weight, mag;
  Stomp::AngularCoordinate::Sphere galaxy_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_galaxy_radec) galaxy_sphere = Stomp::AngularCoordinate::Equatorial;

  uint32_t matched_galaxies = 0;
  uint32_t n_galaxies = 0;
  uint16_t galmap_resolution = resolution;
  if (FLAGS_resolution != -1) {
    galmap_resolution = FLAGS_resolution;
  } 
  std::cout << "Reading galaxies map from " << FLAGS_gal_file << "...\n";
  while (!gal_file.eof()) {
    gal_file >> theta >> phi >> weight >> mag;

    /*std::cout << "Sucsessfully Read in galaxy hpixnum: " << hpixnum << " "
              << "superpixnum: " << superpixnum << " "
              << "resolution: " << resolution << " "
              << "UnmaskedFraction: " << unmasked << " "
              << "galaxies: " << galaxies << " "
              << "n_objects: " << n_objects << " "
	      << "...\n";
    */

    if (!gal_file.eof()) {
      if (mag>FLAGS_mag_min && mag<FLAGS_mag_max) {
	n_galaxies++;
	Stomp::AngularCoordinate ang(theta,phi,
				     galaxy_sphere);
	Stomp::Pixel tmp_pix(ang, FLAGS_resolution, 1.0);
	SysDictIterator sys_iter = sysmap.find(tmp_pix.Pixnum());     
	if (sys_iter != sysmap.end()) {
	  sys_iter->second.AddGalaxies(1);
	  matched_galaxies++;
	}
      }
    }
    //if (hpixnum==1728 && superpixnum==6158)
    //  break;
  }
  gal_file.close();
  std::cout << "Done.  Matched " << matched_galaxies << "/" << n_galaxies <<
    " galaxies...\n";

  std::cout << "Assigning galaxies to density bins...\n";
  SysDensity seeing_density(FLAGS_seeing_min, FLAGS_seeing_max,
			    FLAGS_n_bins);
  SysDensity extinction_density(FLAGS_extinction_min, FLAGS_extinction_max,
				FLAGS_n_bins);
  SysDensity sky_density(FLAGS_sky_min, FLAGS_sky_max,
			 FLAGS_n_bins);
  SysDensity odds_density(FLAGS_odds_min, FLAGS_odds_max,
			  FLAGS_n_bins);

  uint32_t kept_seeing = 0;
  uint32_t kept_extinction = 0;
  uint32_t kept_sky = 0;
  uint32_t kept_odds = 0;
  for (SysDictIterator iter=sysmap.begin();iter!=sysmap.end();++iter) {
    if (seeing_density.AddToBin(iter->second.MeanSeeing(),
				iter->second.TotalGalaxies(),
				iter->second.UnmaskedArea())) kept_seeing++;
    if (extinction_density.AddToBin(iter->second.MeanExtinction(),
				    iter->second.TotalGalaxies(),
				    iter->second.UnmaskedArea()))
      kept_extinction++;
    if (sky_density.AddToBin(iter->second.MeanSky(),
			     iter->second.TotalGalaxies(),
			     iter->second.UnmaskedArea())) kept_sky++;
    if (odds_density.AddToBin(iter->second.MeanOdds(),
			      iter->second.TotalGalaxies(),
			      iter->second.UnmaskedArea())) kept_odds++;
  }
  std::cout << "Done; " << kept_seeing << "/" << kept_extinction << "/" <<
    kept_sky << "/" << kept_odds << "/" << n_keep << ".\n";

  // Write out the results.
  std::string seeing_file_name = FLAGS_output_tag + "SeeingDensity.ascii";
  std::cout << "Writing seeing density to " << seeing_file_name << "\n";
  std::ofstream seeing_file(seeing_file_name.c_str());
  double total_area = 0.0;
  for (BinIterator iter=seeing_density.Begin();
       iter!=seeing_density.End();++iter) {
    total_area += iter->BinArea();
    seeing_file << iter->BinCenter() << " " <<
      iter->BinDensity() << " " << iter->BinDensityError() << " " <<
      iter->BinArea() << " " << total_area << "\n";
  }
  seeing_file.close();

  std::string extinction_file_name = FLAGS_output_tag + "ExtinctionDensity.ascii";
  std::cout << "Writing extinction density to " << extinction_file_name << "\n";
  std::ofstream extinction_file(extinction_file_name.c_str());
  total_area = 0.0;
  for (BinIterator iter=extinction_density.Begin();
       iter!=extinction_density.End();++iter) {
    total_area += iter->BinArea();
    extinction_file << iter->BinCenter() << " " <<
      iter->BinDensity() << " " << iter->BinDensityError() << " " <<
      iter->BinArea() << " " << total_area << "\n";
  }
  extinction_file.close();

  std::string sky_file_name = FLAGS_output_tag + "SkyDensity.ascii";
  std::cout << "Writing sky density to " << sky_file_name << "\n";
  std::ofstream sky_file(sky_file_name.c_str());
  total_area = 0.0;
  for (BinIterator iter=sky_density.Begin();
       iter!=sky_density.End();++iter) {
    total_area += iter->BinArea();
    sky_file << iter->BinCenter() << " " <<
      iter->BinDensity() << " " << iter->BinDensityError() << " " <<
      iter->BinArea() << " " << total_area << "\n";
  }
  sky_file.close();

  std::string odds_file_name = FLAGS_output_tag + "OddsDensity.ascii";
  std::cout << "Writing odds density to " << odds_file_name << "\n";
  std::ofstream odds_file(odds_file_name.c_str());
  total_area = 0.0;
  for (BinIterator iter=odds_density.Begin();
       iter!=odds_density.End();++iter) {
    total_area += iter->BinArea();
    odds_file << iter->BinCenter() << " " <<
      iter->BinDensity() << " " << iter->BinDensityError() << " " <<
      iter->BinArea() << " " << total_area << "\n";
  }
  odds_file.close();

  std::cout << "Done.\n";

  return 0;
}
