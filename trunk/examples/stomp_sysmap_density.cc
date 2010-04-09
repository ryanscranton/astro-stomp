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
    SysPixel(uint32_t resolution, uint32_t hpixnum, uint32_t superpixnum,
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

typedef std::map<const uint32_t, Stomp::SysPixel> SysDict;
typedef SysDict::iterator SysDictIterator;

// Define our command-line flags.
DEFINE_string(output_tag, "",
              "Output files will be e.g. SeeingDensity_OUTPUT_TAG.");
DEFINE_string(galmap_file, "",
              "Name of the ASCII file containing the galaxy map.");
DEFINE_string(sysmap_file, "",
              "Name of the ASCII file containing the systematics map.");
DEFINE_int32(n_bins, 20,
	     "Number of bins to use when calculating densities.");
DEFINE_double(seeing_min, 0.5, "Minimum allowable seeing value");
DEFINE_double(seeing_max, 2.0, "Maximum allowable seeing value");
DEFINE_double(extinction_min, 0.0, "Minimum allowable extinction value");
DEFINE_double(extinction_max, 1.0, "Maximum allowable extinction value");
DEFINE_double(sky_min, 0.0, "Minimum allowable sky brightness value");
DEFINE_double(sky_max, 24.0, "Maximum allowable sky brightness value");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --output_tag=<Output file suffix>";
  usage += " --galmap_file=<GalMap ASCII>";
  usage += " --sysmap_file=<SysMap ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_output_tag.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_galmap_file.empty()) {
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
  uint32_t hpixnum, superpixnum, resolution, n_objects;

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
  std::ifstream galmap_file(FLAGS_galmap_file.c_str());
  double galaxies;

  uint32_t matched_pixels = 0;
  std::cout << "Reading galaxies map from " << FLAGS_galmap_file << "...\n";
  while (!galmap_file.eof()) {
    galmap_file >> hpixnum >> superpixnum >> resolution >> unmasked >>
      galaxies >> n_objects;

    if (!galmap_file.eof()) {
      Stomp::Pixel tmp_pix(resolution, hpixnum, superpixnum);
      SysDictIterator sys_iter = sysmap.find(tmp_pix.Pixnum());
      if (sys_iter != sysmap.end()) {
	sys_iter->second.AddGalaxies(galaxies);
	matched_pixels++;
      }
    }
  }
  galmap_file.close();
  std::cout << "Done.  Matched " << matched_pixels << "/" << n_keep <<
    " pixels...\n";

  std::cout << "Assigning galaxies to gensity bins...\n";
  SysDensity seeing_density(FLAGS_seeing_min, FLAGS_seeing_max,
			    FLAGS_n_bins);
  SysDensity extinction_density(FLAGS_extinction_min, FLAGS_extinction_max,
				FLAGS_n_bins);
  SysDensity sky_density(FLAGS_sky_min, FLAGS_sky_max,
			 FLAGS_n_bins);

  uint32_t kept_seeing = 0;
  uint32_t kept_extinction = 0;
  uint32_t kept_sky = 0;
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
  }
  std::cout << "Done; " << kept_seeing << "/" << kept_extinction << "/" <<
    kept_sky << "/" << n_keep << ".\n";

  // Write out the results.
  std::string seeing_file_name = "SeeingDensity_" + FLAGS_output_tag;
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

  std::string extinction_file_name = "ExtinctionDensity_" + FLAGS_output_tag;
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

  std::string sky_file_name = "SkyDensity_" + FLAGS_output_tag;
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

  std::cout << "Done.\n";

  return 0;
}
