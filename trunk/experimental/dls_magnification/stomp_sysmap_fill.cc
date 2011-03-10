#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

/*
Takes as input two systematics maps and attempts to write the systematics of a lower resolution map onto a higher one (and maybe the other direction). boolen variables set which systematics field to overwrite between the maps. Output is a systematics map a the resolution of the input sysmap.
 */

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
    inline void SetSysParameters(double mean_seeing, double mean_extinction,
				 double mean_sky, double mean_odds,
				 int n_objects) {
      sum_seeing_ =  mean_seeing*n_objects;
      sum_extinction_ = mean_extinction*n_objects;
      sum_sky_ = mean_sky*n_objects;
      sum_odds_ = mean_odds*n_objects;
      n_obj_ = n_objects;
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
DEFINE_string(output_file, "",
              "Output systematics file");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --input_file=<input systematics file>";
  usage += " --output_file=<output systematics file>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);
  
  if (FLAGS_input_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }
  if (FLAGS_output_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  SysDict sysmap;
  Stomp::PixelVector coverage_map;
  Stomp::PixelVector empty_pixels;

  std::cout << "Parsing " << FLAGS_input_file.c_str() << "...\n";
  std::ifstream input_file(FLAGS_input_file.c_str());
  uint32_t hpixnum, superpixnum, pixnum, x, y;
  double unmasked_fraction, seeing, extinction, sky_noise, odds;
  int resolution, n_objects;
  int n_empty_pixels = 0;
  long int n_objects_total = 0;
  long int n_pixels = 0;
  double unmasked_total = 0;

  if (input_file) {
    while (!input_file.eof()) {
      input_file >> hpixnum >> superpixnum >> resolution >> unmasked_fraction >>
	seeing >> extinction >> sky_noise >> odds >> n_objects;
      n_objects_total += n_objects;
      unmasked_total += unmasked_fraction;
      n_pixels++;
      Stomp::Pixel::HPix2XY(static_cast<uint32_t>(resolution), hpixnum, 
			    superpixnum, x, y);
      Stomp::Pixel tmp_pix(x, y, resolution, unmasked_fraction);
      Stomp::SysPixel pix(tmp_pix.Resolution(), tmp_pix.Pixnum(),
			  tmp_pix.Weight(), seeing, extinction, sky_noise,
			  odds, n_objects);
      sysmap[tmp_pix.Pixnum()] = pix;
      coverage_map.push_back(tmp_pix);
      if (n_objects == 0 && unmasked_fraction>0.50) {
	empty_pixels.push_back(tmp_pix);
	n_empty_pixels++;
      }
    }
  }
  std::cout << "Need to Fill " << n_empty_pixels << " empty pixels...\n";
  
  for (Stomp::PixelIterator iter=empty_pixels.begin();
       iter!=empty_pixels.end();++iter) {
    resolution = iter->Resolution();
    Stomp::PixelVector child_pixels;
    int super_resolution = resolution/2;
    int n_objects_hold = 0;

    bool run_loop = true;
    while (run_loop) {
      if (super_resolution < 128) {
	std::cout << "WARNING::Resolution required to corse (<128).\n"
		  << "\t Breaking Loop and leaving pixel empty.\n";
	break;
      }
      iter->SetToSuperPix(static_cast<uint32_t>(super_resolution));
      iter->SubPix(static_cast<uint32_t>(resolution), child_pixels);
      double sum_seeing = 0;
      double sum_extinction = 0;
      double sum_sky =0;
      double sum_odds = 0;
      int n_objects = 0;
      double unmasked_total = 0;

      for (Stomp::PixelIterator child_iter = child_pixels.begin();
	   child_iter != child_pixels.end(); ++child_iter) {
	SysDictIterator child_sys_iter = sysmap.find(child_iter->Pixnum());
	if (child_sys_iter != sysmap.end()) {
	  if (child_sys_iter->second.NObjects() > 0) {
	    sum_seeing += child_sys_iter->second.UnmaskedFraction()*
	      child_sys_iter->second.MeanSeeing();
	    sum_extinction += child_sys_iter->second.UnmaskedFraction()*
	      child_sys_iter->second.MeanExtinction();
	    sum_sky += child_sys_iter->second.UnmaskedFraction()*
	      child_sys_iter->second.MeanSky();
	    sum_odds += child_sys_iter->second.UnmaskedFraction()*
	      child_sys_iter->second.MeanOdds();
	    unmasked_total += child_sys_iter->second.UnmaskedFraction();
	    n_objects += child_sys_iter->second.NObjects();
	  }
	}
      }
      if (n_objects >= 2) {
	for (Stomp::PixelIterator child_iter = child_pixels.begin();
	     child_iter != child_pixels.end(); ++child_iter) {
	  SysDictIterator child_sys_iter = sysmap.find(child_iter->Pixnum());
	  if (child_sys_iter != sysmap.end()) {
	    child_sys_iter->second.SetSysParameters(sum_seeing/unmasked_total,
						    sum_extinction/
						    unmasked_total,
						    sum_sky/unmasked_total,
						    sum_odds/unmasked_total,
						    ceil(n_objects/
							 unmasked_total));
	    run_loop = false;
	  }
	}
      }
      n_objects_hold = n_objects;
      super_resolution /= 2;
    }
  }

  std::cout << "Writing systematics map to " << FLAGS_output_file << "\n";

  std::ofstream output_file(FLAGS_output_file.c_str());
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
        sys_iter->second.MeanSky() << " " <<
	sys_iter->second.MeanOdds() << " " <<
	sys_iter->second.NObjects() << "\n";
    }
  }

  return 0;
}

