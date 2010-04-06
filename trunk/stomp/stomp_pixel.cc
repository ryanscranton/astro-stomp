// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the central class for the library: the hierarchical
// pixelization that makes all of the rest of the spatial classes work.  This
// class defines all of the core Pixel operations: converting angular position
// on the sphere into pixel index, increasing and decreasing pixel resolution,
// finding neighboring pixels and so on.

#include "stomp_core.h"
#include "stomp_pixel.h"
#include "stomp_angular_coordinate.h"
#include "stomp_angular_bin.h"

namespace Stomp {

Pixel::Pixel() {
  level_ = 0;
  y_ = 0;
  x_ = 0;
  weight_ = 0.0;
}

Pixel::Pixel(const uint32_t input_resolution,
	     const uint32_t input_pixnum,
	     const double input_weight) {
  if ((input_resolution < HPixResolution) ||
      (input_resolution%2 != 0) ||
      (input_resolution > MaxPixelResolution))
    std::cout << "Invalid resolution value.\n ";

  if (input_pixnum > MaxPixnum)
    std::cout << "Invalid pixel index value.\n ";

  SetResolution(input_resolution);
  y_ = input_pixnum/(Nx0*Resolution());
  x_ = input_pixnum - Nx0*Resolution()*y_;
  weight_ = input_weight;
}

Pixel::Pixel(const uint32_t input_x,
	     const uint32_t input_y,
	     const uint32_t input_resolution,
	     const double input_weight) {
  if ((input_resolution < HPixResolution) ||
      (input_resolution%2 != 0) ||
      (input_resolution > MaxPixelResolution)) {
    std::cout << "Invalid resolution value: " << input_resolution << ".\n";
    exit(2);
  }

  if (input_x > Nx0*input_resolution)
    std::cout << "Invalid x index value.\n";

  if (input_y > Ny0*input_resolution)
    std::cout << "Invalid y index value.\n";

  SetResolution(input_resolution);
  x_ = input_x;
  y_ = input_y;
  weight_ = input_weight;
}

Pixel::Pixel(AngularCoordinate& ang, const uint32_t input_resolution,
	     const double input_weight) {

  if ((input_resolution < HPixResolution) ||
      (input_resolution % 2 != 0) ||
      (input_resolution > MaxPixelResolution)) {
    std::cout << "Invalid resolution value: " << input_resolution << ".\n";
    exit(2);
  }

  SetResolution(input_resolution);

  double eta = (ang.Eta() - EtaOffSet)*DegToRad;

  if (eta <= 0.0) eta += 2.0*Pi;

  eta /= 2.0*Pi;

  x_ = static_cast<uint32_t>(Nx0*Resolution()*eta);

  double lambda = (90.0 - ang.Lambda())*DegToRad;

  if (lambda >= Pi) {
    y_ = Ny0*Resolution() - 1;
  } else {
    y_ = static_cast<uint32_t>(Ny0*Resolution()*
			       ((1.0 - cos(lambda))/2.0));
  }

  weight_ = input_weight;
}

Pixel::~Pixel() {
  x_ = y_ = 0;
  level_ = 0;
  weight_ = 0.0;
}

bool Pixel::operator<(Pixel& pix) {
  if (this->Resolution() == pix.Resolution()) {
    if (this->PixelY() == pix.PixelY()) {
      return (this->PixelX() < pix.PixelX() ? true : false);
    } else {
      return (this->PixelY() < pix.PixelY() ? true : false);
    }
  } else {
    return (this->Resolution() < pix.Resolution() ? true : false);
  }
}

bool Pixel::operator==(Pixel& pix) {
  return (this->Resolution() == pix.Resolution() &&
	  this->PixelY() == pix.PixelY() &&
	  this->PixelX() == pix.PixelX() ? true : false);
}

bool Pixel::operator!=(Pixel& pix) {
  return (this->Resolution() != pix.Resolution() ||
	  this->PixelY() != pix.PixelY() ||
	  this->PixelX() != pix.PixelX() ? true : false);
}

Pixel Pixel::operator++() {
  if (x_ == Nx0*Resolution() - 1) {
    x_ = 0;
    y_++;
  } else {
    if (y_ < Ny0*Resolution()) x_++;
  }
  return *this;
}

void Pixel::SetPixnumFromAng(AngularCoordinate& ang) {

  double eta = (ang.Eta() - EtaOffSet)*DegToRad;

  if (eta <= 0.0) eta += 2.0*Pi;

  eta /= 2.0*Pi;
  x_ = static_cast<uint32_t>(Nx0*Resolution()*eta);

  double lambda = (90.0 - ang.Lambda())*DegToRad;

  if (lambda >= Pi) {
    y_ = Ny0*Resolution() - 1;
  } else {
    y_ = static_cast<uint32_t>(Ny0*Resolution()*
                                    ((1.0 - cos(lambda))/2.0));
  }
}

void Pixel::SetResolution(uint32_t resolution) {
  level_ = MostSignificantBit(resolution);
  x_ = 0;
  y_ = 0;
}

void Pixel::SetLevel(uint8_t level) {
  level_ = level;
  x_ = 0;
  y_ = 0;
}

void Pixel::SetPixnumFromXY(uint32_t x, uint32_t y) {
  x_ = x;
  y_ = y;
}

uint8_t Pixel::Level() {
  return level_;
}

uint32_t Pixel::Resolution() {
  return static_cast<uint32_t>(1 << level_);
}

double Pixel::Weight() {
  return weight_;
}

void Pixel::SetWeight(double weight) {
  weight_ = weight;
}

void Pixel::ReverseWeight() {
  weight_ *= -1.0;
}

void Pixel::InvertWeight() {
  weight_ = 1.0/weight_;
}

uint32_t Pixel::PixelX() {
  return x_;
}

uint32_t Pixel::PixelY() {
  return y_;
}

bool Pixel::SetToSuperPix(uint32_t lo_resolution) {
  bool success = false;
  if (Resolution() < lo_resolution) {
    std::cout << "Illegal resolution value: " << lo_resolution <<
      " < " << Resolution();
  } else {
    x_ /= Resolution()/lo_resolution;
    y_ /= Resolution()/lo_resolution;
    level_ = MostSignificantBit(lo_resolution);

    success = true;
  }

  return success;
}

bool Pixel::SetToLevel(uint8_t lo_level) {
  bool success = false;
  if (Level() < lo_level) {
    std::cout << "Illegal level value: " << lo_level << " < " << Level();
  } else {
    x_ /= Resolution()/Level2Resolution(lo_level);
    y_ /= Resolution()/Level2Resolution(lo_level);
    level_ = lo_level;

    success = true;
  }

  return success;
}

void Pixel::SubPix(uint32_t hi_resolution, PixelVector& pix) {
  if (!pix.empty()) pix.clear();

  uint32_t x_min, x_max, y_min, y_max;
  SubPix(hi_resolution, x_min, x_max, y_min, y_max);

  for (uint32_t y=y_min;y<=y_max;y++) {
    for (uint32_t x=x_min;x<=x_max;x++) {
      Pixel tmp_pix(x, y, hi_resolution, weight_);
      pix.push_back(tmp_pix);
    }
  }
}

void Pixel::SubPix(uint32_t hi_resolution, uint32_t& x_min,
		   uint32_t& x_max, uint32_t& y_min,
		   uint32_t& y_max) {
  if (Resolution() == hi_resolution) {
    y_min = y_max = y_;
    x_min = x_max = x_;
  } else {
    uint32_t tmp_res;

    x_min = x_max = x_;
    y_min = y_max = y_;
    for (tmp_res=Resolution();tmp_res<hi_resolution;tmp_res*=2) {
      x_min *= 2;
      y_min *= 2;

      x_max = 2*x_max + 1;
      y_max = 2*y_max + 1;
    }
  }
}

void Pixel::CohortPix(Pixel& pix_a, Pixel& pix_b, Pixel& pix_c) {
  uint32_t super_x, super_y, x1, x2, x3, x4, y1, y2, y3, y4;

  pix_a.SetResolution(Resolution());
  pix_b.SetResolution(Resolution());
  pix_c.SetResolution(Resolution());

  super_x = x_/2;
  super_y = y_/2;

  x1 = 2*super_x;
  y1 = 2*super_y;

  x2 = 2*super_x + 1;
  y2 = 2*super_y;

  x3 = 2*super_x;
  y3 = 2*super_y + 1;

  x4 = 2*super_x + 1;
  y4 = 2*super_y + 1;

  if ((x_ == x1) && (y_ == y1)) {
    pix_a.SetPixnumFromXY(x2,y2);
    pix_b.SetPixnumFromXY(x3,y3);
    pix_c.SetPixnumFromXY(x4,y4);
  }
  if ((x_ == x2) && (y_ == y2)) {
    pix_a.SetPixnumFromXY(x1,y1);
    pix_b.SetPixnumFromXY(x3,y3);
    pix_c.SetPixnumFromXY(x4,y4);
  }
  if ((x_ == x3) && (y_ == y3)) {
    pix_b.SetPixnumFromXY(x1,y1);
    pix_a.SetPixnumFromXY(x2,y2);
    pix_c.SetPixnumFromXY(x4,y4);
  }
  if ((x_ == x4) && (y_ == y4)) {
    pix_c.SetPixnumFromXY(x1,y1);
    pix_a.SetPixnumFromXY(x2,y2);
    pix_b.SetPixnumFromXY(x3,y3);
  }
}

bool Pixel::FirstCohort() {
  return (2*(x_/2) == x_ && 2*(y_/2) == y_ ? true : false);
}

double Pixel::Area() {
  return HPixArea*HPixResolution*HPixResolution/
    (Resolution()*Resolution());
}

uint32_t Pixel::SuperPix(uint32_t lo_resolution) {
  return (Resolution() < lo_resolution ?
	  Nx0*Ny0*lo_resolution*lo_resolution :
	  Nx0*lo_resolution*
	  static_cast<uint32_t>(y_*lo_resolution/Resolution()) +
	  static_cast<uint32_t>(x_*lo_resolution/Resolution()));
}

uint32_t Pixel::Superpixnum() {
  return Nx0*HPixResolution*
    static_cast<uint32_t>(y_*HPixResolution/Resolution()) +
    static_cast<uint32_t>(x_*HPixResolution/Resolution());
}

uint32_t Pixel::HPixnum() {
  return
    Resolution()/HPixResolution*
    (y_ - Resolution()/HPixResolution*
     static_cast<uint32_t>(y_*HPixResolution/Resolution())) +
    (x_ - Resolution()/HPixResolution*
     static_cast<uint32_t>(x_*HPixResolution/Resolution()));
}

uint32_t Pixel::Pixnum() {
  return Nx0*Resolution()*y_ + x_;
}

bool Pixel::Contains(uint32_t pixel_resolution, uint32_t pixel_x,
		     uint32_t pixel_y) {
  return ((pixel_resolution >= Resolution()) &&
	  (pixel_x*Resolution()/pixel_resolution == x_) &&
	  (pixel_y*Resolution()/pixel_resolution == y_) ? true : false);
}

bool Pixel::Contains(Pixel& pix) {
  return ((pix.Resolution() >= Resolution()) &&
	  (pix.PixelX()*Resolution()/pix.Resolution() == x_) &&
	  (pix.PixelY()*Resolution()/pix.Resolution() == y_) ? true : false);
}

bool Pixel::Contains(AngularCoordinate& ang) {
  double eta = (ang.Eta() - EtaOffSet)*DegToRad;
  if (eta <= 0.0) eta += 2.0*Pi;
  eta /= 2.0*Pi;

  if (x_ == static_cast<uint32_t>(Nx0*Resolution()*eta)) {
    double lambda = (90.0 - ang.Lambda())*DegToRad;
    if (lambda >= Pi) {
      return (y_ == Ny0*Resolution() - 1 ? true : false);
    } else {
      return (y_ == static_cast<uint32_t>(Ny0*Resolution()*
					  ((1.0 - cos(lambda))/2.0)) ?
	      true : false);
    }
  } else {
    return false;
  }
}

bool Pixel::WithinBounds(double lon_min, double lon_max,
			 double lat_min, double lat_max,
			 AngularCoordinate::Sphere sphere) {
  double pix_lon_min = 0.0;
  double pix_lon_max = 0.0;
  double pix_lat_min = 0.0;
  double pix_lat_max = 0.0;

  switch (sphere) {
  case AngularCoordinate::Survey:
    pix_lon_min = EtaMin();
    pix_lon_max = EtaMax();
    pix_lat_min = LambdaMin();
    pix_lat_max = LambdaMax();
    break;
  case AngularCoordinate::Equatorial:
    pix_lon_min = RAMin();
    pix_lon_max = RAMax();
    pix_lat_min = DECMin();
    pix_lat_max = DECMax();
    break;
  case AngularCoordinate::Galactic:
    pix_lon_min = GalLonMin();
    pix_lon_max = GalLonMax();
    pix_lat_min = GalLatMin();
    pix_lat_max = GalLatMax();
    break;
  }

  return (DoubleLE(pix_lon_max, lon_max) &&
	  DoubleGE(pix_lon_min, lon_min) &&
	  DoubleLE(pix_lat_max, lat_max) &&
	  DoubleGE(pix_lat_min, lat_min) ? true : false);
}

bool Pixel::IntersectsBounds(double lon_min, double lon_max,
			     double lat_min, double lat_max,
			     AngularCoordinate::Sphere sphere) {
  double pix_lon_min = 0.0;
  double pix_lon_max = 0.0;
  double pix_lat_min = 0.0;
  double pix_lat_max = 0.0;

  switch (sphere) {
  case AngularCoordinate::Survey:
    pix_lon_min = EtaMin();
    pix_lon_max = EtaMax();
    pix_lat_min = LambdaMin();
    pix_lat_max = LambdaMax();
    break;
  case AngularCoordinate::Equatorial:
    pix_lon_min = RAMin();
    pix_lon_max = RAMax();
    pix_lat_min = DECMin();
    pix_lat_max = DECMax();
    break;
  case AngularCoordinate::Galactic:
    pix_lon_min = GalLonMin();
    pix_lon_max = GalLonMax();
    pix_lat_min = GalLatMin();
    pix_lat_max = GalLatMax();
    break;
  }

  // Either the pixel contains the bounds or the bounds contain the pixel.
  bool contained = ((DoubleLE(pix_lon_max, lon_max) &&
		     DoubleGE(pix_lon_min, lon_min) &&
		     DoubleLE(pix_lat_max, lat_max) &&
		     DoubleGE(pix_lat_min, lat_min)) ||
		    (DoubleLE(lon_max, pix_lon_max) &&
		     DoubleGE(lon_min, pix_lon_min) &&
		     DoubleLE(lat_max, pix_lat_max) &&
		     DoubleGE(lat_min, pix_lat_min)));

  // Check to see if any of the corners of the bound are within the pixel.
  bool corner = ((DoubleLE(pix_lon_min, lon_max) &&
		  DoubleGE(pix_lon_min, lon_min) &&
		  DoubleLE(pix_lat_max, lat_max) &&
		  DoubleGE(pix_lat_max, lat_min)) ||
		 (DoubleLE(pix_lon_max, lon_max) &&
		  DoubleGE(pix_lon_max, lon_min) &&
		  DoubleLE(pix_lat_max, lat_max) &&
		  DoubleGE(pix_lat_max, lat_min)) ||
		 (DoubleLE(pix_lon_min, lon_max) &&
		  DoubleGE(pix_lon_min, lon_min) &&
		  DoubleLE(pix_lat_min, lat_max) &&
		  DoubleGE(pix_lat_min, lat_min)) ||
		 (DoubleLE(pix_lon_max, lon_max) &&
		  DoubleGE(pix_lon_max, lon_min) &&
		  DoubleLE(pix_lat_min, lat_max) &&
		  DoubleGE(pix_lat_min, lat_min)));

  // Check the cases where the bounds cut through the pixel horizontally, either
  // passing all the way through or taking a chunk out of the left or right
  // side of the pixel.
  bool lat_middle = (((DoubleGE(pix_lon_min, lon_min) &&
		       DoubleLE(pix_lon_max, lon_max)) ||
		      (DoubleLE(pix_lon_min, lon_min) &&
		       DoubleGE(pix_lon_max, lon_min)) ||
		      (DoubleLE(pix_lon_min, lon_max) &&
		       DoubleGE(pix_lon_max, lon_max))) &&
		     DoubleGE(pix_lat_max, lat_max) &&
		     DoubleLE(pix_lat_min, lat_min));

  // Same as above, but for a vertical slice through the pixel.
  bool lon_middle = (((DoubleGE(pix_lat_min, lat_min) &&
		       DoubleLE(pix_lat_max, lat_max)) ||
		      (DoubleLE(pix_lat_min, lat_min) &&
		       DoubleGE(pix_lat_max, lat_min)) ||
		      (DoubleLE(pix_lat_min, lat_max) &&
		       DoubleGE(pix_lat_max, lat_max))) &&
		     DoubleLE(pix_lon_min, lon_min) &&
		     DoubleGE(pix_lon_max, lon_max));

  // If any of our cases hold, send back True.  False, otherwise.
  return (contained || corner || lon_middle || lat_middle ? true : false);
}

void Pixel::WithinRadius(double theta_max, PixelVector& pix,
			 bool check_full_pixel) {
  AngularBin theta(0.0, theta_max);
  WithinAnnulus(theta, pix, check_full_pixel);
}

void Pixel::WithinAnnulus(double theta_min, double theta_max,
			  PixelVector& pix, bool check_full_pixel) {
  AngularBin theta(theta_min, theta_max);
  WithinAnnulus(theta, pix, check_full_pixel);
}

void Pixel::WithinAnnulus(AngularBin& theta, PixelVector& pix,
			  bool check_full_pixel) {
  if (!pix.empty()) pix.clear();

  uint32_t y_min;
  uint32_t y_max;
  std::vector<uint32_t> x_min;
  std::vector<uint32_t> x_max;

  XYBounds(theta.ThetaMax(), x_min, x_max, y_min, y_max, true);

  uint32_t nx = Nx0*Resolution();
  uint32_t nx_pix;

  for (uint32_t y=y_min,n=0;y<=y_max;y++,n++) {
    if ((x_max[n] < x_min[n]) && (x_min[n] > nx/2)) {
      nx_pix = nx - x_min[n] + x_max[n] + 1;
    } else {
      nx_pix = x_max[n] - x_min[n] + 1;
    }
    if (nx_pix > nx) nx_pix = nx;
    for (uint32_t m=0,x=x_min[n];m<nx_pix;m++,x++) {
      if (x == nx) x = 0;
      Pixel tmp_pix(x, y, Resolution(), 1.0);
      bool within_bounds =
	theta.WithinCosBounds(UnitSphereX()*tmp_pix.UnitSphereX() +
			      UnitSphereY()*tmp_pix.UnitSphereY() +
			      UnitSphereZ()*tmp_pix.UnitSphereZ());
      if (check_full_pixel && within_bounds) {
	if (theta.WithinCosBounds(UnitSphereX()*tmp_pix.UnitSphereX_UL() +
				  UnitSphereY()*tmp_pix.UnitSphereY_UL() +
				  UnitSphereZ()*tmp_pix.UnitSphereZ_UL()) &&
	    theta.WithinCosBounds(UnitSphereX()*tmp_pix.UnitSphereX_UR() +
				  UnitSphereY()*tmp_pix.UnitSphereY_UR() +
				  UnitSphereZ()*tmp_pix.UnitSphereZ_UR()) &&
	    theta.WithinCosBounds(UnitSphereX()*tmp_pix.UnitSphereX_LL() +
				  UnitSphereY()*tmp_pix.UnitSphereY_LL() +
				  UnitSphereZ()*tmp_pix.UnitSphereZ_LL()) &&
	    theta.WithinCosBounds(UnitSphereX()*tmp_pix.UnitSphereX_LR() +
				  UnitSphereY()*tmp_pix.UnitSphereY_LR() +
				  UnitSphereZ()*tmp_pix.UnitSphereZ_LR())) {
	  within_bounds = true;
	} else {
	  within_bounds = false;
	}
      }
      if (within_bounds) pix.push_back(tmp_pix);
    }
  }
}

void Pixel::BoundingRadius(double theta_max, PixelVector& pix) {
  AngularCoordinate ang;
  Ang(ang);

  BoundingRadius(ang, theta_max, pix);
}

void Pixel::BoundingRadius(AngularCoordinate& ang, double theta_max,
			   PixelVector& pix) {
  SetPixnumFromAng(ang);
  if (!pix.empty()) pix.clear();

  uint32_t y_min;
  uint32_t y_max;
  std::vector<uint32_t> x_min;
  std::vector<uint32_t> x_max;

  XYBounds(ang, theta_max, x_min, x_max, y_min, y_max, true);

  uint32_t nx = Nx0*Resolution();
  uint32_t nx_pix;

  for (uint32_t y=y_min,n=0;y<=y_max;y++,n++) {
    if ((x_max[n] < x_min[n]) && (x_min[n] > nx/2)) {
      nx_pix = nx - x_min[n] + x_max[n] + 1;
    } else {
      nx_pix = x_max[n] - x_min[n] + 1;
    }
    if (nx_pix > nx) nx_pix = nx;
    for (uint32_t m=0,x=x_min[n];m<nx_pix;m++,x++) {
      if (x == nx) x = 0;
      Pixel tmp_pix(x, y, Resolution(), 1.0);
      pix.push_back(tmp_pix);
    }
  }
}

void Pixel::XYBounds(double theta, uint32_t& x_min,
		     uint32_t& x_max, uint32_t& y_min,
		     uint32_t& y_max, bool add_buffer) {
  double lammin = Lambda() - theta;
  if (lammin < -90.0) lammin = -90.0;
  double lammax = Lambda() + theta;
  if (lammax > 90.0) lammax = 90.0;

  double sphere_correction  = 1.0;
  if (fabs(lammin) > fabs(lammax)) {
    sphere_correction =
      1.0 + 0.000192312*lammin*lammin -
      1.82764e-08*lammin*lammin*lammin*lammin +
      1.28162e-11*lammin*lammin*lammin*lammin*lammin*lammin;
  } else {
    sphere_correction =
      1.0 + 0.000192312*lammax*lammax -
      1.82764e-08*lammax*lammax*lammax*lammax +
      1.28162e-11*lammax*lammax*lammax*lammax*lammax*lammax;
  }

  uint32_t nx = Nx0*Resolution();
  uint32_t ny = Ny0*Resolution();

  double etamin = Eta() - theta*sphere_correction;
  etamin -= EtaOffSet;
  etamin *= DegToRad;

  if (etamin <= 0.0) etamin = etamin + 2.0*Pi;

  etamin /= 2.0*Pi;
  x_min = static_cast<uint32_t>(nx*etamin);

  lammax = (90.0 - lammax)*DegToRad;

  if (lammax >= Pi) {
    y_min = ny - 1;
  } else {
    y_min = static_cast<uint32_t>(ny*((1.0 - cos(lammax))/2.0));
  }

  double etamax = Eta() + theta*sphere_correction;
  etamax -= EtaOffSet;
  etamax *= DegToRad;

  if (etamax <= 0.0) etamax = etamax + 2.0*Pi;

  etamax /= 2.0*Pi;
  x_max = static_cast<uint32_t>(nx*etamax);

  lammin = (90.0 - lammin)*DegToRad;

  if (lammin >= Pi) {
    y_max = ny - 1;
  } else {
    y_max = static_cast<uint32_t>(ny*((1.0 - cos(lammin))/2.0));
  }

  if (add_buffer) {
    if (x_min == 0) {
      x_min = nx - 1;
    } else {
      x_min--;
    }

    if (x_max == nx - 1) {
      x_max = 0;
    } else {
      x_max++;
    }

    if (y_max < ny - 1) y_max++;
    if (y_min > 0) y_min--;
  }
}

void Pixel::XYBounds(double theta, std::vector<uint32_t>& x_min,
		     std::vector<uint32_t>& x_max,
		     uint32_t& y_min, uint32_t& y_max,
		     bool add_buffer) {
  if (!x_min.empty()) x_min.clear();
  if (!x_max.empty()) x_max.clear();

  double lammin = Lambda() - theta;
  if (lammin < -90.0) lammin = -90.0;
  double lammax = Lambda() + theta;
  if (lammax > 90.0) lammax = 90.0;

  uint32_t ny = Ny0*Resolution();

  lammax = (90.0 - lammax)*DegToRad;

  if (lammax >= Pi) {
    y_min = ny - 1;
  } else {
    y_min = static_cast<uint32_t>(ny*((1.0 - cos(lammax))/2.0));
  }

  lammin = (90.0 - lammin)*DegToRad;

  if (lammin >= Pi) {
    y_max = ny - 1;
  } else {
    y_max = static_cast<uint32_t>(ny*((1.0 - cos(lammin))/2.0));
  }

  if (add_buffer) {
    if (y_max < ny - 1) y_max++;
    if (y_min > 0) y_min--;
  }

  if (!x_min.empty()) x_min.clear();
  if (!x_max.empty()) x_max.clear();

  x_min.reserve(y_max-y_min+1);
  x_max.reserve(y_max-y_min+1);

  uint32_t nx = Nx0*Resolution();

  for (uint32_t y=y_min,n=0;y<=y_max;y++,n++) {
    double lam = 90.0 -
      RadToDeg*acos(1.0 - 2.0*(y+0.5)/(Ny0*Resolution()));

    double sphere_correction = 1.0 +
      lam*lam*(0.000192312 - lam*lam*(1.82764e-08 - 1.28162e-11*lam*lam));

    double etamin = Eta() - theta*sphere_correction;
    etamin -= EtaOffSet;
    etamin *= DegToRad;

    if (etamin <= 0.0) etamin = etamin + 2.0*Pi;

    etamin /= 2.0*Pi;
    x_min.push_back(static_cast<uint32_t>(nx*etamin));

    double etamax = Eta() + theta*sphere_correction;
    etamax -= EtaOffSet;
    etamax *= DegToRad;

    if (etamax <= 0.0) etamax = etamax + 2.0*Pi;

    etamax /= 2.0*Pi;
    x_max.push_back(static_cast<uint32_t>(nx*etamax));

    if (add_buffer) {
      if (x_min[n] == 0) {
	x_min[n] = nx - 1;
      } else {
	x_min[n] -= 1;
      }

      if (x_max[n] == nx - 1) {
	x_max[n] = 0;
      } else {
	x_max[n] += 1;
      }
    }
  }
}

void Pixel::XYBounds(AngularCoordinate& ang, double theta,
		     uint32_t& x_min, uint32_t& x_max,
		     uint32_t& y_min, uint32_t& y_max,
		     bool add_buffer) {
  double lammin = ang.Lambda() - theta;
  if (lammin < -90.0) lammin = -90.0;
  double lammax = ang.Lambda() + theta;
  if (lammax > 90.0) lammax = 90.0;

  double sphere_correction  = 1.0;
  if (fabs(lammin) > fabs(lammax)) {
    sphere_correction =
      1.0 + 0.000192312*lammin*lammin -
      1.82764e-08*lammin*lammin*lammin*lammin +
      1.28162e-11*lammin*lammin*lammin*lammin*lammin*lammin;
  } else {
    sphere_correction =
      1.0 + 0.000192312*lammax*lammax -
      1.82764e-08*lammax*lammax*lammax*lammax +
      1.28162e-11*lammax*lammax*lammax*lammax*lammax*lammax;
  }

  uint32_t nx = Nx0*Resolution();
  uint32_t ny = Ny0*Resolution();

  double etamin = ang.Eta() - theta*sphere_correction;
  etamin -= EtaOffSet;
  etamin *= DegToRad;

  if (etamin <= 0.0) etamin = etamin + 2.0*Pi;

  etamin /= 2.0*Pi;
  x_min = static_cast<uint32_t>(nx*etamin);

  lammax = (90.0 - lammax)*DegToRad;

  if (lammax >= Pi) {
    y_min = ny - 1;
  } else {
    y_min = static_cast<uint32_t>(ny*((1.0 - cos(lammax))/2.0));
  }

  double etamax = ang.Eta() + theta*sphere_correction;
  etamax -= EtaOffSet;
  etamax *= DegToRad;

  if (etamax <= 0.0) etamax = etamax + 2.0*Pi;

  etamax /= 2.0*Pi;
  x_max = static_cast<uint32_t>(nx*etamax);

  lammin = (90.0 - lammin)*DegToRad;

  if (lammin >= Pi) {
    y_max = ny - 1;
  } else {
    y_max = static_cast<uint32_t>(ny*((1.0 - cos(lammin))/2.0));
  }

  if (add_buffer) {
    if (x_min == 0) {
      x_min = nx - 1;
    } else {
      x_min--;
    }

    if (x_max == nx - 1) {
      x_max = 0;
    } else {
      x_max++;
    }

    if (y_max < ny - 1) y_max++;
    if (y_min > 0) y_min--;
  }
}

void Pixel::XYBounds(AngularCoordinate& ang, double theta,
		     std::vector<uint32_t>& x_min,
		     std::vector<uint32_t>& x_max,
		     uint32_t& y_min, uint32_t& y_max,
		     bool add_buffer) {
  if (!x_min.empty()) x_min.clear();
  if (!x_max.empty()) x_max.clear();

  double lammin = ang.Lambda() - theta;
  if (lammin < -90.0) lammin = -90.0;
  double lammax = ang.Lambda() + theta;
  if (lammax > 90.0) lammax = 90.0;

  uint32_t ny = Ny0*Resolution();

  lammax = (90.0 - lammax)*DegToRad;
  if (lammax >= Pi) {
    y_min = ny - 1;
  } else {
    y_min = static_cast<uint32_t>(ny*((1.0 - cos(lammax))/2.0));
  }

  lammin = (90.0 - lammin)*DegToRad;
  if (lammin >= Pi) {
    y_max = ny - 1;
  } else {
    y_max = static_cast<uint32_t>(ny*((1.0 - cos(lammin))/2.0));
  }

  if (add_buffer) {
    if (y_max < ny - 1) y_max++;
    if (y_min > 0) y_min--;
  }

  if (y_max > ny -1) {
    std::cout << "Illegal theta value on y index: Lambda,Eta = " <<
      ang.Lambda() << "," << ang.Eta() << ", theta = " << theta <<
      "\ny = " << y_max << "/" << ny - 1 << "\n";
    exit(2);
  }

  if (y_min > ny -1) {
    std::cout << "Illegal theta value on y index: Lambda,Eta = " <<
      ang.Lambda() << "," << ang.Eta() << ", theta = " << theta <<
      "\ny = " << y_min << "/" << ny - 1 << "\n";
    exit(2);
  }

  if (!x_min.empty()) x_min.clear();
  if (!x_max.empty()) x_max.clear();

  x_min.reserve(y_max-y_min+1);
  x_max.reserve(y_max-y_min+1);

  uint32_t nx = Nx0*Resolution();
  for (uint32_t y=y_min,n=0;y<=y_max;y++,n++) {
    double lam = 90.0 -
      RadToDeg*acos(1.0 - 2.0*(y+0.5)/(Ny0*Resolution()));

    double sphere_correction = 1.0 +
      lam*lam*(0.000192312 - lam*lam*(1.82764e-08 - 1.28162e-11*lam*lam));

    double etamin = ang.Eta() - theta*sphere_correction;
    etamin -= EtaOffSet;
    etamin *= DegToRad;
    if (etamin <= 0.0) etamin = etamin + 2.0*Pi;
    etamin /= 2.0*Pi;
    x_min.push_back(static_cast<uint32_t>(nx*etamin));

    double etamax = ang.Eta() + theta*sphere_correction;
    etamax -= EtaOffSet;
    etamax *= DegToRad;
    if (etamax <= 0.0) etamax = etamax + 2.0*Pi;
    etamax /= 2.0*Pi;
    x_max.push_back(static_cast<uint32_t>(nx*etamax));

    if (add_buffer) {
      if (x_min[n] == 0) {
	x_min[n] = nx - 1;
      } else {
	x_min[n] -= 1;
      }

      if (x_max[n] == nx - 1) {
	x_max[n] = 0;
      } else {
	x_max[n] += 1;
      }
    }

    if (x_max[n] > nx -1) {
      std::cout << "Illegal theta value on x index: Lambda,Eta = " <<
	ang.Lambda() << "," << ang.Eta() << ", theta = " << theta <<
	"\nx = " << x_max[n] << "/" << nx - 1 << "\n";
      exit(2);
    }
    if (x_min[n] > nx -1) {
      std::cout << "Illegal theta value on x index: Lambda,Eta = " <<
	ang.Lambda() << "," << ang.Eta() << ", theta = " << theta <<
	"\nx = " << x_min[n] << "/" << nx - 1 << "\n";
      exit(2);
    }
  }
}

uint8_t Pixel::EtaStep(double theta) {
  double lam = Lambda();

  double deta = 2.5/(Resolution()/4);
  double eta_step = theta;
  eta_step *= 1.0 +
    lam*lam*(0.000192312 - lam*lam*(1.82764e-08 - 1.28162e-11*lam*lam));
  uint8_t etastep = 1;
  while (eta_step > etastep*deta) etastep++;

  return etastep;
}

double Pixel::NearEdgeDistance(AngularCoordinate& ang) {
  // If the test position is within the lambda or eta ranges of the pixel,
  // then we need to do a proper calculation.  Otherwise, we can just return
  // the nearest corner distance.
  double min_edge_distance = -1.0;
  bool inside_bounds = false;
  if (DoubleLE(ang.Lambda(), LambdaMax()) &&
      DoubleGE(ang.Lambda(), LambdaMin())) {
    inside_bounds = true;

    double lam = ang.Lambda();
    double eta_scaling = 1.0 +
      lam*lam*(0.000192312 - lam*lam*(1.82764e-08 - 1.28162e-11*lam*lam));

    min_edge_distance = fabs(ang.Eta() - EtaMin());

    if (min_edge_distance > fabs(ang.Eta() - EtaMax()))
      min_edge_distance = fabs(ang.Eta() - EtaMax());
    min_edge_distance /= eta_scaling;
  }

  if (DoubleLE(ang.Eta(), EtaMax()) &&
      DoubleGE(ang.Eta(), EtaMin())) {
    double lambda_distance = fabs(ang.Lambda() - LambdaMax());
    if (lambda_distance > fabs(ang.Lambda() - LambdaMin()))
      lambda_distance = fabs(ang.Lambda() - LambdaMin());
    if (inside_bounds) {
      if (lambda_distance < min_edge_distance)
	min_edge_distance = lambda_distance;
    } else {
      min_edge_distance = lambda_distance;
    }
    inside_bounds = true;
  }

  if (!inside_bounds) {
    // If we're outside of those bounds, then the nearest part of the pixel
    // should be the near corner.
    min_edge_distance = NearCornerDistance(ang);
  } else {
    // The return value for this function is (sin(theta))^2 rather than just
    // the angle theta.  If we can get by with the small angle approximation,
    // then we use that for speed purposes.
    if (min_edge_distance < 3.0) {
      min_edge_distance = min_edge_distance*DegToRad;
    } else {
      min_edge_distance = sin(min_edge_distance*DegToRad);
    }
    min_edge_distance *= min_edge_distance;
  }

  return min_edge_distance;
}

double Pixel::FarEdgeDistance(AngularCoordinate& ang) {
  // If the test position is within the lambda or eta ranges of the pixel,
  // then we need to do a proper calculation.  Otherwise, we can just return
  // the farthest corner distance.
  double max_edge_distance = -1.0;
  bool inside_bounds = false;
  if (DoubleLE(ang.Lambda(), LambdaMax()) &&
      DoubleGE(ang.Lambda(), LambdaMin())) {
    inside_bounds = true;

    double lam = ang.Lambda();
    double eta_scaling = 1.0 +
      lam*lam*(0.000192312 - lam*lam*(1.82764e-08 - 1.28162e-11*lam*lam));

    max_edge_distance = fabs(ang.Eta() - EtaMin());

    if (max_edge_distance < fabs(ang.Eta() - EtaMax()))
      max_edge_distance = fabs(ang.Eta() - EtaMax());
    max_edge_distance /= eta_scaling;
  }

  if (DoubleLE(ang.Eta(), EtaMax()) &&
      DoubleGE(ang.Eta(), EtaMin())) {
    double lambda_distance = fabs(ang.Lambda() - LambdaMax());
    if (lambda_distance > fabs(ang.Lambda() - LambdaMin()))
      lambda_distance = fabs(ang.Lambda() - LambdaMin());
    if (inside_bounds) {
      if (lambda_distance > max_edge_distance)
	max_edge_distance = lambda_distance;
    } else {
      max_edge_distance = lambda_distance;
    }
    inside_bounds = true;
  }

  if (!inside_bounds) {
    // If we're outside of those bounds, then the farthest part of the pixel
    // should be the far corner.
    max_edge_distance = FarCornerDistance(ang);
  } else {
    // The return value for this function is (sin(theta))^2 rather than just
    // the angle theta.  If we can get by with the small angle approximation,
    // then we use that for speed purposes.
    if (max_edge_distance < 3.0) {
      max_edge_distance = max_edge_distance*DegToRad;
    } else {
      max_edge_distance = sin(max_edge_distance*DegToRad);
    }
    max_edge_distance *= max_edge_distance;
  }
  return max_edge_distance;
}

double Pixel::NearCornerDistance(AngularCoordinate& ang) {
  double costheta_final = ang.UnitSphereX()*UnitSphereX_UL() +
    ang.UnitSphereY()*UnitSphereY_UL() +
    ang.UnitSphereZ()*UnitSphereZ_UL();

  double costheta = ang.UnitSphereX()*UnitSphereX_UR() +
    ang.UnitSphereY()*UnitSphereY_UR() +
    ang.UnitSphereZ()*UnitSphereZ_UR();
  if (costheta > costheta_final) costheta_final = costheta;

  costheta = ang.UnitSphereX()*UnitSphereX_LR() +
    ang.UnitSphereY()*UnitSphereY_LR() +
    ang.UnitSphereZ()*UnitSphereZ_LR();
  if (costheta > costheta_final) costheta_final = costheta;

  costheta = ang.UnitSphereX()*UnitSphereX_LL() +
    ang.UnitSphereY()*UnitSphereY_LL() +
    ang.UnitSphereZ()*UnitSphereZ_LL();
  if (costheta > costheta_final) costheta_final = costheta;

  // The return value for this function is (sin(theta))^2 rather than just
  // the angle theta for computing speed purposes.
  return 1.0 - costheta_final*costheta_final;
}

double Pixel::FarCornerDistance(AngularCoordinate& ang) {
  double costheta_final = ang.UnitSphereX()*UnitSphereX_UL() +
    ang.UnitSphereY()*UnitSphereY_UL() +
    ang.UnitSphereZ()*UnitSphereZ_UL();

  double costheta = ang.UnitSphereX()*UnitSphereX_UR() +
    ang.UnitSphereY()*UnitSphereY_UR() +
    ang.UnitSphereZ()*UnitSphereZ_UR();
  if (costheta < costheta_final) costheta_final = costheta;

  costheta = ang.UnitSphereX()*UnitSphereX_LR() +
    ang.UnitSphereY()*UnitSphereY_LR() +
    ang.UnitSphereZ()*UnitSphereZ_LR();
  if (costheta < costheta_final) costheta_final = costheta;

  costheta = ang.UnitSphereX()*UnitSphereX_LL() +
    ang.UnitSphereY()*UnitSphereY_LL() +
    ang.UnitSphereZ()*UnitSphereZ_LL();
  if (costheta < costheta_final) costheta_final = costheta;

  // The return value for this function is (sin(theta))^2 rather than just
  // the angle theta for computing speed purposes.
  return 1.0 - costheta_final*costheta_final;
}

bool Pixel::EdgeDistances(AngularCoordinate& ang, double& min_edge_distance,
			  double& max_edge_distance) {
  bool inside_bounds = false;
  double lambda_min = LambdaMin();
  double lambda_max = LambdaMax();
  double eta_min = EtaMin();
  double eta_max = EtaMax();
  double lam = ang.Lambda();
  double eta = ang.Eta();

  if (Stomp::DoubleLE(lam, lambda_max) && Stomp::DoubleGE(lam, lambda_min)) {
    inside_bounds = true;

    double eta_scaling = 1.0 +
      lam*lam*(0.000192312 - lam*lam*(1.82764e-08 - 1.28162e-11*lam*lam));

    min_edge_distance = fabs(eta - eta_min);
    max_edge_distance = fabs(eta - eta_min);

    if (min_edge_distance > fabs(eta - eta_max))
      min_edge_distance = fabs(eta - eta_max);
    if (max_edge_distance < fabs(eta - eta_max))
      max_edge_distance = fabs(eta - eta_max);

    min_edge_distance /= eta_scaling;
    max_edge_distance /= eta_scaling;
  }

  if (Stomp::DoubleLE(eta, eta_max) && Stomp::DoubleGE(eta, eta_min)) {
    if (inside_bounds) {
      if (min_edge_distance > fabs(lam - lambda_min))
	min_edge_distance = fabs(lam - lambda_min);
      if (min_edge_distance > fabs(lam - lambda_max))
	min_edge_distance = fabs(lam - lambda_max);

      if (max_edge_distance < fabs(lam - lambda_min))
	max_edge_distance = fabs(lam - lambda_min);
      if (max_edge_distance < fabs(lam - lambda_max))
	max_edge_distance = fabs(lam - lambda_max);
    } else {
      min_edge_distance = fabs(lam - lambda_min);
      max_edge_distance = fabs(lam - lambda_max);

      if (min_edge_distance > fabs(lam - lambda_max))
	min_edge_distance = fabs(lam - lambda_max);
      if (max_edge_distance < fabs(lam - lambda_max))
	max_edge_distance = fabs(lam - lambda_max);
    }
    inside_bounds = true;
  }

  if (inside_bounds) {
    // The return value for this function is (sin(theta))^2 rather than just
    // the angle theta.  If we can get by with the small angle approximation,
    // then we use that for speed purposes.
    if (min_edge_distance < 3.0) {
      min_edge_distance = min_edge_distance*Stomp::DegToRad;
    } else {
      min_edge_distance = sin(min_edge_distance*Stomp::DegToRad);
    }
    min_edge_distance *= min_edge_distance;

    if (max_edge_distance < 3.0) {
      max_edge_distance = max_edge_distance*Stomp::DegToRad;
    } else {
      max_edge_distance = sin(max_edge_distance*Stomp::DegToRad);
    }
    max_edge_distance *= max_edge_distance;
  }
  return inside_bounds;
}

bool Pixel::IsWithinRadius(AngularCoordinate& ang, double theta_max,
			   bool check_full_pixel) {
  AngularBin theta(0.0, theta_max);

  return IsWithinAnnulus(ang, theta, check_full_pixel);
}

bool Pixel::IsWithinRadius(Pixel& pix, double theta_max,
			   bool check_full_pixel) {
  AngularBin theta(0.0, theta_max);
  AngularCoordinate ang;
  pix.Ang(ang);

  return IsWithinAnnulus(ang, theta, check_full_pixel);
}

bool Pixel::IsWithinAnnulus(AngularCoordinate& ang, double theta_min,
			    double theta_max, bool check_full_pixel) {
  AngularBin theta(theta_min, theta_max);

  return IsWithinAnnulus(ang, theta, check_full_pixel);
}

bool Pixel::IsWithinAnnulus(Pixel& pix, double theta_min,
			    double theta_max, bool check_full_pixel) {
  AngularBin theta(theta_min, theta_max);
  AngularCoordinate ang;
  pix.Ang(ang);

  return IsWithinAnnulus(ang, theta, check_full_pixel);
}

bool Pixel::IsWithinAnnulus(AngularCoordinate& ang, AngularBin& theta,
			    bool check_full_pixel) {
  bool within_bounds = theta.WithinCosBounds(ang.UnitSphereX()*UnitSphereX() +
					     ang.UnitSphereY()*UnitSphereY() +
					     ang.UnitSphereZ()*UnitSphereZ());
  if (within_bounds && check_full_pixel) {
    if (theta.WithinCosBounds(ang.UnitSphereX()*UnitSphereX_UL() +
			      ang.UnitSphereY()*UnitSphereY_UL() +
			      ang.UnitSphereZ()*UnitSphereZ_UL()) &&
	theta.WithinCosBounds(ang.UnitSphereX()*UnitSphereX_UR() +
			      ang.UnitSphereY()*UnitSphereY_UR() +
			      ang.UnitSphereZ()*UnitSphereZ_UR()) &&
	theta.WithinCosBounds(ang.UnitSphereX()*UnitSphereX_LL() +
			      ang.UnitSphereY()*UnitSphereY_LL() +
			      ang.UnitSphereZ()*UnitSphereZ_LL()) &&
	theta.WithinCosBounds(ang.UnitSphereX()*UnitSphereX_LR() +
			      ang.UnitSphereY()*UnitSphereY_LR() +
			      ang.UnitSphereZ()*UnitSphereZ_LR())) {
      within_bounds = true;
    } else {
      within_bounds = false;
    }
  }
  return within_bounds;
}

bool Pixel::IsWithinAnnulus(Pixel& pix, AngularBin& theta,
			    bool check_full_pixel) {
  AngularCoordinate ang;
  pix.Ang(ang);

  return IsWithinAnnulus(ang, theta, check_full_pixel);
}

int8_t Pixel::IntersectsAnnulus(AngularCoordinate& ang,
				double theta_min, double theta_max) {
  AngularBin theta(theta_min, theta_max);

  return IntersectsAnnulus(ang, theta);
}

int8_t Pixel::IntersectsAnnulus(Pixel& pix, double theta_min,
				double theta_max) {
  AngularBin theta(theta_min, theta_max);
  AngularCoordinate ang;
  pix.Ang(ang);

  return IntersectsAnnulus(ang, theta);
}

int8_t Pixel::IntersectsAnnulus(AngularCoordinate& ang, AngularBin& theta) {
  // By default, we assume that there is no intersection between the disk
  // and the pixel.
  int8_t intersects_annulus = 0;

  double costheta =
    ang.UnitSphereX()*UnitSphereX() +
    ang.UnitSphereY()*UnitSphereY() +
    ang.UnitSphereZ()*UnitSphereZ();
  double costheta_max = costheta;
  double costheta_min = costheta;

  double costheta_ul =
    ang.UnitSphereX()*UnitSphereX_UL() +
    ang.UnitSphereY()*UnitSphereY_UL() +
    ang.UnitSphereZ()*UnitSphereZ_UL();
  if (costheta_ul > costheta_max) costheta_max = costheta_ul;
  if (costheta_ul < costheta_min) costheta_min = costheta_ul;

  double costheta_ll =
    ang.UnitSphereX()*UnitSphereX_LL() +
    ang.UnitSphereY()*UnitSphereY_LL() +
    ang.UnitSphereZ()*UnitSphereZ_LL();
  if (costheta_ll > costheta_max) costheta_max = costheta_ll;
  if (costheta_ll < costheta_min) costheta_min = costheta_ll;

  double costheta_ur =
    ang.UnitSphereX()*UnitSphereX_UR() +
    ang.UnitSphereY()*UnitSphereY_UR() +
    ang.UnitSphereZ()*UnitSphereZ_UR();
  if (costheta_ur > costheta_max) costheta_max = costheta_ur;
  if (costheta_ur < costheta_min) costheta_min = costheta_ur;

  double costheta_lr =
    ang.UnitSphereX()*UnitSphereX_LR() +
    ang.UnitSphereY()*UnitSphereY_LR() +
    ang.UnitSphereZ()*UnitSphereZ_LR();
  if (costheta_lr > costheta_max) costheta_max = costheta_lr;
  if (costheta_lr < costheta_min) costheta_min = costheta_lr;

  double near_corner_distance = 1.0 - costheta_max*costheta_max;
  double far_corner_distance = 1.0 - costheta_min*costheta_min;
  double near_edge_distance;
  double far_edge_distance;

  if (!EdgeDistances(ang, near_edge_distance, far_edge_distance)) {
    near_edge_distance = near_corner_distance;
    far_edge_distance = far_corner_distance;
  }

  bool contains_center = Contains(ang);

  // First we tackle the inner disk.
  int8_t intersects_inner_disk = 0;

  if (Stomp::DoubleGE(near_edge_distance, theta.Sin2ThetaMin())) {
    // If this is true, it means that the distance between the nearest edge
    // and the annulus center is greater than the inner annulus radius.  This
    // means that the inner disk is either completely inside or outside the
    // pixel.  Checking the center should tell us which is which.
    if (contains_center) {
      if (Stomp::DoubleEQ(theta.ThetaMin(), 0.0)) {
	// If the inner disk has radius = 0, then we've got a disk, not an
	// annulus.  Since the inner disk has no area, then it's counted as
	// outside the pixel and the intersect status will be determined only
	// by the outer disk.
	intersects_inner_disk = 0;
      } else {
	// disk is inside pixel.
	intersects_inner_disk = -1;
      }
    } else {
      // disk is outside pixel.
      intersects_inner_disk = 0;
    }
  } else {
    // If the distance to the nearest edge is less than the inner annulus
    // radius, then there is some intersection between the two; either the
    // disk intersects the pixel edge or it completely contains the pixel.
    // Checking the corners will tell is which is which.
    if (Stomp::DoubleGE(costheta_ul, theta.CosThetaMax()) &&
	Stomp::DoubleGE(costheta_ur, theta.CosThetaMax()) &&
	Stomp::DoubleGE(costheta_ll, theta.CosThetaMax()) &&
	Stomp::DoubleGE(costheta_lr, theta.CosThetaMax())) {
      // pixel is inside the disk.
      intersects_inner_disk = 1;
    } else {
      // pixel intersects the disk.
      intersects_inner_disk = -1;
    }
  }

  // Now, the outer disk.
  int8_t intersects_outer_disk = 0;

  if (Stomp::DoubleGE(near_edge_distance, theta.Sin2ThetaMax())) {
    // If this is true, it means that the distance between the nearest edge
    // and the annulus center is greater than the outer annulus radius.  This
    // means that the outer disk is either completely inside or outside the
    // pixel.  Checking the center should tell us which is which.
    if (contains_center) {
      // disk is inside pixel.
      intersects_outer_disk = -1;
    } else {
      // disk is outside pixel.
      intersects_outer_disk = 0;
    }
  } else {
    // If the distance to the nearest edge is less than the outer annulus
    // radius, then there is some intersection between the two; either the
    // disk intersects the pixel edge or it completely contains the pixel.
    // Checking the corners will tell is which is which.
    if (Stomp::DoubleGE(costheta_ul, theta.CosThetaMin()) &&
	Stomp::DoubleGE(costheta_ur, theta.CosThetaMin()) &&
	Stomp::DoubleGE(costheta_ll, theta.CosThetaMin()) &&
	Stomp::DoubleGE(costheta_lr, theta.CosThetaMin())) {
      // pixel is inside the disk.
      intersects_outer_disk = 1;
    } else {
      // pixel intersects the disk.
      intersects_outer_disk = -1;
    }
  }

  // Now we deal with cases.
  if ((intersects_inner_disk == 1) && (intersects_outer_disk == 1)) {
    // This means that the pixel is contained by the inner and outer disks.
    // Hence, the pixel is in the hole of the annulus.
    intersects_annulus = 0;
  }

  if ((intersects_inner_disk == -1) && (intersects_outer_disk == 1)) {
    // This means that the inner disk shares some area with the pixel and the
    // outer disk contains it completely.
    intersects_annulus = -1;
  }

  if ((intersects_inner_disk == 0) && (intersects_outer_disk == 1)) {
    // The inner disk is outside the pixel, but the outer disk contains it.
    // Hence, the pixel is fully inside the annulus.
    intersects_annulus = 1;
  }


  if ((intersects_inner_disk == 1) && (intersects_outer_disk == -1)) {
    // This should be impossible.  Raise an error.
    std::cout << "Impossible annulus intersection: " <<
      intersects_inner_disk << ", " << intersects_outer_disk << ".  Bailing.\n";
    exit(2);
  }

  if ((intersects_inner_disk == -1) && (intersects_outer_disk == -1)) {
    // There's partial overlap with both the inner and outer disks.
    intersects_annulus = -1;
  }

  if ((intersects_inner_disk == 0) && (intersects_outer_disk == -1)) {
    // The inner disk is outside, but the outer intersects, so there's some
    // intersection between the pixel and annulus.
    intersects_annulus = -1;
  }


  if ((intersects_inner_disk == 1) && (intersects_outer_disk == 0)) {
    // This should be impossible.  Raise an error.
    std::cout << "Impossible annulus intersection: " <<
      intersects_inner_disk << ", " << intersects_outer_disk << ".  Bailing.\n";
    exit(2);
  }

  if ((intersects_inner_disk == -1) && (intersects_outer_disk == 0)) {
    // This should be impossible.  Raise an error.
    std::cout << "Impossible annulus intersection: " <<
      intersects_inner_disk << ", " << intersects_outer_disk << ".  Bailing.\n";
    exit(2);
  }

  if ((intersects_inner_disk == 0) && (intersects_outer_disk == 0)) {
    // The inner disk is outside the pixel, and so is the outer disk.
    // Hence, the pixel is fully outside the annulus.
    intersects_annulus = 0;
  }

  return intersects_annulus;
}

int8_t Pixel::IntersectsAnnulus(Pixel& pix, AngularBin& theta) {
  AngularCoordinate ang;
  pix.Ang(ang);

  return IntersectsAnnulus(ang, theta);
}

uint32_t Pixel::Stripe(uint32_t input_resolution) {
  if ((input_resolution%2 != 0) || (input_resolution < 4)) {
    std::cout << "Illegal resolution in Stripe() call!\nExiting...\n";
    exit(1);
  }

  double stripe_width = 360.0/(Nx0*input_resolution);
  int32_t stripe = static_cast<int32_t>((Eta() + 32.5)/stripe_width) + 10;

  double etamin = stripe_width*(stripe - 10) - 32.5 -
      stripe_width/2.0 + 0.0000001;
  double etamax = stripe_width*(stripe - 10) - 32.5 +
      stripe_width/2.0 - 0.0000001;
  if (Eta() < etamin) stripe++;
  if (Eta() > etamax) stripe++;
  if (stripe < 0) stripe += Nx0*input_resolution;

  return static_cast<uint32_t>(stripe);
}

double Pixel::RA() {
  double ra, dec;
  AngularCoordinate::SurveyToEquatorial(Lambda(), Eta(), ra, dec);

  return ra;
}

double Pixel::DEC() {
  double ra, dec;
  AngularCoordinate::SurveyToEquatorial(Lambda(), Eta(), ra, dec);

  return dec;
}

double Pixel::GalLon() {
  double gal_lon, gal_lat;
  AngularCoordinate::SurveyToGalactic(Lambda(), Eta(), gal_lon, gal_lat);

  return gal_lon;
}

double Pixel::GalLat() {
  double gal_lon, gal_lat;
  AngularCoordinate::SurveyToGalactic(Lambda(), Eta(), gal_lon, gal_lat);

  return gal_lat;
}

void Pixel::Ang(AngularCoordinate& ang) {
  ang.SetSurveyCoordinates(90.0 - RadToDeg*
			   acos(1.0-2.0*(y_+0.5)/(Ny0*Resolution())),
			   RadToDeg*(2.0*Pi*(x_+0.5))/
			   (Nx0*Resolution()) + EtaOffSet);
}

double Pixel::Lambda() {
  return 90.0 -
    RadToDeg*acos(1.0 - 2.0*(y_+0.5)/(Ny0*Resolution()));
}

double Pixel::Eta() {
  double eta =
    RadToDeg*2.0*Pi*(x_+0.5)/(Nx0*Resolution()) +
    EtaOffSet;
  return (eta >= 180.0 ? eta - 360.0 : eta);
}

double Pixel::UnitSphereX() {
  return -1.0*sin(Lambda()*DegToRad);
}

double Pixel::UnitSphereY() {
  return cos(Lambda()*DegToRad)*
    cos(Eta()*DegToRad+EtaPole);
}

double Pixel::UnitSphereZ() {
  return cos(Lambda()*DegToRad)*
    sin(Eta()*DegToRad+EtaPole);
}

double Pixel::LambdaMin() {
  return 90.0 -
    RadToDeg*acos(1.0 - 2.0*(y_+1)/(Ny0*Resolution()));
}

double Pixel::LambdaMax() {
  return 90.0 -
    RadToDeg*acos(1.0 - 2.0*y_/(Ny0*Resolution()));
}

double Pixel::EtaMin() {
  double etamin =
    RadToDeg*2.0*Pi*x_/(Nx0*Resolution()) +
    EtaOffSet;
  return (etamin >= 180.0 ? etamin - 360.0 : etamin);
}

double Pixel::EtaMax() {
  double etamax =
    RadToDeg*2.0*Pi*(x_+1)/(Nx0*Resolution()) +
    EtaOffSet;
  return (etamax >= 180.0 ? etamax - 360.0 : etamax);
}

double Pixel::EtaMaxContinuous() {
  double etamax = EtaMax();
  double etamin = EtaMin();
  return (etamax > etamin ? etamax : etamax + 360.0);
}

bool Pixel::SurveyContinuous() {
  return (EtaMaxContinuous() > 180.0 ? false : true);
}

double Pixel::DECMin() {
  double ra = 0.0, dec = 0.0;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMin(), ra, dec);
  double dec_min = dec;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMax(), ra, dec);
  if (dec < dec_min) dec_min = dec;

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMin(), ra, dec);
  if (dec < dec_min) dec_min = dec;

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMax(), ra, dec);
  if (dec < dec_min) dec_min = dec;

  return dec_min;
}

double Pixel::DECMax() {
  double ra = 0.0, dec = 0.0;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMin(), ra, dec);
  double dec_max = dec;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMax(), ra, dec);
  if (dec > dec_max) dec_max = dec;

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMin(), ra, dec);
  if (dec > dec_max) dec_max = dec;

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMax(), ra, dec);
  if (dec > dec_max) dec_max = dec;

  return dec_max;
}

double Pixel::RAMin() {
  bool crosses_meridian = false;
  double ra = 0.0, dec = 0.0;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMin(), ra, dec);
  double ra_min = ra;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMax(), ra, dec);
  if ((ra > 300.0 && ra_min < 60.0) || (ra < 60.0 && ra_min > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra < ra_min) ra_min = ra;
  } else {
    if (ra > 300.0) ra_min = ra;
  }

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMin(), ra, dec);
  if ((ra > 300.0 && ra_min < 60.0) || (ra < 60.0 && ra_min > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra < ra_min) ra_min = ra;
  } else {
    if (ra_min < 60.0) ra_min = ra;
    if (ra > 300.0 && ra < ra_min) ra_min = ra;
  }

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMax(), ra, dec);
  if ((ra > 300.0 && ra_min < 60.0) || (ra < 60.0 && ra_min > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra < ra_min) ra_min = ra;
  } else {
    if (ra_min < 60.0) ra_min = ra;
    if (ra > 300.0 && ra < ra_min) ra_min = ra;
  }

  return ra_min;
}

double Pixel::RAMax() {
  bool crosses_meridian = false;
  double ra = 0.0, dec = 0.0;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMin(), ra, dec);
  double ra_max = ra;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMax(), ra, dec);
  if ((ra > 300.0 && ra_max < 60.0) || (ra < 60.0 && ra_max > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra > ra_max) ra_max = ra;
  } else {
    if (ra < 60.0) ra_max = ra;
  }

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMin(), ra, dec);
  if ((ra > 300.0 && ra_max < 60.0) || (ra < 60.0 && ra_max > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra > ra_max) ra_max = ra;
  } else {
    if (ra_max > 300.0) ra_max = ra;
    if (ra < 60.0 && ra > ra_max) ra_max = ra;
  }

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMax(), ra, dec);
  if ((ra > 300.0 && ra_max < 60.0) || (ra < 60.0 && ra_max > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra > ra_max) ra_max = ra;
  } else {
    if (ra_max > 300.0) ra_max = ra;
    if (ra < 60.0 && ra > ra_max) ra_max = ra;
  }

  return ra_max;
}

double Pixel::RAMaxContinuous() {
  bool crosses_meridian = false;
  double ra = 0.0, dec = 0.0;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMin(), ra, dec);
  double ra_max = ra;

  AngularCoordinate::SurveyToEquatorial(LambdaMin(), EtaMax(), ra, dec);
  if ((ra > 300.0 && ra_max < 60.0) || (ra < 60.0 && ra_max > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra > ra_max) ra_max = ra;
  } else {
    if (ra < 60.0) ra_max = ra;
  }

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMin(), ra, dec);
  if ((ra > 300.0 && ra_max < 60.0) || (ra < 60.0 && ra_max > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra > ra_max) ra_max = ra;
  } else {
    if (ra_max > 300.0) ra_max = ra;
    if (ra < 60.0 && ra > ra_max) ra_max = ra;
  }

  AngularCoordinate::SurveyToEquatorial(LambdaMax(), EtaMax(), ra, dec);
  if ((ra > 300.0 && ra_max < 60.0) || (ra < 60.0 && ra_max > 300.0))
    crosses_meridian = true;
  if (!crosses_meridian) {
    if (ra > ra_max) ra_max = ra;
  } else {
    if (ra_max > 300.0) ra_max = ra;
    if (ra < 60.0 && ra > ra_max) ra_max = ra;
  }

  if (crosses_meridian && ra_max < 60.0) ra_max += 360.0;

  return ra_max;
}

bool Pixel::EquatorialContinuous() {
  return (RAMaxContinuous() > 360.0 ? false : true);
}

double Pixel::GalLatMin() {
  double gal_lon = 0.0, gal_lat = 0.0;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMin(), gal_lon, gal_lat);
  double gal_lat_min = gal_lat;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMax(), gal_lon, gal_lat);
  if (gal_lat < gal_lat_min) gal_lat_min = gal_lat;

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMin(), gal_lon, gal_lat);
  if (gal_lat < gal_lat_min) gal_lat_min = gal_lat;

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMax(), gal_lon, gal_lat);
  if (gal_lat < gal_lat_min) gal_lat_min = gal_lat;

  return gal_lat_min;
}

double Pixel::GalLatMax() {
  double gal_lon = 0.0, gal_lat = 0.0;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMin(), gal_lon, gal_lat);
  double gal_lat_max = gal_lat;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMax(), gal_lon, gal_lat);
  if (gal_lat > gal_lat_max) gal_lat_max = gal_lat;

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMin(), gal_lon, gal_lat);
  if (gal_lat > gal_lat_max) gal_lat_max = gal_lat;

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMax(), gal_lon, gal_lat);
  if (gal_lat > gal_lat_max) gal_lat_max = gal_lat;

  return gal_lat_max;
}

double Pixel::GalLonMin() {
  bool crosses_meridian = false;
  double gal_lon = 0.0, gal_lat = 0.0;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMin(), gal_lon, gal_lat);
  double gal_lon_min = gal_lon;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMax(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_min < 60.0) ||
      (gal_lon < 60.0 && gal_lon_min > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon < gal_lon_min) gal_lon_min = gal_lon;
  } else {
    if (gal_lon > 300.0) gal_lon_min = gal_lon;
  }

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMin(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_min < 60.0) ||
      (gal_lon < 60.0 && gal_lon_min > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon < gal_lon_min) gal_lon_min = gal_lon;
  } else {
    if (gal_lon_min < 60.0) gal_lon_min = gal_lon;
    if (gal_lon > 300.0 && gal_lon < gal_lon_min) gal_lon_min = gal_lon;
  }

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMax(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_min < 60.0) ||
      (gal_lon < 60.0 && gal_lon_min > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon < gal_lon_min) gal_lon_min = gal_lon;
  } else {
    if (gal_lon_min < 60.0) gal_lon_min = gal_lon;
    if (gal_lon > 300.0 && gal_lon < gal_lon_min) gal_lon_min = gal_lon;
  }

  return gal_lon_min;
}

double Pixel::GalLonMax() {
  bool crosses_meridian = false;
  double gal_lon = 0.0, gal_lat = 0.0;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMin(), gal_lon, gal_lat);
  double gal_lon_max = gal_lon;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMax(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_max < 60.0) ||
      (gal_lon < 60.0 && gal_lon_max > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  } else {
    if (gal_lon < 60.0) gal_lon_max = gal_lon;
  }

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMin(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_max < 60.0) ||
      (gal_lon < 60.0 && gal_lon_max > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  } else {
    if (gal_lon_max > 300.0) gal_lon_max = gal_lon;
    if (gal_lon < 60.0 && gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  }

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMax(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_max < 60.0) ||
      (gal_lon < 60.0 && gal_lon_max > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  } else {
    if (gal_lon_max > 300.0) gal_lon_max = gal_lon;
    if (gal_lon < 60.0 && gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  }

  return gal_lon_max;
}

double Pixel::GalLonMaxContinuous() {
  bool crosses_meridian = false;
  double gal_lon = 0.0, gal_lat = 0.0;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMin(), gal_lon, gal_lat);
  double gal_lon_max = gal_lon;

  AngularCoordinate::SurveyToGalactic(LambdaMin(), EtaMax(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_max < 60.0) ||
      (gal_lon < 60.0 && gal_lon_max > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  } else {
    if (gal_lon < 60.0) gal_lon_max = gal_lon;
  }

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMin(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_max < 60.0) ||
      (gal_lon < 60.0 && gal_lon_max > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  } else {
    if (gal_lon_max > 300.0) gal_lon_max = gal_lon;
    if (gal_lon < 60.0 && gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  }

  AngularCoordinate::SurveyToGalactic(LambdaMax(), EtaMax(), gal_lon, gal_lat);
  if ((gal_lon > 300.0 && gal_lon_max < 60.0) ||
      (gal_lon < 60.0 && gal_lon_max > 300.0)) crosses_meridian = true;
  if (!crosses_meridian) {
    if (gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  } else {
    if (gal_lon_max > 300.0) gal_lon_max = gal_lon;
    if (gal_lon < 60.0 && gal_lon > gal_lon_max) gal_lon_max = gal_lon;
  }

  if (crosses_meridian) gal_lon_max += 360.0;

  return gal_lon_max;
}

bool Pixel::GalacticContinuous() {
  return (GalLonMaxContinuous() > 360.0 ? false : true);
}

bool Pixel::ContinuousBounds(AngularCoordinate::Sphere sphere) {
  return (sphere == AngularCoordinate::Survey ? SurveyContinuous() :
	  (sphere == AngularCoordinate::Equatorial ? EquatorialContinuous() :
	   GalacticContinuous()));
}

double Pixel::UnitSphereX_UL() {
  return -1.0*sin(LambdaMax()*DegToRad);
}

double Pixel::UnitSphereY_UL() {
  return cos(LambdaMax()*DegToRad)*
    cos(EtaMin()*DegToRad+EtaPole);
}

double Pixel::UnitSphereZ_UL() {
  return cos(LambdaMax()*DegToRad)*
    sin(EtaMin()*DegToRad+EtaPole);
}

double Pixel::UnitSphereX_UR() {
  return -1.0*sin(LambdaMax()*DegToRad);
}

double Pixel::UnitSphereY_UR() {
  return cos(LambdaMax()*DegToRad)*
    cos(EtaMax()*DegToRad+EtaPole);
}

double Pixel::UnitSphereZ_UR() {
  return cos(LambdaMax()*DegToRad)*
    sin(EtaMax()*DegToRad+EtaPole);
}

double Pixel::UnitSphereX_LL() {
  return -1.0*sin(LambdaMin()*DegToRad);
}

double Pixel::UnitSphereY_LL() {
  return cos(LambdaMin()*DegToRad)*
    cos(EtaMin()*DegToRad+EtaPole);
}

double Pixel::UnitSphereZ_LL() {
  return cos(LambdaMin()*DegToRad)*
    sin(EtaMin()*DegToRad+EtaPole);
}

double Pixel::UnitSphereX_LR() {
  return -1.0*sin(LambdaMin()*DegToRad);
}

double Pixel::UnitSphereY_LR() {
  return cos(LambdaMin()*DegToRad)*
    cos(EtaMax()*DegToRad+EtaPole);
}

double Pixel::UnitSphereZ_LR() {
  return cos(LambdaMin()*DegToRad)*
    sin(EtaMax()*DegToRad+EtaPole);
}

void Pixel::Iterate(bool wrap_pixel) {
  if (x_ == Nx0*Resolution() - 1) {
    x_ = 0;
    if (!wrap_pixel) y_++;
  } else {
    x_++;
  }
}

uint32_t Pixel::PixelX0() {
  return static_cast<uint32_t>(x_*HPixResolution/Resolution())*
    Resolution()/HPixResolution;
}

uint32_t Pixel::PixelY0() {
  return static_cast<uint32_t>(y_*HPixResolution/Resolution())*
    Resolution()/HPixResolution;
}

uint32_t Pixel::PixelX1() {
  return PixelX0() + Resolution()/HPixResolution;
}

uint32_t Pixel::PixelY1() {
  return PixelY0() + Resolution()/HPixResolution;
}

void Pixel::GenerateRandomPoints(AngularVector& ang, uint32_t n_point) {
  if (!ang.empty()) ang.clear();
  ang.reserve(n_point);

  MTRand mtrand;

  mtrand.seed();

  double z_min = sin(LambdaMin()*DegToRad);
  double z_max = sin(LambdaMax()*DegToRad);
  double eta_min = EtaMin();
  double eta_max = EtaMax();
  double z = 0.0;
  double lambda = 0.0;
  double eta = 0.0;

  for (uint32_t m=0;m<n_point;m++) {
    z = z_min + mtrand.rand(z_max - z_min);
    lambda = asin(z)*RadToDeg;
    eta = eta_min + mtrand.rand(eta_max - eta_min);

    AngularCoordinate tmp_ang;
    tmp_ang.SetSurveyCoordinates(lambda, eta);
    ang.push_back(tmp_ang);
  }
}

uint8_t Pixel::Resolution2Level(const uint32_t resolution) {
  return MostSignificantBit(resolution);
}

uint32_t Pixel::Level2Resolution(const uint8_t level) {
  return static_cast<uint32_t>(1 << level);
}

void Pixel::Ang2Pix(uint32_t input_resolution, AngularCoordinate& ang,
		    uint32_t& output_pixnum) {
  double lambda = ang.Lambda();
  double eta = ang.Eta();
  uint32_t nx = Nx0*input_resolution;
  uint32_t ny = Ny0*input_resolution;

  eta -= EtaOffSet;

  eta *= DegToRad;

  if (eta <= 0.0) eta += 2.0*Pi;

  eta /= 2.0*Pi;
  uint32_t i = static_cast<uint32_t>(nx*eta);

  lambda = (90.0 - lambda)*DegToRad;

  uint32_t j;
  if (lambda >= Pi) {
    j = ny - 1;
  } else {
    j = static_cast<uint32_t>(ny*((1.0 - cos(lambda))/2.0));
  }

  output_pixnum = nx*j + i;
}

void Pixel::Pix2Ang(uint32_t input_resolution, uint32_t input_pixnum,
		    AngularCoordinate& ang) {
  uint32_t nx = Nx0*input_resolution;
  uint32_t ny = Ny0*input_resolution;

  uint32_t y = input_pixnum/nx;
  uint32_t x = input_pixnum - nx*y;

  ang.SetSurveyCoordinates(90.0 - RadToDeg*acos(1.0-2.0*(y+0.5)/ny),
			   RadToDeg*(2.0*Pi*(x+0.5))/nx +
			   EtaOffSet);
}

void Pixel::Pix2HPix(uint32_t input_resolution,
                     uint32_t input_pixnum,
                     uint32_t& output_hpixnum,
                     uint32_t& output_superpixnum) {
  uint32_t nx = Nx0*input_resolution;

  uint32_t y = input_pixnum/nx;
  uint32_t x = input_pixnum - nx*y;

  uint32_t hnx = input_resolution/HPixResolution;

  uint32_t x0 = x/hnx;
  uint32_t y0 = y/hnx;

  x -= x0*hnx;
  y -= y0*hnx;

  output_hpixnum = hnx*y + x;
  output_superpixnum = Nx0*HPixResolution*y0 + x0;
}

void Pixel::HPix2Pix(uint32_t input_resolution,
                     uint32_t input_hpixnum,
                     uint32_t input_superpixnum,
                     uint32_t& output_pixnum) {
  uint32_t nx = Nx0*input_resolution;
  uint32_t hnx = input_resolution/HPixResolution;

  uint32_t y0 = input_superpixnum/(Nx0*HPixResolution);
  uint32_t x0 = input_superpixnum - y0*Nx0*HPixResolution;

  y0 *= input_resolution/HPixResolution;
  x0 *= input_resolution/HPixResolution;

  uint32_t y = input_hpixnum/hnx;
  uint32_t x = input_hpixnum - hnx*y;

  output_pixnum = nx*(y+y0) + x+x0;
}

void Pixel::SuperPix(uint32_t hi_resolution, uint32_t hi_pixnum,
		     uint32_t lo_resolution, uint32_t& lo_pixnum) {
  if (hi_resolution < lo_resolution) {
    std::cout << "Can't go from low resolution to higher resolution.\n ";
    exit(1);
  } else {
    uint32_t nx_hi = Nx0*hi_resolution;
    uint32_t nx_lo = Nx0*lo_resolution;

    uint32_t ratio = hi_resolution/lo_resolution;

    uint32_t j = hi_pixnum/nx_hi;
    uint32_t i = hi_pixnum - nx_hi*j;

    i /= ratio;
    j /= ratio;

    lo_pixnum = nx_lo*j + i;
  }
}

void Pixel::NextSubPix(uint32_t input_resolution, uint32_t input_pixnum,
		       uint32_t& sub_pixnum1,
		       uint32_t& sub_pixnum2,
		       uint32_t& sub_pixnum3,
		       uint32_t& sub_pixnum4) {
  uint32_t nx_hi = 2*Nx0*input_resolution;
  uint32_t nx_lo = Nx0*input_resolution;

  uint32_t j = input_pixnum/nx_lo;
  uint32_t i = input_pixnum - nx_lo*j;

  sub_pixnum1 = nx_hi*(2*j) + 2*i;
  sub_pixnum2 = nx_hi*(2*j) + 2*i + 1;
  sub_pixnum3 = nx_hi*(2*j + 1) + 2*i;
  sub_pixnum4 = nx_hi*(2*j + 1) + 2*i + 1;
}

void Pixel::SubPix(uint32_t lo_resolution, uint32_t lo_pixnum,
		   uint32_t hi_resolution, uint32_t& x_min,
		   uint32_t& x_max, uint32_t& y_min,
		   uint32_t& y_max) {
  uint32_t nx_hi = Nx0*hi_resolution;
  if (lo_resolution == hi_resolution) {
    y_min = lo_pixnum/nx_hi;
    y_max = lo_pixnum/nx_hi;
    x_min = lo_pixnum - nx_hi*y_min;
    x_max = lo_pixnum - nx_hi*y_max;
  } else {
    uint32_t tmp_pixnum, pixnum1, pixnum2, pixnum3, pixnum4;
    uint32_t tmp_res;

    tmp_pixnum = lo_pixnum;
    for (tmp_res=lo_resolution;tmp_res<hi_resolution;tmp_res*=2) {
      NextSubPix(tmp_res, tmp_pixnum, pixnum1, pixnum2, pixnum3, pixnum4);
      tmp_pixnum = pixnum1;
    }

    y_min = tmp_pixnum/nx_hi;
    x_min = tmp_pixnum - nx_hi*y_min;

    tmp_pixnum = lo_pixnum;
    for (tmp_res=lo_resolution;tmp_res<hi_resolution;tmp_res*=2) {
      NextSubPix(tmp_res, tmp_pixnum, pixnum1, pixnum2, pixnum3, pixnum4);
      tmp_pixnum = pixnum4;
    }

    y_max = tmp_pixnum/nx_hi;
    x_max = tmp_pixnum - nx_hi*y_max;
  }
}

void Pixel::AreaIndex(uint32_t input_resolution,
		      double lammin, double lammax,
		      double etamin, double etamax,
		      uint32_t& x_min, uint32_t& x_max,
		      uint32_t& y_min, uint32_t& y_max) {
  uint32_t nx = Nx0*input_resolution;
  etamin -= EtaOffSet;
  etamin *= DegToRad;
  if (etamin <= 0.0) etamin = etamin + 2.0*Pi;
  etamin /= 2.0*Pi;
  x_min = static_cast<uint32_t>(nx*etamin);

  etamax -= EtaOffSet;
  etamax *= DegToRad;
  if (etamax <= 0.0) etamax = etamax + 2.0*Pi;
  etamax /= 2.0*Pi;
  x_max = static_cast<uint32_t>(nx*etamax);

  uint32_t ny = Ny0*input_resolution;
  lammax = (90.0 - lammax)*DegToRad;
  if (lammax >= Pi) {
    y_min = ny - 1;
  } else {
    y_min = static_cast<uint32_t>(ny*((1.0 - cos(lammax))/2.0));
  }

  lammin = (90.0 - lammin)*DegToRad;
  if (lammin >= Pi) {
    y_max = ny - 1;
  } else {
    y_max = static_cast<uint32_t>(ny*((1.0 - cos(lammin))/2.0));
  }
}

uint8_t Pixel::Pix2EtaStep(uint32_t input_resolution, uint32_t input_pixnum,
		       double theta) {
  uint32_t nx = Nx0*input_resolution;
  uint32_t ny = Ny0*input_resolution;

  uint32_t y = input_pixnum/nx;
  double lam = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.5)/ny);

  double deta = 2.5/(input_resolution/4);
  double eta_step = theta;
  eta_step *= 1.0 +
    lam*lam*(0.000192312 - lam*lam*(1.82764e-08 - 1.28162e-11*lam*lam));
  uint8_t etastep = 1;
  while (eta_step > etastep*deta) etastep++;

  return etastep;
}

void Pixel::PixelBound(uint32_t input_resolution, uint32_t input_pixnum,
		       double& lammin, double& lammax, double& etamin,
		       double& etamax) {
  uint32_t nx = Nx0*input_resolution;
  uint32_t ny = Ny0*input_resolution;

  uint32_t y = input_pixnum/nx;
  uint32_t x = input_pixnum - nx*y;

  lammin = 90.0 - RadToDeg*acos(1.0 - 2.0*(y+1)/ny);
  lammax = 90.0 - RadToDeg*acos(1.0 - 2.0*y/ny);
  etamin = RadToDeg*2.0*Pi*(x+0.0)/nx + EtaOffSet;
  if (etamin >= 180.0) etamin = etamin - 360.0;
  etamax = RadToDeg*2.0*Pi*(x+1.0)/nx + EtaOffSet;
  if (etamax >= 180.0) etamax = etamax - 360.0;
}

void Pixel::CohortPix(uint32_t input_resolution, uint32_t input_pixnum,
		      uint32_t& co_pixnum1,
		      uint32_t& co_pixnum2,
		      uint32_t& co_pixnum3) {
  uint32_t tmp_pixnum, pixnum1, pixnum2, pixnum3, pixnum4;

  SuperPix(input_resolution, input_pixnum, input_resolution/2, tmp_pixnum);

  NextSubPix(input_resolution/2, tmp_pixnum,
	     pixnum1, pixnum2, pixnum3, pixnum4);

  if (input_pixnum == pixnum1) {
    co_pixnum1 = pixnum2;
    co_pixnum2 = pixnum3;
    co_pixnum3 = pixnum4;
  }
  if (input_pixnum == pixnum2) {
    co_pixnum1 = pixnum1;
    co_pixnum2 = pixnum3;
    co_pixnum3 = pixnum4;
  }
  if (input_pixnum == pixnum3) {
    co_pixnum1 = pixnum1;
    co_pixnum2 = pixnum2;
    co_pixnum3 = pixnum4;
  }
  if (input_pixnum == pixnum4) {
    co_pixnum1 = pixnum1;
    co_pixnum2 = pixnum2;
    co_pixnum3 = pixnum3;
  }
}

double Pixel::PixelArea(uint32_t resolution) {
  return HPixArea*HPixResolution*HPixResolution/
    (resolution*resolution);
}

double Pixel::HPixelArea(uint32_t resolution) {
  return HPixArea*HPixResolution*HPixResolution/
    (resolution*resolution);
}

void Pixel::XY2Pix(uint32_t resolution, uint32_t x, uint32_t y,
		   uint32_t& pixnum) {
  pixnum = Nx0*resolution*y + x;
}

void Pixel::Pix2XY(uint32_t resolution, uint32_t pixnum,
		   uint32_t& x, uint32_t& y) {
  y = pixnum/(Nx0*resolution);
  x = pixnum - Nx0*resolution*y;
}

bool Pixel::LocalOrder(Pixel pix_a, Pixel pix_b) {
  if (pix_a.Resolution() == pix_b.Resolution()) {
    if (pix_a.PixelY() == pix_b.PixelY()) {
      return (pix_a.PixelX() < pix_b.PixelX() ? true : false);
    } else {
      return (pix_a.PixelY() < pix_b.PixelY() ? true : false);
    }
  } else {
    return (pix_a.Resolution() < pix_b.Resolution() ? true : false);
  }
}

bool Pixel::SuperPixelBasedOrder(Pixel pix_a, Pixel pix_b) {
  if (pix_a.Superpixnum() == pix_b.Superpixnum()) {
    if (pix_a.Resolution() == pix_b.Resolution()) {
      return (pix_a.HPixnum() < pix_b.HPixnum() ? true : false);
    } else {
      return (pix_a.Resolution() < pix_b.Resolution() ? true : false);
    }
  } else {
    return (pix_a.Superpixnum() < pix_b.Superpixnum() ? true : false);
  }
}

bool Pixel::SuperPixelOrder(Pixel pix_a, Pixel pix_b) {
  return (pix_a.Superpixnum() < pix_b.Superpixnum() ? true : false);
}

bool Pixel::WeightedOrder(Pixel pix_a, Pixel pix_b) {
  return (pix_a.Weight() < pix_b.Weight() ? true : false);
}

bool Pixel::WeightMatch(Pixel& pix_a, Pixel& pix_b) {
  return ((pix_b.Weight() < pix_a.Weight() + 0.000001) &&
	  (pix_b.Weight() > pix_a.Weight() - 0.000001) ? true : false);
}

bool Pixel::WeightedPixelMatch(Pixel& pix_a, Pixel& pix_b) {
  return ((pix_a.Resolution() == pix_b.Resolution()) &&
	  (pix_a.PixelX() == pix_b.PixelX()) &&
	  (pix_a.PixelY() == pix_b.PixelY()) &&
	  (pix_b.Weight() < pix_a.Weight() + 0.000001) &&
	  (pix_b.Weight() > pix_a.Weight() - 0.000001) ? true : false);
}

bool Pixel::PixelMatch(Pixel& pix_a, Pixel& pix_b) {
  return ((pix_a.Resolution() == pix_b.Resolution()) &&
	  (pix_a.PixelX() == pix_b.PixelX()) &&
	  (pix_a.PixelY() == pix_b.PixelY()) ? true : false);
}

void Pixel::ResolvePixel(PixelVector& pix, bool ignore_weight) {
  sort(pix.begin(),pix.end(),Pixel::SuperPixelBasedOrder);

  PixelVector tmp_pix;
  PixelVector final_pix;
  uint32_t superpixnum = pix[0].Superpixnum();

  for (PixelIterator iter=pix.begin();iter!=pix.end();++iter) {
    if (superpixnum == iter->Superpixnum()) {
      tmp_pix.push_back(*iter);
    } else {
      ResolveSuperPixel(tmp_pix,ignore_weight);

      for (uint32_t i=0;i<tmp_pix.size();i++)
        final_pix.push_back(tmp_pix[i]);

      tmp_pix.clear();
      superpixnum = iter->Superpixnum();
      tmp_pix.push_back(*iter);
    }
  }

  ResolveSuperPixel(tmp_pix,ignore_weight);

  for (uint32_t i=0;i<tmp_pix.size();i++)
    final_pix.push_back(tmp_pix[i]);

  tmp_pix.clear();

  pix.clear();
  pix.reserve(final_pix.size());

  for (uint32_t i=0;i<final_pix.size();i++) pix.push_back(final_pix[i]);

  final_pix.clear();
}

void Pixel::ResolveSuperPixel(PixelVector& pix, bool ignore_weight) {
  if (ignore_weight)
    for (uint32_t i=0;i<pix.size();i++) pix[i].SetWeight(1.0);

  // First, remove any duplicate pixels or pixels that are contained within
  // other pixels.
  PixelVector unique_pix;
  Pixel::FindUniquePixels(pix, unique_pix);

  uint32_t n_start = pix.size();
  pix.clear();
  pix.reserve(unique_pix.size());

  for (uint32_t i=0;i<unique_pix.size();i++) pix.push_back(unique_pix[i]);

  sort(pix.begin(), pix.end(), Pixel::SuperPixelBasedOrder);
  unique_pix.clear();
  uint32_t n_finish = unique_pix.size();

  // Now we iterate through the pixel list, looking for possible cases where
  // we might combine high resolution pixels into low resolution pixels.  In
  // order to do this, we need to find all of the cohort pixels as well as
  // verify that they all have the same Weight.
  while (n_start != n_finish) {
    n_start = pix.size();

    unique_pix.reserve(pix.size());
    PixelIterator search_begin = pix.begin();
    for (uint32_t i=0;i<pix.size();i++) {
      if ((pix[i].Resolution() > HPixResolution) &&
	  pix[i].FirstCohort()) {
        bool found_cohort = false;
        Pixel pix_a;
        Pixel pix_b;
        Pixel pix_c;

	pix[i].CohortPix(pix_a,pix_b,pix_c);

	PixelPair iter = equal_range(search_begin,pix.end(),pix_a,
				     Pixel::SuperPixelBasedOrder);
	if ((iter.first != iter.second) &&
	    (Pixel::WeightMatch(pix[i],*iter.first))) {
	  found_cohort = true;
	} else {
	  found_cohort = false;
	}

	if (found_cohort) {
          iter = equal_range(search_begin,pix.end(),pix_b,
                             Pixel::SuperPixelBasedOrder);
          if ((iter.first != iter.second) &&
              (Pixel::WeightMatch(pix[i],*iter.first))) {
	    found_cohort = true;
	  } else {
	    found_cohort = false;
	  }
	}

	if (found_cohort) {
          iter = equal_range(search_begin,pix.end(),pix_c,
                             Pixel::SuperPixelBasedOrder);
          if ((iter.first != iter.second) &&
              (Pixel::WeightMatch(pix[i],*iter.first))) {
	    found_cohort = true;
	  } else {
	    found_cohort = false;
	  }
	}

	if (found_cohort) {
	  pix_a = pix[i];
	  pix_a.SetToSuperPix(pix_a.Resolution()/2);
	  unique_pix.push_back(pix_a);
	} else {
	  unique_pix.push_back(pix[i]);
	}
      } else {
	unique_pix.push_back(pix[i]);
      }
      ++search_begin;
    }

    if (unique_pix.size() != pix.size()) {
      std::cout <<
	"Something has gone wrong searching for superpixels. Exiting.\n";
      exit(1);
    }

    for (uint32_t i=0;i<unique_pix.size();i++) pix[i] = unique_pix[i];
    unique_pix.clear();

    Pixel::FindUniquePixels(pix, unique_pix);

    n_finish = unique_pix.size();

    pix.clear();
    pix.reserve(unique_pix.size());
    for (uint32_t i=0;i<unique_pix.size();i++)
      pix.push_back(unique_pix[i]);

    unique_pix.clear();
  }
}

void Pixel::FindUniquePixels(PixelVector& input_pix, PixelVector& unique_pix) {
  if (!unique_pix.empty()) unique_pix.clear();
  sort(input_pix.begin(), input_pix.end(), Pixel::SuperPixelBasedOrder);

  PixelIterator search_end = input_pix.begin();
  unique_pix.push_back(input_pix[0]);
  ++search_end;
  // First, we iterate through the pixels, checking to see if copies of a given
  // pixel exist (we take the first such instance) or if superpixels which
  // contain the given pixel exist (we always take larger pixels).  If neither
  // case is met, then we keep the pixel.
  for (uint32_t i=1;i<input_pix.size();i++) {
    bool keep_pixel = true;

    // Check against copies of the current pixel.
    if (Pixel::WeightedPixelMatch(input_pix[i], input_pix[i-1]))
      keep_pixel = false;

    if ((keep_pixel) && (input_pix[i].Resolution() > HPixResolution)) {
      Pixel tmp_pix = input_pix[i];

      // Check for larger pixels which might contain the current pixel.
      while (tmp_pix.Resolution() > HPixResolution) {
        tmp_pix.SetToSuperPix(tmp_pix.Resolution()/2);

        if (binary_search(input_pix.begin(),search_end,tmp_pix,
                          Pixel::SuperPixelBasedOrder))
          keep_pixel = false;
      }
    }

    if (keep_pixel) unique_pix.push_back(input_pix[i]);
    ++search_end;
  }
}

void Pixel::Ang2HPix(uint32_t input_resolution, AngularCoordinate& ang,
		     uint32_t& output_hpixnum,
		     uint32_t& output_superpixnum) {
  uint32_t nx = Nx0*input_resolution;
  uint32_t ny = Ny0*input_resolution;

  uint32_t hnx = input_resolution/HPixResolution;

  double eta = (ang.Eta() - EtaOffSet)*DegToRad;

  if (eta <= 0.0) eta += 2.0*Pi;

  eta /= 2.0*Pi;
  uint32_t x = static_cast<uint32_t>(nx*eta);

  double lambda = (90.0 - ang.Lambda())*DegToRad;

  uint32_t y;
  if (lambda >= Pi) {
    y = ny - 1;
  } else {
    y = static_cast<uint32_t>(ny*((1.0 - cos(lambda))/2.0));
  }

  uint32_t x0 = x/hnx;
  uint32_t y0 = y/hnx;

  x -= x0*hnx;
  y -= y0*hnx;

  output_hpixnum = nx*y + x;
  output_superpixnum = Nx0*HPixResolution*y0 + x0;
}

void Pixel::HPix2Ang(uint32_t input_resolution, uint32_t input_hpixnum,
		     uint32_t input_superpixnum,
		     AngularCoordinate& ang) {
  uint32_t nx = Nx0*input_resolution;
  uint32_t ny = Ny0*input_resolution;

  uint32_t hnx = input_resolution/HPixResolution;

  uint32_t y0 = input_superpixnum/(Nx0*HPixResolution);
  uint32_t x0 = input_superpixnum - y0*Nx0*HPixResolution;

  y0 *= hnx;
  x0 *= hnx;

  uint32_t y = input_hpixnum/hnx;
  uint32_t x = input_hpixnum - hnx*y;

  ang.SetSurveyCoordinates(90.0-RadToDeg*acos(1.0-2.0*(y+y0+0.5)/ny),
			   RadToDeg*(2.0*Pi*(x+x0+0.5))/nx +
			   EtaOffSet);
}

void Pixel::XY2HPix(uint32_t input_resolution, uint32_t x,
		    uint32_t y, uint32_t& output_hpixnum,
		    uint32_t& output_superpixnum) {
  uint32_t hnx = input_resolution/HPixResolution;

  uint32_t x0 = x/hnx;
  uint32_t y0 = y/hnx;

  x -= x0;
  y -= y0;

  output_hpixnum = hnx*y + x;
  output_superpixnum = Nx0*HPixResolution*y0 + x0;
}

void Pixel::HPix2XY(uint32_t input_resolution, uint32_t input_hpixnum,
		    uint32_t input_superpixnum, uint32_t& x, uint32_t& y) {
  uint32_t hnx = input_resolution/HPixResolution;

  uint32_t y0 = input_superpixnum/(Nx0*HPixResolution);
  uint32_t x0 = input_superpixnum - y0*Nx0*HPixResolution;

  uint32_t tmp_y = input_hpixnum/hnx;
  uint32_t tmp_x = input_hpixnum - hnx*tmp_y;

  x = tmp_x + x0*hnx;
  y = tmp_y + y0*hnx;
}

void Pixel::SuperHPix(uint32_t hi_resolution, uint32_t hi_hpixnum,
		      uint32_t lo_resolution, uint32_t& lo_hpixnum) {
  if (hi_resolution < lo_resolution) {
    std::cout << "Can't go from low resolution to higher resolution.\n ";
    exit(1);
  } else {
    uint32_t nx_hi = hi_resolution/HPixResolution;
    uint32_t nx_lo = lo_resolution/HPixResolution;

    uint32_t ratio = hi_resolution/lo_resolution;

    uint32_t y = hi_hpixnum/nx_hi;
    uint32_t x = hi_hpixnum - nx_hi*y;

    x /= ratio;
    y /= ratio;

    lo_hpixnum = nx_lo*y + x;
  }
}

void Pixel::NextSubHPix(uint32_t input_resolution, uint32_t input_hpixnum,
			uint32_t& sub_hpixnum1,
			uint32_t& sub_hpixnum2,
			uint32_t& sub_hpixnum3,
			uint32_t& sub_hpixnum4) {
  uint32_t nx_hi = 2*input_resolution/HPixResolution;
  uint32_t nx_lo = input_resolution/HPixResolution;

  uint32_t y = input_hpixnum/nx_lo;
  uint32_t x = input_hpixnum - nx_lo*y;

  sub_hpixnum1 = nx_hi*(2*y) + 2*x;
  sub_hpixnum2 = nx_hi*(2*y) + 2*x + 1;
  sub_hpixnum3 = nx_hi*(2*y + 1) + 2*x;
  sub_hpixnum4 = nx_hi*(2*y + 1) + 2*x + 1;
}

void Pixel::SubHPix(uint32_t lo_resolution, uint32_t lo_hpixnum,
		    uint32_t lo_superpixnum, uint32_t hi_resolution,
		    uint32_t& x_min, uint32_t& x_max,
		    uint32_t& y_min, uint32_t& y_max) {
  uint32_t tmp_x, tmp_y;

  if (lo_resolution == hi_resolution) {
    HPix2XY(lo_resolution,lo_hpixnum,lo_superpixnum,tmp_x,tmp_y);

    y_min = tmp_y;
    y_max = tmp_y;
    x_min = tmp_x;
    x_max = tmp_x;
  } else {
    uint32_t tmp_hpixnum, hpixnum1, hpixnum2, hpixnum3, hpixnum4;
    uint32_t tmp_res;

    tmp_hpixnum = lo_hpixnum;
    for (tmp_res=lo_resolution;tmp_res<hi_resolution;tmp_res*=2) {
      NextSubHPix(tmp_res, tmp_hpixnum, hpixnum1,
		  hpixnum2, hpixnum3, hpixnum4);
      tmp_hpixnum = hpixnum1;
    }

    HPix2XY(hi_resolution,tmp_hpixnum,lo_superpixnum,tmp_x,tmp_y);

    y_min = tmp_y;
    x_min = tmp_x;

    tmp_hpixnum = lo_hpixnum;
    for (tmp_res=lo_resolution;tmp_res<hi_resolution;tmp_res*=2) {
      NextSubHPix(tmp_res, tmp_hpixnum, hpixnum1,
		  hpixnum2, hpixnum3, hpixnum4);
      tmp_hpixnum = hpixnum4;
    }

    HPix2XY(hi_resolution,tmp_hpixnum,lo_superpixnum,tmp_x,tmp_y);

    y_max = tmp_y;
    x_max = tmp_x;
  }
}

void Pixel::HPixelBound(uint32_t input_resolution, uint32_t input_hpixnum,
			uint32_t input_superpixnum,
			double& lammin, double& lammax,
			double& etamin, double& etamax) {
  uint32_t nx = Nx0*input_resolution;
  uint32_t ny = Ny0*input_resolution;

  uint32_t hnx = input_resolution/HPixResolution;

  uint32_t y0 = input_superpixnum/(Nx0*HPixResolution);
  uint32_t x0 = input_superpixnum - y0*Nx0*HPixResolution;

  y0 *= hnx;
  x0 *= hnx;

  uint32_t y = input_hpixnum/hnx;
  uint32_t x = input_hpixnum - hnx*y;

  lammin = 90.0 - RadToDeg*acos(1.0 - 2.0*(y+y0+1)/ny);
  lammax = 90.0 - RadToDeg*acos(1.0 - 2.0*(y+y0)/ny);
  etamin =
    RadToDeg*2.0*Pi*(x+x0+0.0)/nx + EtaOffSet;
  if (etamin >= 180.0) etamin = etamin - 360.0;
  etamax =
    RadToDeg*2.0*Pi*(x+x0+1.0)/nx + EtaOffSet;
  if (etamax >= 180.0) etamax = etamax - 360.0;
}

void Pixel::CohortHPix(uint32_t input_resolution, uint32_t input_hpixnum,
		       uint32_t& co_hpixnum1,
		       uint32_t& co_hpixnum2,
		       uint32_t& co_hpixnum3) {
  uint32_t tmp_hpixnum, hpixnum1, hpixnum2, hpixnum3, hpixnum4;

  SuperHPix(input_resolution, input_hpixnum, input_resolution/2, tmp_hpixnum);

  NextSubHPix(input_resolution/2, tmp_hpixnum,
	      hpixnum1, hpixnum2, hpixnum3, hpixnum4);

  if (input_hpixnum == hpixnum1) {
    co_hpixnum1 = hpixnum2;
    co_hpixnum2 = hpixnum3;
    co_hpixnum3 = hpixnum4;
  }
  if (input_hpixnum == hpixnum2) {
    co_hpixnum1 = hpixnum1;
    co_hpixnum2 = hpixnum3;
    co_hpixnum3 = hpixnum4;
  }
  if (input_hpixnum == hpixnum3) {
    co_hpixnum1 = hpixnum1;
    co_hpixnum2 = hpixnum2;
    co_hpixnum3 = hpixnum4;
  }
  if (input_hpixnum == hpixnum4) {
    co_hpixnum1 = hpixnum1;
    co_hpixnum2 = hpixnum2;
    co_hpixnum3 = hpixnum3;
  }
}

uint8_t Pixel::HPix2EtaStep(uint32_t input_resolution, uint32_t input_hpixnum,
			    uint32_t input_superpixnum, double theta) {

  uint32_t ny = Ny0*input_resolution;
  uint32_t hnx = input_resolution/HPixResolution;

  uint32_t y0 = input_superpixnum/(Nx0*HPixResolution);
  uint32_t x0 = input_superpixnum - y0*Nx0*HPixResolution;

  y0 *= hnx;
  x0 *= hnx;

  uint32_t y = input_hpixnum/hnx;
  double lam = 90.0-RadToDeg*acos(1.0-2.0*(y+y0+0.5)/ny);
  double deta = 2.5/(input_resolution/4);

  double eta_step = theta;
  eta_step *=  1.0 +
    lam*lam*(0.000192312 - lam*lam*(1.82764e-08 - 1.28162e-11*lam*lam));
  uint8_t etastep = 1;

  while (eta_step > etastep*deta) etastep++;

  return etastep;
}

} // end namespace Stomp
