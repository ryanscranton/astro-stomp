// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//

#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"
#include "stomp_map.h"
#include "stomp_footprint.h"

namespace Stomp {

FootprintBound::FootprintBound() {
  area_ = 0.0;
  pixel_area_ = 0.0;
  lammin_ = etamin_ = 200.0;
  lammax_ = etamax_ = -200.0;
  x_min_ = x_max_ = y_min_ = y_max_ = 0;
  max_resolution_level_ = Stomp::MaxPixelLevel;
  found_starting_resolution_ = false;
  found_xy_bounds_ = false;
}

FootprintBound::~FootprintBound() {
  if (!pix_.empty()) pix_.clear();

  area_ = 0.0;
  pixel_area_ = 0.0;
  lammin_ = etamin_ = 200.0;
  lammax_ = etamax_ = -200.0;
  x_min_ = x_max_ = y_min_ = y_max_ = 0;
  max_resolution_level_ = 0;
  found_starting_resolution_ = false;
  found_xy_bounds_ = false;
}

bool FootprintBound::CheckPoint(AngularCoordinate& ang) {
  return true;
}

bool FootprintBound::FindAngularBounds() {
  lammin_ = -90.0;
  lammax_ = 90.0;
  etamin_ = -180.0;
  etamax_ = 180.0;

  return true;
}

bool FootprintBound::FindArea() {
  return HPixArea*MaxSuperpixnum;
}

uint8_t FootprintBound::FindStartingResolutionLevel() {
  uint16_t max_resolution = static_cast<uint16_t>(1 << max_resolution_level_);

  double min_area = Pixel::PixelArea(max_resolution);

  if (area_ < 10.0*min_area) {
    return -1;
  }

  uint16_t starting_resolution = HPixResolution;

  // We want to start things off with the coarsest possible resolution to
  // save time, but we have to be careful that we're not so coarse that we
  // miss parts of the footprint.  This finds the resolution that has pixels
  // about 1/100th the area of the footprint.
  while (area_/Pixel::PixelArea(starting_resolution) <= 100.0)
    starting_resolution *= 2;

  // We've checked against the possibility that our maximum resolution is too
  // coarse to resolve the area above.  If it happens that our starting
  // resolution is higher than our maximum resolution (i.e. the area is
  // between 10 and 100 times the area of our maximum resolution pixel), we
  // reset the starting resolution to the maximum resolution.
  if (starting_resolution > max_resolution)
    starting_resolution = max_resolution;

  if (starting_resolution % 2 == 0) found_starting_resolution_ = true;

  return MostSignificantBit(starting_resolution);
}

bool FootprintBound::FindXYBounds(const uint8_t resolution_level) {
  uint16_t resolution = static_cast<uint16_t>(1 << resolution_level);
  uint32_t nx = Nx0*resolution, ny = Ny0*resolution;
  Pixel::AreaIndex(resolution,lammin_,lammax_,etamin_,etamax_,
		   x_min_, x_max_, y_min_, y_max_);

  // Checking top border
  bool found_pixel = true;
  bool boundary_failure = false;

  Pixel tmp_pix;
  tmp_pix.SetResolution(resolution);

  uint8_t n_iter = 0;
  uint8_t max_iter = 20;
  while (found_pixel && n_iter < max_iter) {
    found_pixel = false;
    uint32_t y = y_max_, nx_pix;

    if ((x_max_ < x_min_) && (x_min_ > nx/2)) {
      nx_pix = nx - x_min_ + x_max_ + 1;
    } else {
      nx_pix = x_max_ - x_min_ + 1;
    }

    for (uint32_t m=0,x=x_min_;m<nx_pix;m++,x++) {
      if (x == nx) x = 0;

      tmp_pix.SetPixnumFromXY(x,y);

      // This if statement checks positions within the pixel against the
      // footprint bound.
      if (ScorePixel(tmp_pix) < -0.000001) {
        found_pixel = true;
        m = nx_pix + 1;
      }
    }

    if (found_pixel) {
      // The exception to that case is if we've already reached the maximum
      // y index for the pixels.  In that case, we're just done.
      if (y_max_ < ny - 1) {
        y_max_++;
      } else {
        found_pixel = false;
      }
    }
    n_iter++;
  }
  if (n_iter == max_iter) boundary_failure = true;

  // Checking bottom border
  found_pixel = true;
  n_iter = 0;
  while (!boundary_failure && found_pixel && n_iter < max_iter) {
    found_pixel = false;
    uint32_t y = y_min_, nx_pix;

    if ((x_max_ < x_min_) && (x_min_ > nx/2)) {
      nx_pix = nx - x_min_ + x_max_ + 1;
    } else {
      nx_pix = x_max_ - x_min_ + 1;
    }

    for (uint32_t m=0,x=x_min_;m<nx_pix;m++,x++) {
      if (x == nx) x = 0;

      tmp_pix.SetPixnumFromXY(x,y);

      // This if statement checks positions within the pixel against the
      // footprint bound.
      if (ScorePixel(tmp_pix) < -0.000001) {
        found_pixel = true;
        m = nx_pix + 1;
      }
    }

    if (found_pixel) {
      // The exception to that case is if we've already reached the minimum
      // y index for the pixels.  In that case, we're just done.
      if (y_min_ > 0) {
        y_min_--;
      } else {
        found_pixel = false;
      }
    }
    n_iter++;
  }

  if (n_iter == max_iter) boundary_failure = true;

  // Checking left border
  found_pixel = true;
  n_iter = 0;
  while (!boundary_failure && found_pixel && n_iter < max_iter) {
    found_pixel = false;
    uint32_t x = x_min_;

    for (uint32_t y=y_min_;y<=y_max_;y++) {

      tmp_pix.SetPixnumFromXY(x,y);

      // This if statement checks positions within the pixel against the
      // footprint bound.
      if (ScorePixel(tmp_pix) < -0.000001) {
        found_pixel = true;
        y = y_max_ + 1;
      }
    }

    if (found_pixel) {
      if (x_min_ == 0) {
        x_min_ = nx - 1;
      } else {
        x_min_--;
      }
    }
    n_iter++;
  }
  if (n_iter == max_iter) boundary_failure = true;

  // Checking right border
  found_pixel = true;
  n_iter = 0;
  while (!boundary_failure && found_pixel && n_iter < max_iter) {
    found_pixel = false;
    uint32_t x = x_max_;

    for (uint32_t y=y_min_;y<=y_max_;y++) {

      tmp_pix.SetPixnumFromXY(x,y);

      // This if statement checks positions within the pixel against the
      // footprint bound.
      if (ScorePixel(tmp_pix) < -0.000001) {
        found_pixel = true;
        y = y_max_ + 1;
      }
    }

    if (found_pixel) {
      if (x_max_ == nx - 1) {
        x_max_ = 0;
      } else {
        x_max_++;
      }
    }
    n_iter++;
  }

  if (n_iter == max_iter) boundary_failure = true;

  return !boundary_failure;
}

bool FootprintBound::Pixelize() {

  if (!pix_.empty()) pix_.clear();
  pixel_area_ = 0.0;

  uint8_t starting_resolution_level = FindStartingResolutionLevel();

  if ((starting_resolution_level < HPixLevel) ||
      (starting_resolution_level > max_resolution_level_)) {
    return false;
  }

  // We need to be careful around the poles since the pixels there get
  // very distorted.
  if ((lammin_ > 85.0) || (lammax_ < -85.0))
    starting_resolution_level = MostSignificantBit(512);

  if (FindXYBounds(starting_resolution_level)) {

    PixelVector resolve_pix, previous_pix;

    for (uint8_t resolution_level=starting_resolution_level;
         resolution_level<=max_resolution_level_;resolution_level++) {
      uint16_t resolution = static_cast<uint16_t>(1 << resolution_level);

      unsigned n_keep = 0;
      uint32_t nx = Nx0*resolution;
      Pixel tmp_pix;
      tmp_pix.SetResolution(resolution);

      double score;
      AngularCoordinate ang;

      if (resolution_level == starting_resolution_level) {
        resolve_pix.clear();
        previous_pix.clear();

        uint32_t nx_pix;
        if ((x_max_ < x_min_) && (x_min_ > nx/2)) {
          nx_pix = nx - x_min_ + x_max_ + 1;
        } else {
          nx_pix = x_max_ - x_min_ + 1;
        }

        for (uint32_t y=y_min_;y<=y_max_;y++) {
          for (uint32_t m=0,x=x_min_;m<nx_pix;m++,x++) {
            if (x==nx) x = 0;
            tmp_pix.SetPixnumFromXY(x,y);

            score = ScorePixel(tmp_pix);

            if (score < -0.99999) {
              tmp_pix.SetWeight(weight_);
              AddToPixelizedArea(resolution);
              pix_.push_back(tmp_pix);
              n_keep++;
            } else {
              if (score < -0.00001) {
                tmp_pix.SetWeight(score);
                resolve_pix.push_back(tmp_pix);
              }
            }
            previous_pix.push_back(tmp_pix);
          }
        }
      } else {
        if (resolve_pix.size() == 0) {
          std::cout << "Missed all pixels in initial search; trying again...\n";
          for (PixelIterator iter=previous_pix.begin();
               iter!=previous_pix.end();++iter) {
            PixelVector sub_pix;
            iter->SubPix(resolution,sub_pix);
            for (PixelIterator sub_iter=sub_pix.begin();
                 sub_iter!=sub_pix.end();++sub_iter)
              resolve_pix.push_back(*sub_iter);
          }
        }

        previous_pix.clear();

	previous_pix.reserve(resolve_pix.size());
        for (PixelIterator iter=resolve_pix.begin();
             iter!=resolve_pix.end();++iter) previous_pix.push_back(*iter);

        resolve_pix.clear();

        uint32_t x_min, x_max, y_min, y_max;

        for (PixelIterator iter=previous_pix.begin();
             iter!=previous_pix.end();++iter) {

          iter->SubPix(resolution,x_min,x_max,y_min,y_max);

          for (uint32_t y=y_min;y<=y_max;y++) {
            for (uint32_t x=x_min;x<=x_max;x++) {
              tmp_pix.SetPixnumFromXY(x,y);

              score = ScorePixel(tmp_pix);

              if (score < -0.99999) {
                tmp_pix.SetWeight(weight_);
                AddToPixelizedArea(resolution);
                pix_.push_back(tmp_pix);
                n_keep++;
              } else {
                if (score < -0.00001) {
                  tmp_pix.SetWeight(score);
                  resolve_pix.push_back(tmp_pix);
                }
              }
            }
          }
        }
      }
    }

    previous_pix.clear();

    uint16_t max_resolution =
      static_cast<uint16_t>(1 << max_resolution_level_);
    if (area_ > pixel_area_) {
      sort(resolve_pix.begin(),resolve_pix.end(),Pixel::WeightedOrder);

      uint32_t n=0;
      double ur_weight = resolve_pix[n].Weight();
      while ((n < resolve_pix.size()) &&
             ((area_ > pixel_area_) ||
              ((resolve_pix[n].Weight() < ur_weight + 0.1) &&
               (resolve_pix[n].Weight() > ur_weight - 0.1)))) {
        ur_weight = resolve_pix[n].Weight();
        resolve_pix[n].SetWeight(weight_);
        AddToPixelizedArea(max_resolution);
        pix_.push_back(resolve_pix[n]);
        n++;
      }
    }

    Pixel::ResolvePixel(pix_);

    return true;
  } else {
    return false;
  }
}

double FootprintBound::ScorePixel(Pixel& pix) {

  double inv_nx = 1.0/static_cast<double>(Nx0*pix.Resolution());
  double inv_ny = 1.0/static_cast<double>(Ny0*pix.Resolution());
  double x = static_cast<double>(pix.PixelX());
  double y = static_cast<double>(pix.PixelY());

  double lammid = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.5)*inv_ny);
  double lammin = 90.0 - RadToDeg*acos(1.0-2.0*(y+1.0)*inv_ny);
  double lammax = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.0)*inv_ny);
  double lam_quart = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.75)*inv_ny);
  double lam_three = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.25)*inv_ny);

  double etamid = RadToDeg*(2.0*Pi*(x+0.5))*inv_nx + EtaOffSet;
  if (DoubleGE(etamid, 180.0)) etamid -= 360.0;
  if (DoubleLE(etamid, -180.0)) etamid += 360.0;

  double etamin = RadToDeg*(2.0*Pi*(x+0.0))*inv_nx + EtaOffSet;
  if (DoubleGE(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLE(etamin, -180.0)) etamin += 360.0;

  double etamax = RadToDeg*(2.0*Pi*(x+1.0))*inv_nx + EtaOffSet;
  if (DoubleGE(etamax, 180.0)) etamax -= 360.0;
  if (DoubleLE(etamax, -180.0)) etamax += 360.0;

  double eta_quart = RadToDeg*(2.0*Pi*(x+0.25))*inv_nx + EtaOffSet;
  if (DoubleGE(eta_quart, 180.0)) eta_quart -= 360.0;
  if (DoubleLE(eta_quart, -180.0)) eta_quart += 360.0;

  double eta_three = RadToDeg*(2.0*Pi*(x+0.75))*inv_nx + EtaOffSet;
  if (DoubleGE(eta_three, 180.0)) eta_three -= 360.0;
  if (DoubleLE(eta_three, -180.0)) eta_three += 360.0;

  double score = 0.0;

  AngularCoordinate ang(lammid, etamid, AngularCoordinate::Survey);
  if (CheckPoint(ang)) score -= 4.0;

  ang.SetSurveyCoordinates(lam_quart,etamid);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,etamid);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lammid,eta_quart);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lammid,eta_quart);
  if (CheckPoint(ang)) score -= 3.0;

  ang.SetSurveyCoordinates(lam_quart,eta_quart);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,eta_quart);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_quart,eta_three);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,eta_three);
  if (CheckPoint(ang)) score -= 3.0;

  ang.SetSurveyCoordinates(lammid,etamax);
  if (CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammid,etamin);
  if (CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammax,etamid);
  if (CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammin,etamid);
  if (CheckPoint(ang)) score -= 2.0;

  ang.SetSurveyCoordinates(lammax,etamax);
  if (CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammax,etamin);
  if (CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammin,etamax);
  if (CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammin,etamin);
  if (CheckPoint(ang)) score -= 1.0;

  return score/40.0;
}

Map::Map* FootprintBound::ExportMap() {
  return new Map::Map(pix_);
}

void FootprintBound::ExportMap(Map& stomp_map) {
  stomp_map.Clear();
  stomp_map.Initialize(pix_);
}

void FootprintBound::SetMaxResolution(uint16_t resolution) {
  max_resolution_level_ = Stomp::MostSignificantBit(resolution);
}

void FootprintBound::SetArea(double input_area) {
  area_ = input_area;
}

double FootprintBound::Area() {
  return area_;
}

void FootprintBound::AddToPixelizedArea(uint16_t resolution) {
  pixel_area_ += Pixel::PixelArea(resolution);
}

double FootprintBound::Weight() {
    return weight_;
}

void FootprintBound::SetWeight(double input_weight) {
  weight_ = input_weight;
}

double FootprintBound::PixelizedArea() {
  return pixel_area_;
}

uint32_t FootprintBound::NPixel() {
  return pix_.size();
}

void FootprintBound::SetAngularBounds(double lammin, double lammax,
				      double etamin, double etamax) {
  lammin_ = lammin;
  lammax_ = lammax;
  etamin_ = etamin;
  etamax_ = etamax;
}

double FootprintBound::LambdaMin() {
  return lammin_;
}

double FootprintBound::LambdaMax() {
  return lammax_;
}

double FootprintBound::EtaMin() {
  return etamin_;
}

double FootprintBound::EtaMax() {
  return etamax_;
}

uint32_t FootprintBound::XMin() {
  return x_min_;
}

uint32_t FootprintBound::XMax() {
  return x_max_;
}

uint32_t FootprintBound::YMin() {
  return y_min_;
}

uint32_t FootprintBound::YMax() {
  return y_max_;
}

PixelIterator FootprintBound::Begin() {
  return pix_.begin();
}

PixelIterator FootprintBound::End() {
  return pix_.end();
}

void FootprintBound::Clear() {
  pix_.clear();
}

CircleBound::CircleBound(const AngularCoordinate& ang,
                         double radius, double weight) {
  SetWeight(weight);

  ang_ = ang;
  radius_ = radius;
  sin2radius_ = sin(radius*DegToRad)*sin(radius*DegToRad);

  FindArea();
  FindAngularBounds();
  SetMaxResolution();
}

CircleBound::~CircleBound() {
  Clear();
  radius_ = sin2radius_ = 0.0;
}

bool CircleBound::FindAngularBounds() {

  double lammin = ang_.Lambda() - radius_;
  if (DoubleLE(lammin, -90.0)) lammin = -90.0;

  double lammax = ang_.Lambda() + radius_;
  if (DoubleGE(lammax, 90.0)) lammax = 90.0;

  // double eta_multiplier =
  // AngularCoordinate::EtaMultiplier(0.5*(lammax+lammin));
  double eta_multiplier = 1.0;

  double etamin = ang_.Eta() - radius_*eta_multiplier;
  if (DoubleGT(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLT(etamin, -180.0)) etamin += 360.0;

  double etamax = ang_.Eta() + radius_*eta_multiplier;
  if (DoubleGT(etamax, 180.0)) etamax -= 360.0;
  if (DoubleLT(etamax, -180.0)) etamax += 360.0;

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool CircleBound::FindArea() {
  SetArea((1.0 -
           cos(radius_*DegToRad))*2.0*Pi*StradToDeg);
  return true;
}

bool CircleBound::CheckPoint(AngularCoordinate& ang) {

  double costheta =
      ang.UnitSphereX()*ang_.UnitSphereX() +
      ang.UnitSphereY()*ang_.UnitSphereY() +
      ang.UnitSphereZ()*ang_.UnitSphereZ();

  if (1.0-costheta*costheta <= sin2radius_ + 1.0e-10) return true;

  return false;
}

WedgeBound::WedgeBound(const AngularCoordinate& ang, double radius,
		       double position_angle_min, double position_angle_max,
		       double weight, AngularCoordinate::Sphere sphere) {
  SetWeight(weight);

  ang_ = ang;
  radius_ = radius;
  sin2radius_ = sin(radius*DegToRad)*sin(radius*DegToRad);
  if (DoubleLT(position_angle_min, position_angle_max)) {
    position_angle_max_ = position_angle_max;
    position_angle_min_ = position_angle_min;
  } else {
    position_angle_max_ = position_angle_min;
    position_angle_min_ = position_angle_max;
  }
  sphere_ = sphere;

  FindArea();
  FindAngularBounds();
  SetMaxResolution();
}

WedgeBound::~WedgeBound() {
  Clear();
  radius_ = sin2radius_ = position_angle_min_ = position_angle_max_ = 0.0;
}

bool WedgeBound::FindAngularBounds() {
  // We should be able to define our bounds based on the three points that
  // define the wedge.
  double lammin = ang_.Lambda();
  double lammax = ang_.Lambda();
  double etamin = ang_.Eta();
  double etamax = ang_.Eta();

  // Now we define a new point directly north of the center of the circle.
  AngularCoordinate start_ang;
  switch(sphere_) {
  case AngularCoordinate::Survey:
    start_ang.SetSurveyCoordinates(ang_.Lambda()+radius_, ang_.Eta());
    break;
  case AngularCoordinate::Equatorial:
    start_ang.SetEquatorialCoordinates(ang_.RA(), ang_.DEC()+radius_);
    break;
  case AngularCoordinate::Galactic:
    start_ang.SetGalacticCoordinates(ang_.GalLon(), ang_.GalLat()+radius_);
    break;
  }

  // Using that point as a reference, we can rotate to our minimum and
  // maximum position angles.
  AngularCoordinate new_ang;
  start_ang.Rotate(ang_, position_angle_min_, new_ang);
  
  if (DoubleLE(new_ang.Lambda(), lammin)) lammin = new_ang.Lambda();
  if (DoubleGE(new_ang.Lambda(), lammax)) lammax = new_ang.Lambda();

  if (DoubleLE(new_ang.Eta(), etamin)) etamin = new_ang.Eta();
  if (DoubleGE(new_ang.Eta(), etamax)) etamax = new_ang.Eta();

  start_ang.Rotate(ang_, position_angle_max_, new_ang);

  if (DoubleLE(new_ang.Lambda(), lammin)) lammin = new_ang.Lambda();
  if (DoubleGE(new_ang.Lambda(), lammax)) lammax = new_ang.Lambda();

  if (DoubleLE(new_ang.Eta(), etamin)) etamin = new_ang.Eta();
  if (DoubleGE(new_ang.Eta(), etamax)) etamax = new_ang.Eta();

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool WedgeBound::FindArea() {
  double circle_area =
    (1.0 - cos(radius_*DegToRad))*2.0*Pi*StradToDeg;
  SetArea(circle_area*(position_angle_max_ - position_angle_min_)/360.0);
  return true;
}

bool WedgeBound::CheckPoint(AngularCoordinate& ang) {

  double costheta =
    ang.UnitSphereX()*ang_.UnitSphereX() +
    ang.UnitSphereY()*ang_.UnitSphereY() +
    ang.UnitSphereZ()*ang_.UnitSphereZ();

  if (1.0-costheta*costheta <= sin2radius_ + 1.0e-10) {
    double position_angle = ang_.PositionAngle(ang);
    if (DoubleGE(position_angle, position_angle_min_) &&
	DoubleLE(position_angle, position_angle_max_))
      return true;
  }
  return false;
}

PolygonBound::PolygonBound(AngularVector& ang, double weight) {

  SetWeight(weight);

  for (AngularIterator iter=ang.begin();iter!=ang.end();++iter)
    ang_.push_back(*iter);

  n_vert_ = ang_.size();

  x_.reserve(n_vert_);
  y_.reserve(n_vert_);
  z_.reserve(n_vert_);
  dot_.reserve(n_vert_);

  for (uint32_t i=0;i<n_vert_;i++) {

    std::vector<double> tmp_x, tmp_y, tmp_z;

    for (uint32_t j=0;j<n_vert_;j++) {
      tmp_x.push_back(ang_[j].UnitSphereX());
      tmp_y.push_back(ang_[j].UnitSphereY());
      tmp_z.push_back(ang_[j].UnitSphereZ());
    }

    for (uint32_t j=0;j<n_vert_;j++) {
      if (j == n_vert_ - 1) {
        x_.push_back(tmp_y[j]*tmp_z[0] - tmp_y[0]*tmp_z[j]);
        y_.push_back(tmp_z[j]*tmp_x[0] - tmp_z[0]*tmp_x[j]);
        z_.push_back(tmp_x[j]*tmp_y[0] - tmp_x[0]*tmp_y[j]);
      } else {
        x_.push_back(tmp_y[j]*tmp_z[j+1] - tmp_y[j+1]*tmp_z[j]);
        y_.push_back(tmp_z[j]*tmp_x[j+1] - tmp_z[j+1]*tmp_x[j]);
        z_.push_back(tmp_x[j]*tmp_y[j+1] - tmp_x[j+1]*tmp_y[j]);
      }

      double amplitude = sqrt(x_[j]*x_[j] + y_[j]*y_[j] + z_[j]*z_[j]);

      x_[j] /= amplitude;
      y_[j] /= amplitude;
      z_[j] /= amplitude;

      dot_.push_back(1.0); // This assumes that we're not at constant DEC.
    }
  }

  FindArea();
  FindAngularBounds();
  SetMaxResolution();
}

PolygonBound::~PolygonBound() {
  ang_.clear();
  x_.clear();
  y_.clear();
  z_.clear();
  dot_.clear();
  n_vert_ = 0;
}

bool PolygonBound::FindAngularBounds() {

  double lammin = 100.0, lammax = -100.0, etamin = 200.0, etamax = -200.0;

  for (uint32_t i=0;i<n_vert_;i++) {
    if (ang_[i].Lambda() < lammin) lammin = ang_[i].Lambda();
    if (ang_[i].Lambda() > lammax) lammax = ang_[i].Lambda();
    if (ang_[i].Eta() < etamin) etamin = ang_[i].Eta();
    if (ang_[i].Eta() > etamax) etamax = ang_[i].Eta();
  }

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool PolygonBound::FindArea() {

  double sum = 0.0;

  for (uint32_t j=0,k=1;j<n_vert_;j++,k++) {
    if (k == n_vert_) k = 0;

    double cm = (-x_[j]*x_[k] - y_[j]*y_[k] - z_[j]*z_[k]);

    sum += acos(cm);
  }

  double tmp_area = (sum - (n_vert_ - 2)*Pi)*StradToDeg;

  if (tmp_area > 4.0*Pi*StradToDeg) {
    std::cout << "Polygon area is over half the sphere.  This is bad.\n";
    return false;
  }

  SetArea(tmp_area);

  return true;
}

bool PolygonBound::CheckPoint(AngularCoordinate& ang) {
  bool in_polygon = true;

  uint32_t n=0;
  while ((n < n_vert_) && (in_polygon)) {

    in_polygon = false;
    double dot = 1.0 - x_[n]*ang.UnitSphereX() -
        y_[n]*ang.UnitSphereY() - z_[n]*ang.UnitSphereZ();
    if (DoubleLE(dot_[n],0.0)) {
      if (DoubleLE(fabs(dot_[n]), dot)) in_polygon = true;
    } else {
      if (DoubleGE(dot_[n], dot)) in_polygon = true;
    }
    n++;
  }

  return in_polygon;
}

} // end namespace Stomp
