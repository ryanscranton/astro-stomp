// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the classes used to store point-like data on the
// sphere.  This can be simply a location on the sky (the AngularCoordinate
// class) or a location along with additional information about the object at
// that position (WeightedAngularCoordinate).

#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"
#include "stomp_util.h"

namespace Stomp {

AngularCoordinate::AngularCoordinate(double theta, double phi,
				     Sphere sphere) {
  switch (sphere) {
  case Survey:
    if (Stomp::DoubleGE(theta, 90.0)) theta = 90.0;
    if (Stomp::DoubleLE(theta, -90.0)) theta = -90.0;
    if (phi > 180.0) phi -= 360.0;
    if (phi < -180.0) phi += 360.0;

    theta *= Stomp::DegToRad;
    phi *= Stomp::DegToRad;

    us_x_ = -1.0*sin(theta);
    us_y_ = cos(theta)*cos(phi+Stomp::EtaPole);
    us_z_ = cos(theta)*sin(phi+Stomp::EtaPole);
    break;
  case Equatorial:
    if (Stomp::DoubleGE(phi, 90.0)) phi = 90.0;
    if (Stomp::DoubleLE(phi, -90.0)) phi = -90.0;
    if (theta > 360.0) theta -= 360.0;
    if (theta < 0.0) theta += 360.0;

    theta *= Stomp::DegToRad;
    phi *= Stomp::DegToRad;

    us_x_ = cos(theta-Stomp::Node)*cos(phi);
    us_y_ = sin(theta-Stomp::Node)*cos(phi);
    us_z_ = sin(phi);
    break;
  case Galactic:
    if (Stomp::DoubleGE(phi, 90.0)) phi = 90.0;
    if (Stomp::DoubleLE(phi, -90.0)) phi = -90.0;
    if (theta > 360.0) theta -= 360.0;
    if (theta < 0.0) theta += 360.0;

    double ra, dec;
    GalacticToEquatorial(theta, phi, ra, dec);

    ra *= Stomp::DegToRad;
    dec *= Stomp::DegToRad;

    us_x_ = cos(ra-Stomp::Node)*cos(dec);
    us_y_ = sin(ra-Stomp::Node)*cos(dec);
    us_z_ = sin(dec);
    break;
  }
}

AngularCoordinate::AngularCoordinate(double unit_sphere_x,
				     double unit_sphere_y,
				     double unit_sphere_z) {
  double r_norm = sqrt(unit_sphere_x*unit_sphere_x +
		       unit_sphere_y*unit_sphere_y +
		       unit_sphere_z*unit_sphere_z);

  us_x_ = unit_sphere_x/r_norm;
  us_y_ = unit_sphere_y/r_norm;
  us_z_ = unit_sphere_z/r_norm;
}

AngularCoordinate::~AngularCoordinate() {
  us_x_ = 0.0;
  us_y_ = 0.0;
  us_z_ = 0.0;
}

void AngularCoordinate::SetSurveyCoordinates(double lambda, double eta) {
  if (Stomp::DoubleGE(lambda, 90.0)) lambda = 90.0;
  if (Stomp::DoubleLE(lambda, -90.0)) lambda = -90.0;
  if (eta > 180.0) eta -= 360.0;
  if (eta < -180.0) eta += 360.0;

  eta *= Stomp::DegToRad;
  lambda *= Stomp::DegToRad;

  us_x_ = -1.0*sin(lambda);
  us_y_ = cos(lambda)*cos(eta+Stomp::EtaPole);
  us_z_ = cos(lambda)*sin(eta+Stomp::EtaPole);
}

void AngularCoordinate::SetEquatorialCoordinates(double ra, double dec) {
  if (Stomp::DoubleGE(dec, 90.0)) dec = 90.0;
  if (Stomp::DoubleLE(dec, -90.0)) dec = -90.0;
  if (ra > 360.0) ra -= 360.0;
  if (ra < 0.0) ra += 360.0;

  ra *= Stomp::DegToRad;
  dec *= Stomp::DegToRad;

  us_x_ = cos(ra-Stomp::Node)*cos(dec);
  us_y_ = sin(ra-Stomp::Node)*cos(dec);
  us_z_ = sin(dec);
}

void AngularCoordinate::SetGalacticCoordinates(double gal_lon, double gal_lat) {
  if (Stomp::DoubleGE(gal_lat, 90.0)) gal_lat = 90.0;
  if (Stomp::DoubleLE(gal_lat, -90.0)) gal_lat = -90.0;
  if (gal_lon > 360.0) gal_lon -= 360.0;
  if (gal_lon < 0.0) gal_lon += 360.0;

  double ra, dec;
  GalacticToEquatorial(gal_lon, gal_lat, ra, dec);

  ra *= Stomp::DegToRad;
  dec *= Stomp::DegToRad;

  us_x_ = cos(ra-Stomp::Node)*cos(dec);
  us_y_ = sin(ra-Stomp::Node)*cos(dec);
  us_z_ = sin(dec);
}

void AngularCoordinate::SetUnitSphereCoordinates(double unit_sphere_x,
						 double unit_sphere_y,
						 double unit_sphere_z) {
  double r_norm = sqrt(unit_sphere_x*unit_sphere_x +
		       unit_sphere_y*unit_sphere_y +
		       unit_sphere_z*unit_sphere_z);

  us_x_ = unit_sphere_x/r_norm;
  us_y_ = unit_sphere_y/r_norm;
  us_z_ = unit_sphere_z/r_norm;
}

double AngularCoordinate::Lambda() {
  return -1.0*asin(us_x_)*Stomp::RadToDeg;
}

double AngularCoordinate::Eta() {
  double eta = (atan2(us_z_, us_y_) - Stomp::EtaPole)*Stomp::RadToDeg;
  return (Stomp::DoubleLT(eta, 180.0) && Stomp::DoubleGT(eta, -180.0) ? eta :
	  (Stomp::DoubleGE(eta, 180) ? eta - 360.0 : eta + 360.0)) ;
}

double AngularCoordinate::RA() {
  double ra = (atan2(us_y_, us_x_) + Stomp::Node)*Stomp::RadToDeg;
  return (Stomp::DoubleLT(ra, 360.0) && Stomp::DoubleGT(ra, 0.0) ? ra :
	  (Stomp::DoubleGE(ra, 360) ? ra - 360.0 : ra + 360.0)) ;
}

double AngularCoordinate::DEC() {
  return asin(us_z_)*Stomp::RadToDeg;
}

double AngularCoordinate::GalLon() {
  double gal_lon, gal_lat;
  EquatorialToGalactic(RA(), DEC(), gal_lon, gal_lat);

  return gal_lon;
}

double AngularCoordinate::GalLat() {
  double gal_lon, gal_lat;
  EquatorialToGalactic(RA(), DEC(), gal_lon, gal_lat);

  return gal_lat;
}

double AngularCoordinate::UnitSphereX() {
  return us_x_;
}

double AngularCoordinate::UnitSphereY() {
  return us_y_;
}

double AngularCoordinate::UnitSphereZ() {
  return us_z_;
}

double AngularCoordinate::AngularDistance(AngularCoordinate& ang) {
  return acos(us_x_*ang.UnitSphereX() + us_y_*ang.UnitSphereY() +
	      us_z_*ang.UnitSphereZ())*Stomp::RadToDeg;
}

double AngularCoordinate::AngularDistance(AngularCoordinate* ang) {
  return acos(us_x_*ang->UnitSphereX() + us_y_*ang->UnitSphereY() +
	      us_z_*ang->UnitSphereZ())*Stomp::RadToDeg;
}

double AngularCoordinate::DotProduct(AngularCoordinate& ang) {
  return us_x_*ang.UnitSphereX() + us_y_*ang.UnitSphereY() +
    us_z_*ang.UnitSphereZ();
}

double AngularCoordinate::DotProduct(AngularCoordinate* ang) {
  return us_x_*ang->UnitSphereX() + us_y_*ang->UnitSphereY() +
    us_z_*ang->UnitSphereZ();
}

AngularCoordinate AngularCoordinate::CrossProduct(AngularCoordinate& ang) {
  return AngularCoordinate(us_y_*ang.UnitSphereZ() - us_z_*ang.UnitSphereY(),
			   us_x_*ang.UnitSphereZ() - us_z_*ang.UnitSphereX(),
			   us_x_*ang.UnitSphereY() - us_y_*ang.UnitSphereX());
}

AngularCoordinate AngularCoordinate::CrossProduct(AngularCoordinate* ang) {
  return AngularCoordinate(us_y_*ang->UnitSphereZ()-us_z_*ang->UnitSphereY(),
			   us_x_*ang->UnitSphereZ()-us_z_*ang->UnitSphereX(),
			   us_x_*ang->UnitSphereY()-us_y_*ang->UnitSphereX());
}

void AngularCoordinate::GreatCircle(AngularCoordinate& ang,
				    AngularCoordinate& great_circle) {
  great_circle = CrossProduct(ang);
}

double AngularCoordinate::PositionAngle(AngularCoordinate& ang, Sphere sphere) {
  return Stomp::RadToDeg*atan2(SinPositionAngle(ang, sphere),
			       CosPositionAngle(ang, sphere));
}

double AngularCoordinate::PositionAngle(Pixel& pix, Sphere sphere) {
  return Stomp::RadToDeg*atan2(SinPositionAngle(pix, sphere),
			       CosPositionAngle(pix, sphere));
}

double AngularCoordinate::CosPositionAngle(AngularCoordinate& ang,
					   Sphere sphere) {
  double theta = 0.0, phi = 0.0;
  double ang_theta = 0.0, ang_phi = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*Stomp::DegToRad;
    phi = Lambda()*Stomp::DegToRad;

    ang_theta = ang.Eta()*Stomp::DegToRad;
    ang_phi = ang.Lambda()*Stomp::DegToRad;
    break;
  case Equatorial:
    theta = RA()*Stomp::DegToRad;
    phi = DEC()*Stomp::DegToRad;

    ang_theta = ang.RA()*Stomp::DegToRad;
    ang_phi = ang.DEC()*Stomp::DegToRad;
    break;
  case Galactic:
    theta = GalLon()*Stomp::DegToRad;
    phi = GalLat()*Stomp::DegToRad;

    ang_theta = ang.GalLon()*Stomp::DegToRad;
    ang_phi = ang.GalLat()*Stomp::DegToRad;
    break;
  }

  return cos(phi)*tan(ang_phi) - sin(phi)*cos(ang_theta - theta);
}

double AngularCoordinate::CosPositionAngle(Pixel& pix, Sphere sphere) {
  double theta = 0.0, phi = 0.0;
  double pix_theta = 0.0, pix_phi = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*Stomp::DegToRad;
    phi = Lambda()*Stomp::DegToRad;

    pix_theta = pix.Eta()*Stomp::DegToRad;
    pix_phi = pix.Lambda()*Stomp::DegToRad;
    break;
  case Equatorial:
    theta = RA()*Stomp::DegToRad;
    phi = DEC()*Stomp::DegToRad;

    pix_theta = pix.RA()*Stomp::DegToRad;
    pix_phi = pix.DEC()*Stomp::DegToRad;
    break;
  case Galactic:
    theta = GalLon()*Stomp::DegToRad;
    phi = GalLat()*Stomp::DegToRad;

    pix_theta = pix.GalLon()*Stomp::DegToRad;
    pix_phi = pix.GalLat()*Stomp::DegToRad;
    break;
  }

  return cos(phi)*tan(pix_phi) - sin(phi)*cos(pix_theta - theta);
}

double AngularCoordinate::SinPositionAngle(AngularCoordinate& ang,
					   Sphere sphere) {
  double theta = 0.0;
  double ang_theta = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*Stomp::DegToRad;
    ang_theta = ang.Eta()*Stomp::DegToRad;
    break;
  case Equatorial:
    theta = RA()*Stomp::DegToRad;
    ang_theta = ang.RA()*Stomp::DegToRad;
    break;
  case Galactic:
    theta = GalLon()*Stomp::DegToRad;
    ang_theta = ang.GalLon()*Stomp::DegToRad;
    break;
  }

  return sin(ang_theta - theta);
}

double AngularCoordinate::SinPositionAngle(Pixel& pix, Sphere sphere) {
  double theta = 0.0;
  double pix_theta = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*Stomp::DegToRad;
    pix_theta = pix.Eta()*Stomp::DegToRad;
    break;
  case Equatorial:
    theta = RA()*Stomp::DegToRad;
    pix_theta = pix.RA()*Stomp::DegToRad;
    break;
  case Galactic:
    theta = GalLon()*Stomp::DegToRad;
    pix_theta = pix.GalLon()*Stomp::DegToRad;
    break;
  }

  return sin(pix_theta - theta);
}

void AngularCoordinate::Rotate(AngularCoordinate& fixed_ang,
			       double rotation_angle) {
  // Using a quaternion to rotate our current position about the vector
  // represented by the input AngularCoordinate.  Our starting quarternion is
  //
  // q = cos(theta/2) + fixed_ang*sin(theta/2)
  double cos_theta = cos(0.5*rotation_angle*Stomp::DegToRad);
  double sin_theta = sin(0.5*rotation_angle*Stomp::DegToRad);

  double q1 = cos_theta;
  double q2 = sin_theta*fixed_ang.UnitSphereX();
  double q3 = sin_theta*fixed_ang.UnitSphereY();
  double q4 = sin_theta*fixed_ang.UnitSphereZ();

  double new_x =
    2.0*((-q3*q3 - q4*q4)*UnitSphereX() +
	 (q2*q3 - q1*q4)*UnitSphereY() +
	 (q1*q3 - q2*q4)*UnitSphereZ()) + UnitSphereX();
  double new_y =
    2.0*((q1*q4 + q2*q3)*UnitSphereX() +
	 (-q2*q2 - q4*q4)*UnitSphereY() +
	 (q3*q4 - q1*q3)*UnitSphereZ()) + UnitSphereY();
  double new_z =
    2.0*((q2*q4 - q1*q3)*UnitSphereX() +
	 (q1*q3 + q3*q4)*UnitSphereY() +
	 (-q2*q2 - q3*q3)*UnitSphereZ()) + UnitSphereZ();

  SetUnitSphereCoordinates(new_x, new_y, new_z);
}

void AngularCoordinate::Rotate(AngularCoordinate& fixed_ang,
			       double rotation_angle,
			       AngularCoordinate& rotated_ang) {
  // Using a quaternion to rotate our current position about the vector
  // represented by the input AngularCoordinate.  Our starting quarternion is
  //
  // q = cos(theta/2) + fixed_ang*sin(theta/2)
  double cos_theta = cos(0.5*rotation_angle*Stomp::DegToRad);
  double sin_theta = sin(0.5*rotation_angle*Stomp::DegToRad);

  double q1 = cos_theta;
  double q2 = sin_theta*fixed_ang.UnitSphereX();
  double q3 = sin_theta*fixed_ang.UnitSphereY();
  double q4 = sin_theta*fixed_ang.UnitSphereZ();

  double new_x =
    2.0*((-q3*q3 - q4*q4)*UnitSphereX() +
	 (q2*q3 - q1*q4)*UnitSphereY() +
	 (q1*q3 - q2*q4)*UnitSphereZ()) + UnitSphereX();
  double new_y =
    2.0*((q1*q4 + q2*q3)*UnitSphereX() +
	 (-q2*q2 - q4*q4)*UnitSphereY() +
	 (q3*q4 - q1*q3)*UnitSphereZ()) + UnitSphereY();
  double new_z =
    2.0*((q2*q4 - q1*q3)*UnitSphereX() +
	 (q1*q3 + q3*q4)*UnitSphereY() +
	 (-q2*q2 - q3*q3)*UnitSphereZ()) + UnitSphereZ();

  rotated_ang.SetUnitSphereCoordinates(new_x, new_y, new_z);
}

void AngularCoordinate::GalacticToSurvey(double gal_lon, double gal_lat,
                                         double& lambda, double& eta) {
  double ra, dec;

  GalacticToEquatorial(gal_lon,gal_lat,ra,dec);
  EquatorialToSurvey(ra,dec,lambda,eta);
}

void AngularCoordinate::SurveyToGalactic(double lambda, double eta,
                                         double& gal_lon, double& gal_lat) {
  double ra, dec;

  SurveyToEquatorial(lambda, eta, ra, dec);
  EquatorialToGalactic(ra,dec, gal_lon, gal_lat);
  if (gal_lon < 0.0) gal_lon += 360.0;
  if (gal_lon > 360.0) gal_lon -= 360.0;
}

void AngularCoordinate::SurveyToEquatorial(double lambda, double eta,
                                           double& ra, double& dec) {
  lambda *= Stomp::DegToRad;
  eta *= Stomp::DegToRad;

  double x = -1.0*sin(lambda);
  double y = cos(lambda)*cos(eta+Stomp::EtaPole);
  double z = cos(lambda)*sin(eta+Stomp::EtaPole);

  ra = (atan2(y,x) + Stomp::Node)*Stomp::RadToDeg;
  if (ra < 0.0) ra += 360.0;
  if (ra > 360.0) ra -= 360.0;
  dec = asin(z)*Stomp::RadToDeg;
}

void AngularCoordinate::SurveyToXYZ(double lambda, double eta,
                                    double& x, double& y, double& z) {
  lambda *= Stomp::DegToRad;
  eta *= Stomp::DegToRad;

  x = -1.0*sin(lambda);
  y = cos(lambda)*cos(eta+Stomp::EtaPole);
  z = cos(lambda)*sin(eta+Stomp::EtaPole);
}

void AngularCoordinate::EquatorialToSurvey(double ra, double dec,
                                           double& lambda, double& eta) {
  ra *= Stomp::DegToRad;
  dec *= Stomp::DegToRad;

  double x = cos(ra-Stomp::Node)*cos(dec);
  double y = sin(ra-Stomp::Node)*cos(dec);
  double z = sin(dec);

  lambda = -1.0*asin(x)*Stomp::RadToDeg;
  eta = (atan2(z,y) - Stomp::EtaPole)*Stomp::RadToDeg;
  if (eta < -180.0) eta += 360.0;
  if (eta > 180.0) eta -= 360.0;
}

void AngularCoordinate::EquatorialToXYZ(double ra, double dec,
                                        double& x, double& y, double& z) {
  ra *= Stomp::DegToRad;
  dec *= Stomp::DegToRad;

  x = cos(ra-Stomp::Node)*cos(dec);
  y = sin(ra-Stomp::Node)*cos(dec);
  z = sin(dec);
}

void AngularCoordinate::EquatorialToGalactic(double ra, double dec,
                                             double& gal_lon, double& gal_lat) {
  double g_psi = 0.57477043300;
  double stheta = 0.88998808748;
  double ctheta = 0.45598377618;
  double g_phi = 4.9368292465;

  double a = ra*Stomp::DegToRad - g_phi;
  double b = dec*Stomp::DegToRad;

  double sb = sin(b);
  double cb = cos(b);
  double cbsa = cb*sin(a);

  b = -1.0*stheta*cbsa + ctheta*sb;
  if (b > 1.0) b = 1.0;

  double bo = asin(b)*Stomp::RadToDeg;

  a = atan2(ctheta*cbsa + stheta*sb,cb*cos(a));

  double ao = (a+g_psi+4.0*Stomp::Pi)*Stomp::RadToDeg;

  while (ao > 360.0) ao -= 360.0;

  gal_lon = ao;
  gal_lat = bo;
}

void AngularCoordinate::GalacticToEquatorial(double gal_lon, double gal_lat,
                                             double& ra, double& dec) {
  double g_psi = 4.9368292465;
  double stheta = -0.88998808748;
  double ctheta = 0.45598377618;
  double g_phi = 0.57477043300;

  double a = gal_lon*Stomp::DegToRad - g_phi;
  double b = gal_lat*Stomp::DegToRad;

  double sb = sin(b);
  double cb = cos(b);
  double cbsa = cb*sin(a);

  b = -1.0*stheta*cbsa + ctheta*sb;
  if (b > 1.0) b = 1.0;

  double bo = asin(b)*Stomp::RadToDeg;

  a = atan2(ctheta*cbsa + stheta*sb,cb*cos(a));

  double ao = (a+g_psi+4.0*Stomp::Pi)*Stomp::RadToDeg;
  while (ao > 360.0) ao -= 360.0;

  ra = ao;
  dec = bo;
}

void AngularCoordinate::GalacticToXYZ(double gal_lon, double gal_lat,
                                      double& x, double& y, double& z) {
  double ra, dec;
  GalacticToEquatorial(gal_lat, gal_lon, ra, dec);
  ra *= Stomp::DegToRad;
  dec *= Stomp::DegToRad;

  x = cos(ra-Stomp::Node)*cos(dec);
  y = sin(ra-Stomp::Node)*cos(dec);
  z = sin(dec);
}

double AngularCoordinate::EtaMultiplier(double lam) {
  return 1.0 +
    lam*lam*(0.000192312 - lam*lam*(1.82764e-08 + 1.28162e-11*lam*lam));
}

double AngularCoordinate::RAMultiplier(double dec) {
  return 1.0 +
    dec*dec*(0.000192312 - dec*dec*(1.82764e-08 + 1.28162e-11*dec*dec));
}

double AngularCoordinate::GalLonMultiplier(double glat) {
  return 1.0 +
    glat*glat*(0.000192312 - glat*glat*(1.82764e-08 + 1.28162e-11*glat*glat));
}

WeightedAngularCoordinate::WeightedAngularCoordinate() {
  weight_ = 0.0;
}

WeightedAngularCoordinate::WeightedAngularCoordinate(double theta,
						     double phi,
						     double weight,
						     Sphere sphere) {
  switch (sphere) {
  case Survey:
    SetSurveyCoordinates(theta, phi);
    break;
  case Equatorial:
    SetEquatorialCoordinates(theta, phi);
    break;
  case Galactic:
    SetGalacticCoordinates(theta, phi);
    break;
  }

  weight_ = weight;
}

WeightedAngularCoordinate::WeightedAngularCoordinate(double unit_sphere_x,
						     double unit_sphere_y,
						     double unit_sphere_z,
						     double weight) {
  SetUnitSphereCoordinates(unit_sphere_x, unit_sphere_y, unit_sphere_z);

  weight_ = weight;
}

WeightedAngularCoordinate::~WeightedAngularCoordinate() {
  weight_ = 0.0;
}

void WeightedAngularCoordinate::SetWeight(double weight) {
  weight_ = weight;
}

double WeightedAngularCoordinate::Weight() {
  return weight_;
}

void WeightedAngularCoordinate::SetField(const std::string& field_name,
					 double weight) {
  field_[field_name] = weight;
}

double WeightedAngularCoordinate::Field(const std::string& field_name) {
  return (field_.find(field_name) != field_.end() ? field_[field_name] : 0.0);
}

uint16_t WeightedAngularCoordinate::NFields() {
  return field_.size();
}

bool WeightedAngularCoordinate::HasFields() {
  return (field_.size() > 0 ? true : false);
}

void WeightedAngularCoordinate::FieldNames(
  std::vector<std::string>& field_names) {
  field_names.clear();
  for (FieldIterator iter=field_.begin();iter!=field_.end();++iter)
    field_names.push_back(iter->first);
}

FieldIterator WeightedAngularCoordinate::FieldBegin() {
  return field_.begin();
}

FieldIterator WeightedAngularCoordinate::FieldEnd() {
  return field_.end();
}

void WeightedAngularCoordinate::CopyFields(WeightedAngularCoordinate& w_ang) {
  for (FieldIterator iter=w_ang.FieldBegin();iter!=w_ang.FieldEnd();++iter)
    field_[iter->first] = iter->second;
}

void WeightedAngularCoordinate::CopyFields(WeightedAngularCoordinate* w_ang) {
  for (FieldIterator iter=w_ang->FieldBegin();iter!=w_ang->FieldEnd();++iter)
    field_[iter->first] = iter->second;
}

void WeightedAngularCoordinate::CopyFieldToWeight(
  const std::string& field_name) {
  if (field_.find(field_name) != field_.end()) {
    _BackUpWeight();
    weight_ = field_[field_name];
  } else {
    weight_ = 0.0;
  }
}

void WeightedAngularCoordinate::_BackUpWeight() {
  std::string temporary_field_name = "temporary_field_name_DO_NOT_USE";
  if (field_.find(temporary_field_name) == field_.end()) {
    field_[temporary_field_name] = weight_;
  }
}

void WeightedAngularCoordinate::RestoreOriginalWeight() {
  std::string temporary_field_name = "temporary_field_name_DO_NOT_USE";
  if (field_.find(temporary_field_name) != field_.end()) {
    weight_ = field_[temporary_field_name];
  }
}

CosmoCoordinate::CosmoCoordinate() {
  redshift_ = 0.0;
}

CosmoCoordinate::CosmoCoordinate(double theta, double phi, double weight,
				 double redshift, Sphere sphere) {
  switch (sphere) {
  case Survey:
    SetSurveyCoordinates(theta, phi);
    break;
  case Equatorial:
    SetEquatorialCoordinates(theta, phi);
    break;
  case Galactic:
    SetGalacticCoordinates(theta, phi);
    break;
  }

  SetWeight(weight);

  redshift_ = redshift;
}

CosmoCoordinate::CosmoCoordinate(double unit_sphere_x, double unit_sphere_y,
				 double unit_sphere_z, double weight,
				 double redshift) {
  SetUnitSphereCoordinates(unit_sphere_x, unit_sphere_y, unit_sphere_z);

  SetWeight(weight);

  redshift_ = redshift;
}

CosmoCoordinate::~CosmoCoordinate() {
  redshift_ = 0.0;
}

double CosmoCoordinate::ProjectedRadius(AngularCoordinate& ang) {
  return Cosmology::ProjectedDistance(redshift_, AngularDistance(ang));
}

double CosmoCoordinate::ProjectedRadius(AngularCoordinate* ang) {
  return Cosmology::ProjectedDistance(redshift_, AngularDistance(ang));
}

double CosmoCoordinate::DotProduct(CosmoCoordinate& ang) {
  double dot_product = (UnitSphereX()*ang.UnitSphereX() +
			UnitSphereY()*ang.UnitSphereY() +
			UnitSphereZ()*ang.UnitSphereZ());
  return ComovingDistance()*ang.ComovingDistance()*dot_product;
}

double CosmoCoordinate::DotProduct(CosmoCoordinate* ang) {
  double dot_product = (UnitSphereX()*ang->UnitSphereX() +
			UnitSphereY()*ang->UnitSphereY() +
			UnitSphereZ()*ang->UnitSphereZ());
  return ComovingDistance()*ang->ComovingDistance()*dot_product;
}

double CosmoCoordinate::ComovingDistance() {
  return Cosmology::ComovingDistance(redshift_);
}

double CosmoCoordinate::AngularDiameterDistance() {
  return Cosmology::AngularDiameterDistance(redshift_);
}

double CosmoCoordinate::LuminosityDistance() {
  return Cosmology::LuminosityDistance(redshift_);
}

double CosmoCoordinate::Redshift() {
  return redshift_;
}

void CosmoCoordinate::SetRedshift(double redshift) {
  redshift_ = redshift;
}

} // end namespace Stomp
