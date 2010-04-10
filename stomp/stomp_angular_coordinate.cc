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


// set by Sphere type, more scriptable from python
void AngularCoordinate::Set(
		double theta,
		double phi,
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

bool AngularCoordinate::ToAngularVector(std::vector<double>& thetaVec,
					std::vector<double>& phiVec,
					AngularVector& ang,
					Sphere sphere) {
  bool io_success = false;

  if (thetaVec.size() == phiVec.size()) {
    if (!ang.empty()) ang.clear();

    ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      ang.push_back(AngularCoordinate(thetaVec[i], phiVec[i], sphere));

    io_success = true;
  }

  return io_success;
}

bool AngularCoordinate::ToAngularVector(const std::string& input_file,
					AngularVector& ang,
					Sphere sphere,
					uint8_t theta_column,
					uint8_t phi_column) {
  if (!ang.empty()) ang.clear();
  bool io_success = false;

  if (theta_column != phi_column) {
    std::ifstream input_file_str(input_file.c_str());

    uint32_t n_lines = 0;

    if (input_file_str) {
      while (!input_file_str.eof()) {
	if (!input_file_str.eof()) {
	  // This should read each line into a buffer, convert that buffer into
	  // a string and then break that string into a vector of strings.  We
	  // should then be able to access theta and phi by converting the
	  // appropriate elements of that vector to doubles.
	  char line_buffer[1000];
	  std::vector<std::string> line_elements;

	  input_file_str.getline(line_buffer, 1000);
	  std::string line_string(line_buffer);
	  Stomp::Tokenize(line_string, line_elements);
	  n_lines++;

	  if ((line_elements.size() > theta_column) &&
	      (line_elements.size() > phi_column)) {
	    double theta = strtod(line_elements[theta_column].c_str(), NULL);
	    double phi = strtod(line_elements[phi_column].c_str(), NULL);
	    ang.push_back(AngularCoordinate(theta, phi, sphere));
	  }
	}
      }
      input_file_str.close();
      io_success = true;
    }
  }

  return io_success;
}

bool AngularCoordinate::FromAngularVector(AngularVector& ang,
					  std::vector<double>& thetaVec,
					  std::vector<double>& phiVec,
					  Sphere sphere) {
  if (!thetaVec.empty()) thetaVec.clear();
  if (!phiVec.empty()) phiVec.clear();

  thetaVec.reserve(ang.size());
  phiVec.reserve(ang.size());

  bool io_success = false;

  for (AngularIterator iter=ang.begin();iter!=ang.end();++iter) {
    switch (sphere) {
    case Survey:
      thetaVec.push_back(iter->Lambda());
      phiVec.push_back(iter->Eta());
      break;
    case Equatorial:
      thetaVec.push_back(iter->RA());
      phiVec.push_back(iter->DEC());
    break;
    case Galactic:
      thetaVec.push_back(iter->GalLon());
      phiVec.push_back(iter->GalLat());
      break;
    }
  }

  if (thetaVec.size() == ang.size()) io_success = true;

  return io_success;
}

bool AngularCoordinate::FromAngularVector(AngularVector& ang,
					  const std::string& output_file,
					  Sphere sphere) {
  std::ofstream output_file_str(output_file.c_str());

  bool io_success = false;

  if (output_file_str.is_open()) {
    for (AngularIterator iter=ang.begin();iter!=ang.end();++iter) {
      switch (sphere) {
      case Survey:
	output_file_str << iter->Lambda() << " " << iter->Eta() << "\n";
	break;
      case Equatorial:
	output_file_str << iter->RA() << " " << iter->DEC() << "\n";
	break;
      case Galactic:
	output_file_str << iter->GalLon() << " " << iter->GalLat() << "\n";
	break;
      }
    }
    output_file_str.close();
    io_success = true;
  }

  return io_success;
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

bool WeightedAngularCoordinate::ToWAngularVector(std::vector<double>& thetaVec,
						 std::vector<double>& phiVec,
						 std::vector<double>& weightVec,
						 WAngularVector& w_ang,
						 Sphere sphere) {
  bool io_success = false;

  if ((thetaVec.size() == phiVec.size()) &&
      (thetaVec.size() == weightVec.size())) {
    if (!w_ang.empty()) w_ang.clear();

    w_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      w_ang.push_back(WeightedAngularCoordinate(thetaVec[i], phiVec[i],
						weightVec[i], sphere));

    io_success = true;
  }

  return io_success;
}

bool WeightedAngularCoordinate::ToWAngularVector(std::vector<double>& thetaVec,
						 std::vector<double>& phiVec,
						 double weight,
						 WAngularVector& w_ang,
						 Sphere sphere) {
  bool io_success = false;

  if (thetaVec.size() == phiVec.size()) {
    if (!w_ang.empty()) w_ang.clear();

    w_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      w_ang.push_back(WeightedAngularCoordinate(thetaVec[i], phiVec[i],
						weight, sphere));

    io_success = true;
  }

  return io_success;
}

bool WeightedAngularCoordinate::ToWAngularVector(const std::string& input_file,
						 WAngularVector& w_ang,
						 Sphere sphere,
						 uint8_t theta_column,
						 uint8_t phi_column,
						 int8_t weight_column) {

  if (!w_ang.empty()) w_ang.clear();
  bool io_success = false;

  if (theta_column != phi_column) {
    std::ifstream input_file_str(input_file.c_str());

    uint32_t n_lines = 0;
    uint8_t weight_idx = static_cast<uint8_t>(weight_column);

    if (input_file_str) {
      while (!input_file_str.eof()) {
	if (!input_file_str.eof()) {
	  // This should read each line into a buffer, convert that buffer into
	  // a string and then break that string into a vector of strings.  We
	  // should then be able to access theta and phi by converting the
	  // appropriate elements of that vector to doubles.
	  char line_buffer[1000];
	  std::vector<std::string> line_elements;

	  input_file_str.getline(line_buffer, 1000);
	  std::string line_string(line_buffer);
	  Stomp::Tokenize(line_string, line_elements);
	  n_lines++;

	  if ((line_elements.size() > theta_column) &&
	      (line_elements.size() > phi_column)) {
	    double theta = strtod(line_elements[theta_column].c_str(), NULL);
	    double phi = strtod(line_elements[phi_column].c_str(), NULL);
	    double weight = 1.0;
	    if ((weight_column > -1) &&
		(line_elements.size() > weight_idx))
	      weight = strtod(line_elements[weight_idx].c_str(), NULL);
	    w_ang.push_back(WeightedAngularCoordinate(theta, phi,
						      weight, sphere));
	  }
	}
      }
      input_file_str.close();
      io_success = true;
    }
  }

  return io_success;
}

bool WeightedAngularCoordinate::FromWAngularVector(WAngularVector& w_ang,
						   std::vector<double>& theta,
						   std::vector<double>& phi,
						   std::vector<double>& weight,
						   Sphere sphere) {
  if (!theta.empty()) theta.clear();
  if (!phi.empty()) phi.clear();
  if (!weight.empty()) weight.clear();

  theta.reserve(w_ang.size());
  phi.reserve(w_ang.size());
  weight.reserve(w_ang.size());

  bool io_success = false;

  for (WAngularIterator iter=w_ang.begin();iter!=w_ang.end();++iter) {
    weight.push_back(iter->Weight());
    switch (sphere) {
    case Survey:
      theta.push_back(iter->Lambda());
      phi.push_back(iter->Eta());
      break;
    case Equatorial:
      theta.push_back(iter->RA());
      phi.push_back(iter->DEC());
    break;
    case Galactic:
      theta.push_back(iter->GalLon());
      phi.push_back(iter->GalLat());
      break;
    }
  }

  if (theta.size() == w_ang.size()) io_success = true;

  return io_success;
}

bool WeightedAngularCoordinate::FromWAngularVector(WAngularVector& w_ang,
						   const std::string& out_file,
						   Sphere sphere) {
  std::ofstream output_file_str(out_file.c_str());

  bool io_success = false;

  if (output_file_str.is_open()) {
    for (WAngularIterator iter=w_ang.begin();iter!=w_ang.end();++iter) {
      switch (sphere) {
      case Survey:
	output_file_str << iter->Lambda() << " " << iter->Eta() << " ";
	break;
      case Equatorial:
	output_file_str << iter->RA() << " " << iter->DEC() << " ";
	break;
      case Galactic:
	output_file_str << iter->GalLon() << " " << iter->GalLat() << " ";
	break;
      }
      output_file_str << iter->Weight() << "\n";
    }
    output_file_str.close();
    io_success = true;
  }

  return io_success;
}

bool WeightedAngularCoordinate::AddField(WAngularVector& w_ang,
					 const std::string& field_name,
					 std::vector<double>& field_value) {
  bool add_success = false;

  if (w_ang.size() == field_value.size()) {
    for (uint32_t i=0;i<w_ang.size();i++)
      w_ang[i].SetField(field_name, field_value[i]);
    add_success = true;
  }

  return add_success;
}

CosmoCoordinate::CosmoCoordinate() {
  redshift_ = 0.0;
}

CosmoCoordinate::CosmoCoordinate(double theta, double phi, double redshift,
				 double weight, Sphere sphere) {
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
				 double unit_sphere_z, double redshift,
				 double weight) {
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

bool CosmoCoordinate::ToCosmoVector(std::vector<double>& thetaVec,
				    std::vector<double>& phiVec,
				    std::vector<double>& redshiftVec,
				    std::vector<double>& weightVec,
				    CosmoVector& z_ang,
				    Sphere sphere) {
  bool io_success = false;

  if ((thetaVec.size() == phiVec.size()) &&
      (thetaVec.size() == redshiftVec.size()) &&
      (thetaVec.size() == weightVec.size())) {
    if (!z_ang.empty()) z_ang.clear();

    z_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      z_ang.push_back(CosmoCoordinate(thetaVec[i], phiVec[i], redshiftVec[i],
				      weightVec[i], sphere));

    io_success = true;
  }

  return io_success;
}

bool CosmoCoordinate::ToCosmoVector(std::vector<double>& thetaVec,
				    std::vector<double>& phiVec,
				    std::vector<double>& redshiftVec,
				    double weight,
				    CosmoVector& z_ang,
				    Sphere sphere) {
  bool io_success = false;

  if ((thetaVec.size() == phiVec.size()) &&
      (thetaVec.size() == redshiftVec.size())) {
    if (!z_ang.empty()) z_ang.clear();

    z_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      z_ang.push_back(CosmoCoordinate(thetaVec[i], phiVec[i], redshiftVec[i],
				      weight, sphere));

    io_success = true;
  }

  return io_success;
}

bool CosmoCoordinate::ToCosmoVector(const std::string& input_file,
				    CosmoVector& z_ang,
				    Sphere sphere,
				    uint8_t theta_column,
				    uint8_t phi_column,
				    uint8_t redshift_column,
				    int8_t weight_column) {

  if (!z_ang.empty()) z_ang.clear();
  bool io_success = false;

  if (theta_column != phi_column) {
    std::ifstream input_file_str(input_file.c_str());

    uint32_t n_lines = 0;
    uint8_t weight_idx = static_cast<uint8_t>(weight_column);

    if (input_file_str) {
      while (!input_file_str.eof()) {
	if (!input_file_str.eof()) {
	  // This should read each line into a buffer, convert that buffer into
	  // a string and then break that string into a vector of strings.  We
	  // should then be able to access theta and phi by converting the
	  // appropriate elements of that vector to doubles.
	  char line_buffer[1000];
	  std::vector<std::string> line_elements;

	  input_file_str.getline(line_buffer, 1000);
	  std::string line_string(line_buffer);
	  Stomp::Tokenize(line_string, line_elements);
	  n_lines++;

	  if ((line_elements.size() > theta_column) &&
	      (line_elements.size() > phi_column) &&
	      (line_elements.size() > redshift_column)) {
	    double theta = strtod(line_elements[theta_column].c_str(), NULL);
	    double phi = strtod(line_elements[phi_column].c_str(), NULL);
	    double redshift =
	      strtod(line_elements[redshift_column].c_str(), NULL);
	    double weight = 1.0;
	    if ((weight_column > -1) &&
		(line_elements.size() > weight_idx))
	      weight = strtod(line_elements[weight_idx].c_str(), NULL);
	    z_ang.push_back(CosmoCoordinate(theta, phi, redshift,
					    weight, sphere));
	  }
	}
      }
      input_file_str.close();
      io_success = true;
    }
  }

  return io_success;
}

bool CosmoCoordinate::FromCosmoVector(CosmoVector& z_ang,
				      std::vector<double>& theta,
				      std::vector<double>& phi,
				      std::vector<double>& redshift,
				      std::vector<double>& weight,
				      Sphere sphere) {
  if (!theta.empty()) theta.clear();
  if (!phi.empty()) phi.clear();
  if (!redshift.empty()) redshift.clear();
  if (!weight.empty()) weight.clear();

  theta.reserve(z_ang.size());
  phi.reserve(z_ang.size());
  weight.reserve(z_ang.size());

  bool io_success = false;

  for (CosmoIterator iter=z_ang.begin();iter!=z_ang.end();++iter) {
    redshift.push_back(iter->Redshift());
    weight.push_back(iter->Weight());
    switch (sphere) {
    case Survey:
      theta.push_back(iter->Lambda());
      phi.push_back(iter->Eta());
      break;
    case Equatorial:
      theta.push_back(iter->RA());
      phi.push_back(iter->DEC());
    break;
    case Galactic:
      theta.push_back(iter->GalLon());
      phi.push_back(iter->GalLat());
      break;
    }
  }

  if (theta.size() == z_ang.size()) io_success = true;

  return io_success;
}

bool CosmoCoordinate::FromCosmoVector(CosmoVector& z_ang,
				      const std::string& output_file,
				      Sphere sphere) {
  std::ofstream output_file_str(output_file.c_str());

  bool io_success = false;

  if (output_file_str.is_open()) {
    for (CosmoIterator iter=z_ang.begin();iter!=z_ang.end();++iter) {
      switch (sphere) {
      case Survey:
	output_file_str << iter->Lambda() << " " << iter->Eta() << " ";
	break;
      case Equatorial:
	output_file_str << iter->RA() << " " << iter->DEC() << " ";
	break;
      case Galactic:
	output_file_str << iter->GalLon() << " " << iter->GalLat() << " ";
	break;
      }
      output_file_str << iter->Redshift() << " " << iter->Weight() << "\n";
    }
    output_file_str.close();
    io_success = true;
  }

  return io_success;
}

} // end namespace Stomp
