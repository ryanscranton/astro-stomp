// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)
// vim: set et ts=2 sw=2:

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
				     Sphere sphere, bool radians) {
  if (!radians) {
    theta *= DegToRad;
    phi *= DegToRad;
  }

  double half_pi = 0.5*Pi;
  double two_pi = 2.0*Pi;

  switch (sphere) {
  case Survey:
    if (DoubleGE(theta, half_pi)) theta = half_pi;
    if (DoubleLE(theta, -half_pi)) theta = -half_pi;
    if (phi > Pi) phi -= two_pi;
    if (phi < -Pi) phi += two_pi;

    us_x_ = -1.0*sin(theta);
    us_y_ = cos(theta)*cos(phi+EtaPole);
    us_z_ = cos(theta)*sin(phi+EtaPole);
    break;
  case Equatorial:
    if (DoubleGE(phi, half_pi)) phi = half_pi;
    if (DoubleLE(phi, -half_pi)) phi = -half_pi;
    if (theta > two_pi) theta -= two_pi;
    if (theta < 0.0) theta += two_pi;

    us_x_ = cos(theta-Node)*cos(phi);
    us_y_ = sin(theta-Node)*cos(phi);
    us_z_ = sin(phi);
    break;
  case Galactic:
    if (DoubleGE(phi, half_pi)) phi = half_pi;
    if (DoubleLE(phi, -half_pi)) phi = -half_pi;
    if (theta > two_pi) theta -= two_pi;
    if (theta < 0.0) theta += two_pi;

    double ra, dec;
    GalacticToEquatorial(theta, phi, ra, dec, true);  // input/output in radians

    us_x_ = cos(ra-Node)*cos(dec);
    us_y_ = sin(ra-Node)*cos(dec);
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

#ifdef WITH_NUMPY
Stomp::AngularCoordinate::Sphere AngularCoordinate::SystemFromString(
  const std::string& system)
  throw (const char* ) {
  Stomp::AngularCoordinate::Sphere sys = Stomp::AngularCoordinate::Survey;
  if (system == "eq" || system == "equatorial") {
    sys = Stomp::AngularCoordinate::Equatorial;
  } else if (system == "sdss" || system == "survey") {
    sys = Stomp::AngularCoordinate::Survey;
  } else if (system == "gal" || system == "galactic") {
    sys = Stomp::AngularCoordinate::Galactic;
  } else {
    std::stringstream err;
    err << "Bad coordinate system indicator '" << system << "'\n\tShould be" <<
      "should be 'survey', 'sdss', 'eq', 'equatorial', 'gal', 'galactic'";
    throw err.str().c_str();
  }
  return sys;
}
#endif

void AngularCoordinate::Set(double theta, double phi, Sphere sphere,
			    bool radians) {
  switch (sphere) {
  case Survey:
    SetSurveyCoordinates(theta, phi, radians);
    break;
  case Equatorial:
    SetEquatorialCoordinates(theta, phi, radians);
    break;
  case Galactic:
    SetGalacticCoordinates(theta, phi, radians);
    break;
  }
}

void AngularCoordinate::SetSurveyCoordinates(double lambda, double eta,
					     bool radians) {
  if (!radians) {
    eta *= DegToRad;
    lambda *= DegToRad;
  }

  double half_pi = 0.5*Pi;
  double two_pi = 2.0*Pi;

  if (DoubleGE(lambda, half_pi)) lambda = half_pi;
  if (DoubleLE(lambda, -half_pi)) lambda = -half_pi;
  if (eta > Pi) eta -= two_pi;
  if (eta < -Pi) eta += two_pi;

  us_x_ = -1.0*sin(lambda);
  us_y_ = cos(lambda)*cos(eta+EtaPole);
  us_z_ = cos(lambda)*sin(eta+EtaPole);
}

void AngularCoordinate::SetEquatorialCoordinates(double ra, double dec,
						 bool radians) {
  if (!radians) {
    ra *= DegToRad;
    dec *= DegToRad;
  }

  double half_pi = 0.5*Pi;
  double two_pi = 2.0*Pi;

  if (DoubleGE(dec, half_pi)) dec = half_pi;
  if (DoubleLE(dec, -half_pi)) dec = -half_pi;
  if (ra > two_pi) ra -= two_pi;
  if (ra < 0.0) ra += two_pi;

  us_x_ = cos(ra-Node)*cos(dec);
  us_y_ = sin(ra-Node)*cos(dec);
  us_z_ = sin(dec);
}

void AngularCoordinate::SetGalacticCoordinates(double gal_lon, double gal_lat,
					       bool radians) {
  if (!radians) {
    gal_lon *= DegToRad;
    gal_lat *= DegToRad;
  }

  double half_pi = 0.5*Pi;
  double two_pi = 2.0*Pi;

  if (DoubleGE(gal_lat, half_pi)) gal_lat = half_pi;
  if (DoubleLE(gal_lat, -half_pi)) gal_lat = -half_pi;
  if (gal_lon > two_pi) gal_lon -= two_pi;
  if (gal_lon < 0.0) gal_lon += two_pi;

  double ra, dec;
  GalacticToEquatorial(gal_lon, gal_lat, ra, dec, true); // I/O in radians

  us_x_ = cos(ra-Node)*cos(dec);
  us_y_ = sin(ra-Node)*cos(dec);
  us_z_ = sin(dec);
}

void AngularCoordinate::SetUnitSphereCoordinates(double unit_sphere_x,
						 double unit_sphere_y,
						 double unit_sphere_z) {
  double r_norm = 1.0/sqrt(unit_sphere_x*unit_sphere_x +
			   unit_sphere_y*unit_sphere_y +
			   unit_sphere_z*unit_sphere_z);

  us_x_ = unit_sphere_x*r_norm;
  us_y_ = unit_sphere_y*r_norm;
  us_z_ = unit_sphere_z*r_norm;
}

void AngularCoordinate::SetUnitSphereCoordinates(double unit_sphere_x,
						 double unit_sphere_y,
						 double unit_sphere_z,
						 Sphere sphere) {
  double r_norm = 1.0/sqrt(unit_sphere_x*unit_sphere_x +
			   unit_sphere_y*unit_sphere_y +
			   unit_sphere_z*unit_sphere_z);
  double phi = asin(unit_sphere_z*r_norm);
  double theta = atan2(unit_sphere_y*r_norm, unit_sphere_x*r_norm);

  switch (sphere) {
  case Survey:
    SetSurveyCoordinates(phi, theta, true);
    break;
  case Equatorial:
    SetEquatorialCoordinates(theta, phi, true);
    break;
  case Galactic:
    SetGalacticCoordinates(theta, phi, true);
    break;
  }
}

double AngularCoordinate::Lambda() {
  return LambdaRadians()*RadToDeg;
}

double AngularCoordinate::Eta() {
  return EtaRadians()*RadToDeg;
}

double AngularCoordinate::LambdaRadians() {
  return -1.0*asin(us_x_);
}

double AngularCoordinate::EtaRadians() {
  double eta = (atan2(us_z_, us_y_) - EtaPole);
  return (DoubleLT(eta, Pi) && DoubleGT(eta, -Pi) ? eta :
	  (DoubleGE(eta, Pi) ? eta - 2.0*Pi : eta + 2.0*Pi)) ;
}

double AngularCoordinate::RA() {
  return RARadians()*RadToDeg;
}

double AngularCoordinate::DEC() {
  return DECRadians()*RadToDeg;
}

double AngularCoordinate::RARadians() {
  double ra = (atan2(us_y_, us_x_) + Node);
  return (DoubleLT(ra, 2.0*Pi) && DoubleGT(ra, 0.0) ? ra :
	  (DoubleGE(ra, 2.0*Pi) ? ra - 2.0*Pi : ra + 2.0*Pi)) ;
}

double AngularCoordinate::DECRadians() {
  return asin(us_z_);
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

double AngularCoordinate::GalLonRadians() {
  double gal_lon, gal_lat;
  EquatorialToGalactic(RARadians(), DECRadians(), gal_lon, gal_lat, true);

  return gal_lon;
}

double AngularCoordinate::GalLatRadians() {
  double gal_lon, gal_lat;
  EquatorialToGalactic(RARadians(), DECRadians(), gal_lon, gal_lat, true);

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

double AngularCoordinate::UnitSphereX(Sphere sphere) {
  double theta = 0.0, phi = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*DegToRad;
    phi = Lambda()*DegToRad;
    break;
  case Equatorial:
    theta = RA()*DegToRad;
    phi = DEC()*DegToRad;
    break;
  case Galactic:
    theta = GalLon()*DegToRad;
    phi = GalLat()*DegToRad;
    break;
  }
  return cos(theta)*cos(phi);
}

double AngularCoordinate::UnitSphereY(Sphere sphere) {
  double theta = 0.0, phi = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*DegToRad;
    phi = Lambda()*DegToRad;
    break;
  case Equatorial:
    theta = RA()*DegToRad;
    phi = DEC()*DegToRad;
    break;
  case Galactic:
    theta = GalLon()*DegToRad;
    phi = GalLat()*DegToRad;
    break;
  }
  return sin(theta)*cos(phi);
}

double AngularCoordinate::UnitSphereZ(Sphere sphere) {
  double phi = 0.0;
  switch (sphere) {
  case Survey:
    phi = Lambda()*DegToRad;
    break;
  case Equatorial:
    phi = DEC()*DegToRad;
    break;
  case Galactic:
    phi = GalLat()*DegToRad;
    break;
  }
  return sin(phi);
}

double AngularCoordinate::AngularDistance(AngularCoordinate& ang) {
  double cos_theta = (us_x_*ang.UnitSphereX() + us_y_*ang.UnitSphereY() +
		      us_z_*ang.UnitSphereZ());
  return RadToDeg*asin(sqrt(fabs(1.0 - cos_theta*cos_theta)));
}

double AngularCoordinate::AngularDistance(AngularCoordinate* ang) {
  double cos_theta = (us_x_*ang->UnitSphereX() + us_y_*ang->UnitSphereY() +
		      us_z_*ang->UnitSphereZ());
  return RadToDeg*asin(sqrt(fabs(1.0 - cos_theta*cos_theta)));
}

double AngularCoordinate::DotProduct(AngularCoordinate& ang) {
  return us_x_*ang.UnitSphereX() + us_y_*ang.UnitSphereY() +
    us_z_*ang.UnitSphereZ();
}

double AngularCoordinate::DotProduct(AngularCoordinate* ang) {
  return us_x_*ang->UnitSphereX() + us_y_*ang->UnitSphereY() +
    us_z_*ang->UnitSphereZ();
}

AngularCoordinate AngularCoordinate::CrossProduct(AngularCoordinate& ang,
						  Sphere sphere) {
  return AngularCoordinate(UnitSphereY(sphere)*ang.UnitSphereZ(sphere) -
			   UnitSphereZ(sphere)*ang.UnitSphereY(sphere),
			   UnitSphereX(sphere)*ang.UnitSphereZ(sphere) -
			   UnitSphereZ(sphere)*ang.UnitSphereX(sphere),
			   UnitSphereX(sphere)*ang.UnitSphereY(sphere) -
			   UnitSphereY(sphere)*ang.UnitSphereX(sphere));
}

AngularCoordinate AngularCoordinate::CrossProduct(AngularCoordinate* ang,
						  Sphere sphere) {
  return AngularCoordinate(UnitSphereY(sphere)*ang->UnitSphereZ(sphere) -
			   UnitSphereZ(sphere)*ang->UnitSphereY(sphere),
			   UnitSphereX(sphere)*ang->UnitSphereZ(sphere) -
			   UnitSphereZ(sphere)*ang->UnitSphereX(sphere),
			   UnitSphereX(sphere)*ang->UnitSphereY(sphere) -
			   UnitSphereY(sphere)*ang->UnitSphereX(sphere));
}

void AngularCoordinate::GreatCircle(AngularCoordinate& ang,
				    AngularCoordinate& great_circle,
				    Sphere sphere) {
  great_circle = CrossProduct(ang, sphere);
}

double AngularCoordinate::PositionAngle(AngularCoordinate& ang, Sphere sphere) {
  double pos_angle = RadToDeg*atan2(SinPositionAngle(ang, sphere),
				    CosPositionAngle(ang, sphere));
  return (DoubleGT(pos_angle, 0.0) ? pos_angle : pos_angle + 360.0);
}

double AngularCoordinate::PositionAngle(Pixel& pix, Sphere sphere) {
  double pos_angle = RadToDeg*atan2(SinPositionAngle(pix, sphere),
				    CosPositionAngle(pix, sphere));
  return (DoubleGT(pos_angle, 0.0) ? pos_angle : pos_angle + 360.0);
}

double AngularCoordinate::CosPositionAngle(AngularCoordinate& ang,
					   Sphere sphere) {
  double theta = 0.0, phi = 0.0;
  double ang_theta = 0.0, ang_phi = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*DegToRad;
    phi = Lambda()*DegToRad;

    ang_theta = ang.Eta()*DegToRad;
    ang_phi = ang.Lambda()*DegToRad;
    break;
  case Equatorial:
    theta = RA()*DegToRad;
    phi = DEC()*DegToRad;

    ang_theta = ang.RA()*DegToRad;
    ang_phi = ang.DEC()*DegToRad;
    break;
  case Galactic:
    theta = GalLon()*DegToRad;
    phi = GalLat()*DegToRad;

    ang_theta = ang.GalLon()*DegToRad;
    ang_phi = ang.GalLat()*DegToRad;
    break;
  }

  return cos(phi)*tan(ang_phi) - sin(phi)*cos(ang_theta - theta);
}

double AngularCoordinate::CosPositionAngle(Pixel& pix, Sphere sphere) {
  double theta = 0.0, phi = 0.0;
  double pix_theta = 0.0, pix_phi = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*DegToRad;
    phi = Lambda()*DegToRad;

    pix_theta = pix.Eta()*DegToRad;
    pix_phi = pix.Lambda()*DegToRad;
    break;
  case Equatorial:
    theta = RA()*DegToRad;
    phi = DEC()*DegToRad;

    pix_theta = pix.RA()*DegToRad;
    pix_phi = pix.DEC()*DegToRad;
    break;
  case Galactic:
    theta = GalLon()*DegToRad;
    phi = GalLat()*DegToRad;

    pix_theta = pix.GalLon()*DegToRad;
    pix_phi = pix.GalLat()*DegToRad;
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
    theta = Eta()*DegToRad;
    ang_theta = ang.Eta()*DegToRad;
    break;
  case Equatorial:
    theta = RA()*DegToRad;
    ang_theta = ang.RA()*DegToRad;
    break;
  case Galactic:
    theta = GalLon()*DegToRad;
    ang_theta = ang.GalLon()*DegToRad;
    break;
  }

  return sin(ang_theta - theta);
}

double AngularCoordinate::SinPositionAngle(Pixel& pix, Sphere sphere) {
  double theta = 0.0;
  double pix_theta = 0.0;
  switch (sphere) {
  case Survey:
    theta = Eta()*DegToRad;
    pix_theta = pix.Eta()*DegToRad;
    break;
  case Equatorial:
    theta = RA()*DegToRad;
    pix_theta = pix.RA()*DegToRad;
    break;
  case Galactic:
    theta = GalLon()*DegToRad;
    pix_theta = pix.GalLon()*DegToRad;
    break;
  }

  return sin(pix_theta - theta);
}

void AngularCoordinate::Rotate(AngularCoordinate& fixed_ang,
			       double rotation_angle, Sphere sphere) {
  double new_x, new_y, new_z;
  Rotate(fixed_ang, rotation_angle, sphere, new_x, new_y, new_z);

  SetUnitSphereCoordinates(new_x, new_y, new_z, sphere);
}

void AngularCoordinate::Rotate(AngularCoordinate& fixed_ang,
			       double rotation_angle,
			       AngularCoordinate& rotated_ang,
			       Sphere sphere) {
  double new_x, new_y, new_z;
  Rotate(fixed_ang, rotation_angle, sphere, new_x, new_y, new_z);

  rotated_ang.SetUnitSphereCoordinates(new_x, new_y, new_z, sphere);
}

void AngularCoordinate::Rotate(AngularCoordinate& fixed_ang,
			       double rotation_angle, Sphere sphere,
			       double& new_unit_sphere_x,
			       double& new_unit_sphere_y,
			       double& new_unit_sphere_z) {
  // Use a quaternion to rotate our current position about the vector
  // represented by the input AngularCoordinate.  Our quarternion axis is
  //
  // q = cos(theta/2) + fixed_ang*sin(theta/2)
  //
  // In order to have our rotations match up with the orientation of the
  // PositionAngle method, we need to reverse the sign of the rotation angle.
  double cos_theta = cos(-0.5*rotation_angle*DegToRad);
  double sin_theta = sin(-0.5*rotation_angle*DegToRad);

  double q_w = cos_theta;
  double q_x = sin_theta*fixed_ang.UnitSphereX(sphere);
  double q_y = sin_theta*fixed_ang.UnitSphereY(sphere);
  double q_z = sin_theta*fixed_ang.UnitSphereZ(sphere);

  double q_norm = 1.0/sqrt(q_w*q_w + q_x*q_x + q_y*q_y + q_z*q_z);

  q_w *= q_norm;
  q_x *= q_norm;
  q_y *= q_norm;
  q_z *= q_norm;

  double unit_sphere_x = UnitSphereX(sphere);
  double unit_sphere_y = UnitSphereY(sphere);
  double unit_sphere_z = UnitSphereZ(sphere);

  new_unit_sphere_x =
    2.0*((-q_y*q_y - q_z*q_z)*unit_sphere_x +
	 (q_x*q_y - q_w*q_z)*unit_sphere_y +
	 (q_x*q_z + q_w*q_y)*unit_sphere_z) + unit_sphere_x;
  new_unit_sphere_y =
    2.0*((q_x*q_y + q_w*q_z)*unit_sphere_x +
	 (-q_x*q_x - q_z*q_z)*unit_sphere_y +
	 (q_y*q_z - q_w*q_x)*unit_sphere_z) + unit_sphere_y;
  new_unit_sphere_z =
    2.0*((q_x*q_z - q_w*q_y)*unit_sphere_x +
	 (q_y*q_z + q_w*q_x)*unit_sphere_y +
	 (-q_x*q_x - q_y*q_y)*unit_sphere_z) + unit_sphere_z;
}

void AngularCoordinate::GalacticToSurvey(double gal_lon, double gal_lat,
                                         double& lambda, double& eta,
					 bool radians) {
  if (!radians) {
    gal_lon *= DegToRad;
    gal_lat *= DegToRad;
  }
  double ra, dec;

  GalacticToEquatorial(gal_lon, gal_lat, ra, dec, true);
  EquatorialToSurvey(ra, dec, lambda, eta, true);

  if (!radians) {
    lambda *= RadToDeg;
    eta *= RadToDeg;
  }
}

void AngularCoordinate::SurveyToGalactic(double lambda, double eta,
                                         double& gal_lon, double& gal_lat,
					 bool radians) {
  if (!radians) {
    lambda *= DegToRad;
    eta *= DegToRad;
  }

  double ra, dec;
  SurveyToEquatorial(lambda, eta, ra, dec, true);
  EquatorialToGalactic(ra, dec, gal_lon, gal_lat, true);

  if (!radians) {
    gal_lon *= RadToDeg;
    gal_lat *= RadToDeg;
  }
}

void AngularCoordinate::SurveyToEquatorial(double lambda, double eta,
                                           double& ra, double& dec,
					   bool radians) {
  if (!radians) {
    lambda *= DegToRad;
    eta *= DegToRad;
  }

  double x = -1.0*sin(lambda);
  double y = cos(lambda)*cos(eta+EtaPole);
  double z = cos(lambda)*sin(eta+EtaPole);

  ra = (atan2(y,x) + Node);
  if (ra < 0.0) ra += 2.0*Pi;
  if (ra > 2.0*Pi) ra -= 2.0*Pi;
  dec = asin(z);

  if (!radians) {
    ra *= RadToDeg;
    dec *= RadToDeg;
  }
}

void AngularCoordinate::SurveyToXYZ(double lambda, double eta,
                                    double& x, double& y, double& z,
				    bool radians) {
  if (!radians) {
    lambda *= DegToRad;
    eta *= DegToRad;
  }

  x = -1.0*sin(lambda);
  y = cos(lambda)*cos(eta+EtaPole);
  z = cos(lambda)*sin(eta+EtaPole);
}

void AngularCoordinate::EquatorialToSurvey(double ra, double dec,
                                           double& lambda, double& eta,
					   bool radians) {
  if (!radians) {
    ra *= DegToRad;
    dec *= DegToRad;
  }

  double x = cos(ra-Node)*cos(dec);
  double y = sin(ra-Node)*cos(dec);
  double z = sin(dec);

  lambda = -1.0*asin(x);
  eta = (atan2(z,y) - EtaPole);
  if (eta < -Pi) eta += 2.0*Pi;
  if (eta > Pi) eta -= 2.0*Pi;

  if (!radians) {
    lambda *= RadToDeg;
    eta *= RadToDeg;
  }
}

void AngularCoordinate::EquatorialToXYZ(double ra, double dec,
                                        double& x, double& y, double& z,
					bool radians) {
  if (!radians) {
    ra *= DegToRad;
    dec *= DegToRad;
  }

  x = cos(ra-Node)*cos(dec);
  y = sin(ra-Node)*cos(dec);
  z = sin(dec);
}

void AngularCoordinate::EquatorialToGalactic(double ra, double dec,
                                             double& gal_lon, double& gal_lat,
					     bool radians) {
  if (!radians) {
    ra *= DegToRad;
    dec *= DegToRad;
  }

  double g_psi = 0.57477043300;
  double stheta = 0.88998808748;
  double ctheta = 0.45598377618;
  double g_phi = 4.9368292465;

  double a = ra - g_phi;
  double b = dec;

  double sb = sin(b);
  double cb = cos(b);
  double cbsa = cb*sin(a);

  b = -1.0*stheta*cbsa + ctheta*sb;
  if (b > 1.0) b = 1.0;

  double bo = asin(b);

  a = atan2(ctheta*cbsa + stheta*sb, cb*cos(a));

  double ao = (a+g_psi+4.0*Pi);

  while (ao > 2.0*Pi) ao -= 2.0*Pi;

  gal_lon = ao;
  if (gal_lon < 0.0) gal_lon += 2.0*Pi;
  if (gal_lon > 2.0*Pi) gal_lon -= 2.0*Pi;

  gal_lat = bo;

  if (!radians) {
    gal_lon *= RadToDeg;
    gal_lat *= RadToDeg;
  }
}

void AngularCoordinate::GalacticToEquatorial(double gal_lon, double gal_lat,
                                             double& ra, double& dec,
					     bool radians) {
  if (!radians) {
    gal_lon *= DegToRad;
    gal_lat *= DegToRad;
  }

  double g_psi = 4.9368292465;
  double stheta = -0.88998808748;
  double ctheta = 0.45598377618;
  double g_phi = 0.57477043300;

  double a = gal_lon - g_phi;
  double b = gal_lat;

  double sb = sin(b);
  double cb = cos(b);
  double cbsa = cb*sin(a);

  b = -1.0*stheta*cbsa + ctheta*sb;
  if (b > 1.0) b = 1.0;

  double bo = asin(b);

  a = atan2(ctheta*cbsa + stheta*sb,cb*cos(a));

  double ao = (a+g_psi+4.0*Pi);
  while (ao > 2.0*Pi) ao -= 2.0*Pi;

  ra = ao;
  dec = bo;

  if (!radians) {
    ra *= RadToDeg;
    dec *= RadToDeg;
  }
}

void AngularCoordinate::GalacticToXYZ(double gal_lon, double gal_lat,
                                      double& x, double& y, double& z,
				      bool radians) {
  if (!radians) {
    gal_lon *= DegToRad;
    gal_lat *= DegToRad;
  }

  double ra, dec;
  GalacticToEquatorial(gal_lon, gal_lat, ra, dec, true);

  x = cos(ra-Node)*cos(dec);
  y = sin(ra-Node)*cos(dec);
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
					Sphere sphere, bool radians) {
  bool io_success = false;

  if (thetaVec.size() == phiVec.size()) {
    if (!ang.empty()) ang.clear();

    ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      ang.push_back(AngularCoordinate(thetaVec[i], phiVec[i], sphere, radians));

    io_success = true;
  }

  return io_success;
}

bool AngularCoordinate::ToAngularVector(const std::string& input_file,
					AngularVector& ang,
					Sphere sphere,
					bool radians,
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
	  Tokenize(line_string, line_elements);
	  n_lines++;

	  if ((line_elements.size() > theta_column) &&
	      (line_elements.size() > phi_column)) {
	    double theta = strtod(line_elements[theta_column].c_str(), NULL);
	    double phi = strtod(line_elements[phi_column].c_str(), NULL);
	    ang.push_back(AngularCoordinate(theta, phi, sphere, radians));
	  }
	}
      }
      input_file_str.close();
      io_success = true;
    } else {
      std::cout << "Stomp::AngularCoordinate::ToAngularVector - " << input_file << " does not exist!\n";
    }
  }

  return io_success;
}

bool AngularCoordinate::FromAngularVector(AngularVector& ang,
					  std::vector<double>& thetaVec,
					  std::vector<double>& phiVec,
					  Sphere sphere, bool radians) {
  if (!thetaVec.empty()) thetaVec.clear();
  if (!phiVec.empty()) phiVec.clear();

  thetaVec.reserve(ang.size());
  phiVec.reserve(ang.size());

  bool io_success = false;

  for (AngularIterator iter=ang.begin();iter!=ang.end();++iter) {
    switch (sphere) {
    case Survey:
      if (radians) {
	thetaVec.push_back(iter->Lambda());
	phiVec.push_back(iter->Eta());
      } else {
	thetaVec.push_back(iter->LambdaRadians());
	phiVec.push_back(iter->EtaRadians());
      }
      break;
    case Equatorial:
      if (radians) {
	thetaVec.push_back(iter->RARadians());
	phiVec.push_back(iter->DECRadians());
      } else {
	thetaVec.push_back(iter->RA());
	phiVec.push_back(iter->DEC());
      }
    break;
    case Galactic:
      if (radians) {
	thetaVec.push_back(iter->GalLonRadians());
	phiVec.push_back(iter->GalLatRadians());
      } else {
	thetaVec.push_back(iter->GalLon());
	phiVec.push_back(iter->GalLat());
      }
      break;
    }
  }

  if (thetaVec.size() == ang.size()) io_success = true;

  return io_success;
}

bool AngularCoordinate::FromAngularVector(AngularVector& ang,
					  const std::string& output_file,
					  Sphere sphere, bool radians) {
  std::ofstream output_file_str(output_file.c_str());

  bool io_success = false;

  if (output_file_str.is_open()) {
    for (AngularIterator iter=ang.begin();iter!=ang.end();++iter) {
      switch (sphere) {
      case Survey:
	if (radians) {
	  output_file_str << iter->LambdaRadians() <<
	    " " << iter->EtaRadians() << "\n";
	} else {
	  output_file_str << iter->Lambda() << " " << iter->Eta() << "\n";
	}
	break;
      case Equatorial:
	if (radians) {
	  output_file_str << iter->RARadians() <<
	    " " << iter->DECRadians() << "\n";
	} else {
	  output_file_str << iter->RA() << " " << iter->DEC() << "\n";
	}
	break;
      case Galactic:
	if (radians) {
	  output_file_str << iter->GalLonRadians() <<
	    " " << iter->GalLatRadians() << "\n";
	} else {
	  output_file_str << iter->GalLon() << " " << iter->GalLat() << "\n";
	}
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
						     Sphere sphere,
						     bool radians) {
  switch (sphere) {
  case Survey:
    SetSurveyCoordinates(theta, phi, radians);
    break;
  case Equatorial:
    SetEquatorialCoordinates(theta, phi, radians);
    break;
  case Galactic:
    SetGalacticCoordinates(theta, phi, radians);
    break;
  }

  weight_ = weight;
}

WeightedAngularCoordinate::WeightedAngularCoordinate(double theta,
						     double phi,
						     double weight,
						     FieldDict& fields,
						     Sphere sphere,
						     bool radians) {
  switch (sphere) {
  case Survey:
    SetSurveyCoordinates(theta, phi, radians);
    break;
  case Equatorial:
    SetEquatorialCoordinates(theta, phi, radians);
    break;
  case Galactic:
    SetGalacticCoordinates(theta, phi, radians);
    break;
  }

  weight_ = weight;

  for (FieldIterator iter=fields.begin();iter!=fields.end();++iter)
    SetField(iter->first, iter->second);
}

WeightedAngularCoordinate::WeightedAngularCoordinate(double unit_sphere_x,
						     double unit_sphere_y,
						     double unit_sphere_z,
						     double weight) {
  SetUnitSphereCoordinates(unit_sphere_x, unit_sphere_y, unit_sphere_z);

  weight_ = weight;
}

WeightedAngularCoordinate::WeightedAngularCoordinate(double unit_sphere_x,
						     double unit_sphere_y,
						     double unit_sphere_z,
						     double weight,
						     FieldDict& fields) {
  SetUnitSphereCoordinates(unit_sphere_x, unit_sphere_y, unit_sphere_z);

  weight_ = weight;

  for (FieldIterator iter=fields.begin();iter!=fields.end();++iter)
    SetField(iter->first, iter->second);
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
						 Sphere sphere,
						 bool radians) {
  bool io_success = false;

  if ((thetaVec.size() == phiVec.size()) &&
      (thetaVec.size() == weightVec.size())) {
    if (!w_ang.empty()) w_ang.clear();

    w_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      w_ang.push_back(WeightedAngularCoordinate(thetaVec[i], phiVec[i],
						weightVec[i], sphere, radians));

    io_success = true;
  }

  return io_success;
}

bool WeightedAngularCoordinate::ToWAngularVector(std::vector<double>& thetaVec,
						 std::vector<double>& phiVec,
						 double weight,
						 WAngularVector& w_ang,
						 Sphere sphere,
						 bool radians) {
  bool io_success = false;

  if (thetaVec.size() == phiVec.size()) {
    if (!w_ang.empty()) w_ang.clear();

    w_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      w_ang.push_back(WeightedAngularCoordinate(thetaVec[i], phiVec[i],
						weight, sphere, radians));

    io_success = true;
  }

  return io_success;
}

bool WeightedAngularCoordinate::ToWAngularVector(const std::string& input_file,
						 WAngularVector& w_ang,
						 Sphere sphere,
						 bool radians,
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
	  Tokenize(line_string, line_elements);
	  n_lines++;

	  if ((line_elements.size() > theta_column) &&
	      (line_elements.size() > phi_column)) {
	    double theta = strtod(line_elements[theta_column].c_str(), NULL);
	    double phi = strtod(line_elements[phi_column].c_str(), NULL);
	    double weight = 1.0;
	    if ((weight_column > -1) &&
		(line_elements.size() > weight_idx))
	      weight = strtod(line_elements[weight_idx].c_str(), NULL);
	    w_ang.push_back(WeightedAngularCoordinate(theta, phi, weight,
						      sphere, radians));
	  }
	}
      }
      input_file_str.close();
      io_success = true;
    } else {
      std::cout << "Stomp::WeightedAngularCoordinate::ToWAngularVector - " << input_file << " does not exist!\n";
    }
  }

  return io_success;
}

bool WeightedAngularCoordinate::FromWAngularVector(WAngularVector& w_ang,
						   std::vector<double>& theta,
						   std::vector<double>& phi,
						   std::vector<double>& weight,
						   Sphere sphere,
						   bool radians) {
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
      if (radians) {
	theta.push_back(iter->LambdaRadians());
	phi.push_back(iter->EtaRadians());
      } else {
	theta.push_back(iter->Lambda());
	phi.push_back(iter->Eta());
      }
      break;
    case Equatorial:
      if (radians) {
	theta.push_back(iter->RARadians());
	phi.push_back(iter->DECRadians());
      } else {
	theta.push_back(iter->RA());
	phi.push_back(iter->DEC());
      }
    break;
    case Galactic:
      if (radians) {
	theta.push_back(iter->GalLonRadians());
	phi.push_back(iter->GalLatRadians());
      } else {
	theta.push_back(iter->GalLon());
	phi.push_back(iter->GalLat());
      }
      break;
    }
  }

  if (theta.size() == w_ang.size()) io_success = true;

  return io_success;
}

bool WeightedAngularCoordinate::FromWAngularVector(WAngularVector& w_ang,
						   const std::string& out_file,
						   Sphere sphere,
						   bool radians) {
  std::ofstream output_file_str(out_file.c_str());

  bool io_success = false;

  if (output_file_str.is_open()) {
    for (WAngularIterator iter=w_ang.begin();iter!=w_ang.end();++iter) {
      switch (sphere) {
      case Survey:
	if (radians) {
	  output_file_str << iter->LambdaRadians() <<
	    " " << iter->EtaRadians() << " ";
	} else {
	  output_file_str << iter->Lambda() << " " << iter->Eta() << " ";
	}
	break;
      case Equatorial:
	if (radians) {
	  output_file_str << iter->RARadians() <<
	    " " << iter->DECRadians() << " ";
	} else {
	  output_file_str << iter->RA() << " " << iter->DEC() << " ";
	}
	break;
      case Galactic:
	if (radians) {
	  output_file_str << iter->GalLonRadians() <<
	    " " << iter->GalLatRadians() << " ";
	} else {
	  output_file_str << iter->GalLon() << " " << iter->GalLat() << " ";
	}
	break;
      }
      output_file_str << iter->Weight() << "\n";
    }
    output_file_str.close();
    io_success = true;
  }

  return io_success;
}

bool WeightedAngularCoordinate::ToWAngularVector(const std::string& input_file,
						 WAngularVector& w_ang,
						 FieldColumnDict& field_columns,
						 Sphere sphere,
						 bool radians,
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
	  Tokenize(line_string, line_elements);
	  n_lines++;

	  if ((line_elements.size() > theta_column) &&
	      (line_elements.size() > phi_column)) {
	    double theta = strtod(line_elements[theta_column].c_str(), NULL);
	    double phi = strtod(line_elements[phi_column].c_str(), NULL);
	    double weight = 1.0;
	    if ((weight_column > -1) &&
		(line_elements.size() > weight_idx)) {
	      weight = strtod(line_elements[weight_idx].c_str(), NULL);
	    }

	    FieldDict fields;
	    for (FieldColumnIterator iter=field_columns.begin();
		 iter!=field_columns.end();++iter) {
	      fields[iter->first] = 0.0;
	      if (line_elements.size() > iter->second) {
		fields[iter->first] =
		  strtod(line_elements[iter->second].c_str(), NULL);
	      }
	    }

	    w_ang.push_back(WeightedAngularCoordinate(theta, phi, weight,
						      fields, sphere, radians));
	  }
	}
      }
      input_file_str.close();
      io_success = true;
    } else {
      std::cout << "Stomp::WeightedAngularCoordinate::ToWAngularVector - " << input_file << " does not exist!\n";
    }
  }

  return io_success;
}

bool WeightedAngularCoordinate::FromWAngularVector(WAngularVector& w_ang,
						   FieldColumnDict& field_names,
						   const std::string& out_file,
						   Sphere sphere,
						   bool radians,
						   uint8_t theta_column,
						   uint8_t phi_column,
						   uint8_t weight_column) {
  std::ofstream output_file_str(out_file.c_str());

  bool io_success = false;

  uint8_t max_column = theta_column;
  if (phi_column > max_column) max_column = phi_column;
  if (weight_column > max_column) max_column = weight_column;

  for (FieldColumnIterator iter=field_names.begin();
       iter!=field_names.end();++iter) {
    if (iter->second > max_column) max_column = iter->second;
  }
  max_column++;

  std::vector<double> column_values;
  column_values.reserve(max_column);

  if (output_file_str.is_open()) {
    for (WAngularIterator iter=w_ang.begin();iter!=w_ang.end();++iter) {
      for (uint8_t i=0;i<max_column;i++) column_values[i] = 0.0;
      switch (sphere) {
      case Survey:
	if (radians) {
	  column_values[theta_column] = iter->LambdaRadians();
	  column_values[phi_column] = iter->EtaRadians();
	} else {
	  column_values[theta_column] = iter->Lambda();
	  column_values[phi_column] = iter->Eta();
	}
	break;
      case Equatorial:
	if (radians) {
	  column_values[theta_column] = iter->RARadians();
	  column_values[phi_column] = iter->DECRadians();
	} else {
	  column_values[theta_column] = iter->RA();
	  column_values[phi_column] = iter->DEC();
	}
	break;
      case Galactic:
	if (radians) {
	  column_values[theta_column] = iter->GalLonRadians();
	  column_values[phi_column] = iter->GalLatRadians();
	} else {
	  column_values[theta_column] = iter->GalLon();
	  column_values[phi_column] = iter->GalLat();
	}
	break;
      }
      column_values[weight_column] = iter->Weight();
      for (FieldColumnIterator field_iter=field_names.begin();
	   field_iter!=field_names.end();++field_iter) {
	column_values[field_iter->second] = iter->Field(field_iter->first);
      }
      output_file_str << column_values[0];
      for (uint8_t i=1;i<max_column;i++)
	output_file_str << " " << column_values[i];
      output_file_str << "\n";
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
				 double weight, Sphere sphere, bool radians) {
  switch (sphere) {
  case Survey:
    SetSurveyCoordinates(theta, phi, radians);
    break;
  case Equatorial:
    SetEquatorialCoordinates(theta, phi, radians);
    break;
  case Galactic:
    SetGalacticCoordinates(theta, phi, radians);
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
				    Sphere sphere, bool radians) {
  bool io_success = false;

  if ((thetaVec.size() == phiVec.size()) &&
      (thetaVec.size() == redshiftVec.size()) &&
      (thetaVec.size() == weightVec.size())) {
    if (!z_ang.empty()) z_ang.clear();

    z_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      z_ang.push_back(CosmoCoordinate(thetaVec[i], phiVec[i], redshiftVec[i],
				      weightVec[i], sphere, radians));

    io_success = true;
  }

  return io_success;
}

bool CosmoCoordinate::ToCosmoVector(std::vector<double>& thetaVec,
				    std::vector<double>& phiVec,
				    std::vector<double>& redshiftVec,
				    double weight,
				    CosmoVector& z_ang,
				    Sphere sphere, bool radians) {
  bool io_success = false;

  if ((thetaVec.size() == phiVec.size()) &&
      (thetaVec.size() == redshiftVec.size())) {
    if (!z_ang.empty()) z_ang.clear();

    z_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      z_ang.push_back(CosmoCoordinate(thetaVec[i], phiVec[i], redshiftVec[i],
				      weight, sphere, radians));

    io_success = true;
  }

  return io_success;
}

bool CosmoCoordinate::ToCosmoVector(const std::string& input_file,
				    CosmoVector& z_ang,
				    Sphere sphere,
				    bool radians,
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
	  Tokenize(line_string, line_elements);
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
					    weight, sphere, radians));
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
				      Sphere sphere, bool radians) {
  if (!theta.empty()) theta.clear();
  if (!phi.empty()) phi.clear();
  if (!redshift.empty()) redshift.clear();
  if (!weight.empty()) weight.clear();

  theta.reserve(z_ang.size());
  phi.reserve(z_ang.size());
  redshift.reserve(z_ang.size());
  weight.reserve(z_ang.size());

  bool io_success = false;

  for (CosmoIterator iter=z_ang.begin();iter!=z_ang.end();++iter) {
    redshift.push_back(iter->Redshift());
    weight.push_back(iter->Weight());
    switch (sphere) {
    case Survey:
      if (radians) {
	theta.push_back(iter->LambdaRadians());
	phi.push_back(iter->EtaRadians());
      } else {
	theta.push_back(iter->Lambda());
	phi.push_back(iter->Eta());
      }
      break;
    case Equatorial:
      if (radians) {
	theta.push_back(iter->RARadians());
	phi.push_back(iter->DECRadians());
      } else {
	theta.push_back(iter->RA());
	phi.push_back(iter->DEC());
      }
    break;
    case Galactic:
      if (radians) {
	theta.push_back(iter->GalLonRadians());
	phi.push_back(iter->GalLatRadians());
      } else {
	theta.push_back(iter->GalLon());
	phi.push_back(iter->GalLat());
      }
      break;
    }
  }

  if (theta.size() == z_ang.size()) io_success = true;

  return io_success;
}

bool CosmoCoordinate::FromCosmoVector(CosmoVector& z_ang,
				      const std::string& output_file,
				      Sphere sphere, bool radians) {
  std::ofstream output_file_str(output_file.c_str());

  bool io_success = false;

  if (output_file_str.is_open()) {
    for (CosmoIterator iter=z_ang.begin();iter!=z_ang.end();++iter) {
      switch (sphere) {
      case Survey:
	if (radians) {
	  output_file_str << iter->LambdaRadians() <<
	    " " << iter->EtaRadians() << " ";
	} else {
	  output_file_str << iter->Lambda() << " " << iter->Eta() << " ";
	}
	break;
      case Equatorial:
	if (radians) {
	  output_file_str << iter->RARadians() <<
	    " " << iter->DECRadians() << " ";
	} else {
	  output_file_str << iter->RA() << " " << iter->DEC() << " ";
	}
	break;
      case Galactic:
	if (radians) {
	  output_file_str << iter->GalLonRadians() <<
	    " " << iter->GalLatRadians() << " ";
	} else {
	  output_file_str << iter->GalLon() << " " << iter->GalLat() << " ";
	}
	break;
      }
      output_file_str << iter->Redshift() << " " << iter->Weight() << "\n";
    }
    output_file_str.close();
    io_success = true;
  }

  return io_success;
}

IndexedAngularCoordinate::IndexedAngularCoordinate() {
  index_ = 0.0;
}

IndexedAngularCoordinate::IndexedAngularCoordinate(double theta,
						   double phi,
						   uint32_t index,
						   Sphere sphere,
						   bool radians) {
  switch (sphere) {
  case Survey:
    SetSurveyCoordinates(theta, phi, radians);
    break;
  case Equatorial:
    SetEquatorialCoordinates(theta, phi, radians);
    break;
  case Galactic:
    SetGalacticCoordinates(theta, phi, radians);
    break;
  }

  index_ = index;
}

IndexedAngularCoordinate::IndexedAngularCoordinate(double unit_sphere_x,
						   double unit_sphere_y,
						   double unit_sphere_z,
						   uint32_t index) {
  SetUnitSphereCoordinates(unit_sphere_x, unit_sphere_y, unit_sphere_z);

  index_ = index;
}

IndexedAngularCoordinate::~IndexedAngularCoordinate() {
  index_ = 0.0;
}

void IndexedAngularCoordinate::SetIndex(uint32_t index) {
  index_ = index;
}

uint32_t IndexedAngularCoordinate::Index() {
  return index_;
}

bool IndexedAngularCoordinate::ToIAngularVector(std::vector<double>& thetaVec,
						std::vector<double>& phiVec,
						std::vector<uint32_t>& indexVec,
						IAngularVector& i_ang,
						Sphere sphere, bool radians) {
  bool io_success = false;

  if ((thetaVec.size() == phiVec.size()) &&
      (thetaVec.size() == indexVec.size())) {
    if (!i_ang.empty()) i_ang.clear();

    i_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      i_ang.push_back(IndexedAngularCoordinate(thetaVec[i], phiVec[i],
					       indexVec[i], sphere, radians));

    io_success = true;
  }

  return io_success;
}

bool IndexedAngularCoordinate::ToIAngularVector(std::vector<double>& thetaVec,
						std::vector<double>& phiVec,
						IAngularVector& i_ang,
						Sphere sphere, bool radians) {
  bool io_success = false;

  if (thetaVec.size() == phiVec.size()) {
    if (!i_ang.empty()) i_ang.clear();

    i_ang.reserve(thetaVec.size());

    for (uint32_t i=0;i<thetaVec.size();i++)
      i_ang.push_back(IndexedAngularCoordinate(thetaVec[i], phiVec[i],
					       i, sphere, radians));

    io_success = true;
  }

  return io_success;
}

bool IndexedAngularCoordinate::ToIAngularVector(const std::string& input_file,
						IAngularVector& i_ang,
						Sphere sphere,
						bool radians,
						uint8_t theta_column,
						uint8_t phi_column,
						int8_t index_column) {

  if (!i_ang.empty()) i_ang.clear();
  bool io_success = false;

  if (theta_column != phi_column) {
    std::ifstream input_file_str(input_file.c_str());

    uint32_t n_lines = 0;
    uint8_t index_idx = static_cast<uint8_t>(index_column);

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
	  Tokenize(line_string, line_elements);

	  if ((line_elements.size() > theta_column) &&
	      (line_elements.size() > phi_column)) {
	    double theta = strtod(line_elements[theta_column].c_str(), NULL);
	    double phi = strtod(line_elements[phi_column].c_str(), NULL);
	    uint32_t index = n_lines;
	    if ((index_column > -1) &&
		(line_elements.size() > index_idx))
	      index =
		static_cast<uint32_t>(strtoul(line_elements[index_idx].c_str(),
					      NULL, 10));
	    i_ang.push_back(IndexedAngularCoordinate(theta, phi, index,
						     sphere, radians));
	  }
	  n_lines++;
	}
      }
      input_file_str.close();
      io_success = true;
    }
  }

  return io_success;
}

bool IndexedAngularCoordinate::FromIAngularVector(IAngularVector& i_ang,
						  std::vector<double>& theta,
						  std::vector<double>& phi,
						  std::vector<uint32_t>& index,
						  Sphere sphere, bool radians) {
  if (!theta.empty()) theta.clear();
  if (!phi.empty()) phi.clear();
  if (!index.empty()) index.clear();

  theta.reserve(i_ang.size());
  phi.reserve(i_ang.size());
  index.reserve(i_ang.size());

  bool io_success = false;

  for (IAngularIterator iter=i_ang.begin();iter!=i_ang.end();++iter) {
    index.push_back(iter->Index());
    switch (sphere) {
    case Survey:
      if (radians) {
	theta.push_back(iter->LambdaRadians());
	phi.push_back(iter->EtaRadians());
      } else {
	theta.push_back(iter->Lambda());
	phi.push_back(iter->Eta());
      }
      break;
    case Equatorial:
      if (radians) {
	theta.push_back(iter->RARadians());
	phi.push_back(iter->DECRadians());
      } else {
	theta.push_back(iter->RA());
	phi.push_back(iter->DEC());
      }
    break;
    case Galactic:
      if (radians) {
	theta.push_back(iter->GalLonRadians());
	phi.push_back(iter->GalLatRadians());
      } else {
	theta.push_back(iter->GalLon());
	phi.push_back(iter->GalLat());
      }
      break;
    }
  }

  if (theta.size() == i_ang.size()) io_success = true;

  return io_success;
}

bool IndexedAngularCoordinate::FromIAngularVector(
  IAngularVector& i_ang, const std::string& output_file, Sphere sphere,
  bool radians) {
  std::ofstream output_file_str(output_file.c_str());

  bool io_success = false;

  if (output_file_str.is_open()) {
    for (IAngularIterator iter=i_ang.begin();iter!=i_ang.end();++iter) {
      switch (sphere) {
      case Survey:
	if (radians) {
	  output_file_str << iter->LambdaRadians() <<
	    " " << iter->EtaRadians() << " ";
	} else {
	  output_file_str << iter->Lambda() << " " << iter->Eta() << " ";
	}
	break;
      case Equatorial:
	if (radians) {
	  output_file_str << iter->RARadians() <<
	    " " << iter->DECRadians() << " ";
	} else {
	  output_file_str << iter->RA() << " " << iter->DEC() << " ";
	}
	break;
      case Galactic:
	if (radians) {
	  output_file_str << iter->GalLonRadians() <<
	    " " << iter->GalLatRadians() << " ";
	} else {
	  output_file_str << iter->GalLon() << " " << iter->GalLat() << " ";
	}
	break;
      }
      output_file_str << iter->Index() << "\n";
    }
    output_file_str.close();
    io_success = true;
  }

  return io_success;
}

} // end namespace Stomp
