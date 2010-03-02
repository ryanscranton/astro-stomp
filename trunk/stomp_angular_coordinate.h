// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the classes used to store point-like data on the
// sphere.  This can be simply a location on the sky (the AngularCoordinate
// class) or a location along with additional information about the object at
// that position (WeightedAngularCoordinate).


#ifndef STOMP_ANGULAR_COORDINATE_H
#define STOMP_ANGULAR_COORDINATE_H

#include <math.h>
#include <string>
#include <vector>
#include <map>

namespace Stomp {

class Pixel;  // class declaration in stomp_pixel.h
class AngularCoordinate;
class WeightedAngularCoordinate;

typedef std::vector<AngularCoordinate> AngularVector;
typedef AngularVector::iterator AngularIterator;
typedef std::vector<AngularCoordinate *> AngularPtrVector;
typedef AngularPtrVector::iterator AngularPtrIterator;

typedef std::map<std::string, double> FieldDict;
typedef FieldDict::iterator FieldIterator;

typedef std::vector<WeightedAngularCoordinate> WAngularVector;
typedef WAngularVector::iterator WAngularIterator;
typedef std::vector<WeightedAngularCoordinate *> WAngularPtrVector;
typedef WAngularPtrVector::iterator WAngularPtrIterator;

class AngularCoordinate {
  // Our generic class for handling angular positions.  The idea is that
  // locations on the celestial sphere should be abstract objects from which
  // you can draw whatever angular coordinate pair is necessary for a given
  // use case.  AngularCoordinate's can be instantiated with a particular
  // coordinate system in mind or that can be set later on.
 public:
  enum Sphere {
    Survey,
    Equatorial,
    Galactic
  };
  AngularCoordinate(double theta = 0.0, double phi = 0.0,
                    Sphere sphere = Survey);
  AngularCoordinate(double unit_sphere_x, double unit_sphere_y,
                    double unit_sphere_z);
  ~AngularCoordinate();

  // In addition to the angular coordinate, this class also allows you to
  // extract the X-Y-Z Cartesian coordinates of the angular position on a unit
  // sphere.  This method initializes that functionality, but probably
  // shouldn't ever need to be called explicitly since it's called whenever
  // necessary by the associated methods.
  void InitializeUnitSphere();

  // These three methods let you explicitly set the angular coordinate in one
  // of the supported angular coordinate systems.  Calling these methods resets
  // any other previous values that the AngularCoordinate instance may have
  // had.
  void SetSurveyCoordinates(double lambda, double eta);
  void SetEquatorialCoordinates(double ra, double dec);
  void SetGalacticCoordinates(double gal_lon, double gal_lat);
  void SetUnitSphereCoordinates(double unit_sphere_x, double unit_sphere_y,
				double unit_sphere_z);

  // The basic methods for extracting each of the angular coordinate values.
  // Survey coordinates.
  double Lambda();
  double Eta();

  // Equatorial coordinates
  double RA();
  double DEC();

  // Galactic coordinates
  double GalLon();
  double GalLat();

  // And the associated methods for doing the same with the unit sphere
  // Cartesian coordinates.
  double UnitSphereX();
  double UnitSphereY();
  double UnitSphereZ();

  // Once you have a single angular position on the sphere, very quickly you
  // end up wanting to know the angular distance between your coordinate and
  // another.  This method returns that value in degrees.
  double AngularDistance(AngularCoordinate& ang);
  double AngularDistance(AngularCoordinate* ang);

  // And these two methods return the dot-product and cross-product between
  // the unit vector represented by your angular position and another.
  double DotProduct(AngularCoordinate& ang);
  double DotProduct(AngularCoordinate* ang);
  AngularCoordinate CrossProduct(AngularCoordinate& ang);
  AngularCoordinate CrossProduct(AngularCoordinate* ang);

  // Given another AngularCoordinate, return the AngularCoordinate that defines
  // the great circle running through the current point and the input point.
  // This is different than just the cross-product in that we need to do the
  // calculation in the native coordinate system of the current point.
  void GreatCircle(AngularCoordinate& ang, AngularCoordinate& great_circle);

  // Given either another AngularCoordinate or a Pixel, return the position
  // angle (in degrees, East of North) from the current coordinate to the input
  // location on the sphere.  This angle will, of course, depend on the
  // projection of the current AngularCoordinate.  We also provide variants
  // where the return values are sines and cosines of the position angle.
  double PositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double PositionAngle(Pixel& pix, Sphere sphere = Equatorial);

  double CosPositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double CosPositionAngle(Pixel& pix, Sphere sphere = Equatorial);

  double SinPositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double SinPositionAngle(Pixel& pix, Sphere sphere = Equatorial);

  // Given an input position on the sphere, we can rotate our current position
  // about that point.  The rotation angle will be in degrees, East of
  // North.  In the second variant, we leave our current position static and
  // return a copy of the rotated position.
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle);
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle,
	      AngularCoordinate& rotated_ang);

  // Static methods for when users want to switch between coordinate systems
  // without instantiating the class.
  static void SurveyToGalactic(double lambda, double eta,
                               double& gal_lon, double& gal_lat);
  static void SurveyToEquatorial(double lambda, double eta,
                                 double& ra, double& dec);
  static void EquatorialToSurvey(double ra, double dec,
                                 double& lambda, double& eta);
  static void EquatorialToGalactic(double ra, double dec,
                                   double& gal_lon, double& gal_lat);
  static void GalacticToSurvey(double gal_lon, double gal_lat,
                               double& lambda, double& eta);
  static void GalacticToEquatorial(double gal_lon, double gal_lat,
                                   double& ra, double& dec);
  static void SurveyToXYZ(double lambda, double eta,
			  double& x, double& y, double& z);
  static void EquatorialToXYZ(double ra, double dec,
			      double& x, double& y, double& z);
  static void GalacticToXYZ(double gal_lon, double gal_lat,
			    double& x, double& y, double& z);

  // This is a bit more obscure.  The idea here is that, when you want to find
  // the pixel bounds that subtend a given angular scale about a point on the
  // sphere, finding those bounds in latitude is easier than in longitude.
  // Given a latitude in one of the coordinate systems, the multiplier tells
  // you how many more pixels you should check in the longitude direction
  // relative to the latitude direction.
  static double EtaMultiplier(double lam);
  static double RAMultiplier(double dec);
  static double GalLonMultiplier(double glat);

 private:
  double us_x_, us_y_, us_z_;
};

class WeightedAngularCoordinate : public AngularCoordinate {
  // Sub-class of AngularCoordinate where we attach a weight value to that
  // angular position.

 public:
  WeightedAngularCoordinate();
  WeightedAngularCoordinate(double theta, double phi,
			    double weight, Sphere sphere = Survey);
  WeightedAngularCoordinate(double unit_sphere_x, double unit_sphere_y,
			    double unit_sphere_z, double weight);
  ~WeightedAngularCoordinate();

  // There are two different ways of associating a weight value with the
  // location information we inherit from AngularCoordinate. The first is via
  // the Weight value, which is a simple double precision float.  The second is
  // with the Field values, where we map a string ("magnitude", "ellipticity",
  // etc.) with a double value.  Field allows us to attach a much larger number
  // of values to this location, but does so at the expense of being a slower
  // interface than Weight.
  void SetWeight(double weight);
  double Weight();

  // Associate a Field name with an input value.
  void SetField(const std::string& field_name, double weight);

  // Return the value of the associated Field.  Returns 0.0 if the input Field
  // does not exist.
  double Field(const std::string& field_name);

  // Return some basic facts about the Fields associated with this object.
  uint16_t NFields();
  bool HasFields();

  // Return a vector of Field names.
  void FieldNames(std::vector<std::string>& field_names);

  // Iterators to access the Field values for this object.
  FieldIterator FieldBegin();
  FieldIterator FieldEnd();

  // If we're making a copy of the object, we need to copy over the Field values
  // as well.
  void CopyFields(WeightedAngularCoordinate& w_ang);
  void CopyFields(WeightedAngularCoordinate* w_ang);

  // Accessing the Field values can be considerably slower than the plain
  // Weight value, depending on how many Fields are associated with this
  // object.  If this value is going to be referenced many, many times, it can
  // be useful to temporarily move that Field to the Weight value.
  void CopyFieldToWeight(const std::string& field_name);

  // Internal method to temporarily store the current weight value in a field.
  void _BackUpWeight();

  // Provided that we can find the back-up copy of the original weight,
  // copy that value back into the Weight variable.
  void RestoreOriginalWeight();

 private:
  double weight_;
  FieldDict field_;
};

} // end namespace Stomp

#endif
