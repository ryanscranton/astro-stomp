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
// that position (WeightedAngularCoordinate).  CosmoCoordinate extends the
// latter functionality to also allow for treatment of cosmological coordinates.


#ifndef STOMP_ANGULAR_COORDINATE_H
#define STOMP_ANGULAR_COORDINATE_H

#include <math.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>

namespace Stomp {

class Pixel;  // class declaration in stomp_pixel.h
class AngularCoordinate;
class WeightedAngularCoordinate;
class CosmoCoordinate;
class IndexedAngularCoordinate;

typedef std::vector<AngularCoordinate> AngularVector;
typedef AngularVector::iterator AngularIterator;
typedef std::vector<AngularCoordinate *> AngularPtrVector;
typedef AngularPtrVector::iterator AngularPtrIterator;

typedef std::map<std::string, double> FieldDict;
typedef FieldDict::iterator FieldIterator;

typedef std::map<std::string, uint8_t> FieldColumnDict;
typedef FieldColumnDict::iterator FieldColumnIterator;

typedef std::vector<WeightedAngularCoordinate> WAngularVector;
typedef WAngularVector::iterator WAngularIterator;
typedef std::vector<WeightedAngularCoordinate *> WAngularPtrVector;
typedef WAngularPtrVector::iterator WAngularPtrIterator;

typedef std::vector<CosmoCoordinate> CosmoVector;
typedef CosmoVector::iterator CosmoIterator;
typedef std::vector<CosmoCoordinate *> CosmoPtrVector;
typedef CosmoPtrVector::iterator CosmoPtrIterator;

typedef std::vector<IndexedAngularCoordinate> IAngularVector;
typedef IAngularVector::iterator IAngularIterator;
typedef std::vector<IndexedAngularCoordinate *> IAngularPtrVector;
typedef IAngularPtrVector::iterator IAngularPtrIterator;

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
                    Sphere sphere = Survey, bool radians = false);
  AngularCoordinate(double unit_sphere_x, double unit_sphere_y,
                    double unit_sphere_z);
  ~AngularCoordinate();

#ifdef WITH_NUMPY
  // Convert a string version of the coordinate system to Sphere; python-only.
  static Sphere SystemFromString(const std::string& system)
    throw (const char* );
#endif

  // These three methods let you explicitly set the angular coordinate in one
  // of the supported angular coordinate systems.  Calling these methods resets
  // any other previous values that the AngularCoordinate instance may have
  // had.
  void SetSurveyCoordinates(double lambda, double eta, bool radians = false);
  void SetEquatorialCoordinates(double ra, double dec, bool radians = false);
  void SetGalacticCoordinates(double gal_lon, double gal_lat,
			      bool radians = false);
  void SetUnitSphereCoordinates(double unit_sphere_x, double unit_sphere_y,
				double unit_sphere_z);
  void SetUnitSphereCoordinates(double unit_sphere_x, double unit_sphere_y,
				double unit_sphere_z, Sphere sphere);
  // set by Sphere type, more scriptable from python
  void Set(double theta, double phi, Sphere sphere, bool radians = false);

  // The basic methods for extracting each of the angular coordinate values.
  // Survey coordinates.
  double Lambda();
  double Eta();
  double LambdaRadians();
  double EtaRadians();

  // Equatorial coordinates
  double RA();
  double DEC();
  double RARadians();
  double DECRadians();

  // Galactic coordinates
  double GalLon();
  double GalLat();
  double GalLonRadians();
  double GalLatRadians();

  // And the associated methods for doing the same with the unit sphere
  // Cartesian coordinates.
  double UnitSphereX();
  double UnitSphereY();
  double UnitSphereZ();

  // The previous methods give access to the internal representation that the
  // class uses to store an angular position.  Occasionally, however, it can be
  // useful to provide a Cartesian breakdown of the position for a particular
  // coordinate basis.
  double UnitSphereX(Sphere sphere);
  double UnitSphereY(Sphere sphere);
  double UnitSphereZ(Sphere sphere);

  // Once you have a single angular position on the sphere, very quickly you
  // end up wanting to know the angular distance between your coordinate and
  // another.  This method returns that value in degrees.
  double AngularDistance(AngularCoordinate& ang);
  double AngularDistance(AngularCoordinate* ang);

  // And these two methods return the dot-product and cross-product between
  // the unit vector represented by your angular position and another.
  double DotProduct(AngularCoordinate& ang);
  double DotProduct(AngularCoordinate* ang);
  AngularCoordinate CrossProduct(AngularCoordinate& ang,
				 Sphere sphere = Equatorial);
  AngularCoordinate CrossProduct(AngularCoordinate* ang,
				 Sphere sphere = Equatorial);

  // Given another AngularCoordinate, return the AngularCoordinate that defines
  // the great circle running through the current point and the input point.
  // This is different than just the cross-product in that we need to do the
  // calculation in the native coordinate system of the current point.
  void GreatCircle(AngularCoordinate& ang, AngularCoordinate& great_circle,
		   Sphere sphere = Equatorial);

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
  // return a copy of the rotated position.  The third variant is where we'll
  // actually do the work of the rotation, but the other two methods are the
  // preferred interfaces.
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle,
	      Sphere sphere = Equatorial);
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle,
	      AngularCoordinate& rotated_ang, Sphere sphere = Equatorial);
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle,
	      Sphere sphere, double& unit_sphere_x, double& unit_sphere_y,
	      double& unit_sphere_z);

  // Static methods for when users want to switch between coordinate systems
  // without instantiating the class.
  static void SurveyToGalactic(double lambda, double eta,
                               double& gal_lon, double& gal_lat,
			       bool radians = false);
  static void SurveyToEquatorial(double lambda, double eta,
                                 double& ra, double& dec,
				 bool radians = false);
  static void EquatorialToSurvey(double ra, double dec,
                                 double& lambda, double& eta,
				 bool radians = false);
  static void EquatorialToGalactic(double ra, double dec,
                                   double& gal_lon, double& gal_lat,
				   bool radians = false);
  static void GalacticToSurvey(double gal_lon, double gal_lat,
                               double& lambda, double& eta,
			       bool radians = false);
  static void GalacticToEquatorial(double gal_lon, double gal_lat,
                                   double& ra, double& dec,
				   bool radians = false);
  static void SurveyToXYZ(double lambda, double eta,
			  double& x, double& y, double& z,
			  bool radians = false);
  static void EquatorialToXYZ(double ra, double dec,
			      double& x, double& y, double& z,
			      bool radians = false);
  static void GalacticToXYZ(double gal_lon, double gal_lat,
			    double& x, double& y, double& z,
			    bool radians = false);

  // This is a bit more obscure.  The idea here is that, when you want to find
  // the pixel bounds that subtend a given angular scale about a point on the
  // sphere, finding those bounds in latitude is easier than in longitude.
  // Given a latitude in one of the coordinate systems, the multiplier tells
  // you how many more pixels you should check in the longitude direction
  // relative to the latitude direction.
  static double EtaMultiplier(double lam);
  static double RAMultiplier(double dec);
  static double GalLonMultiplier(double glat);

  // Finally, we provide a number of static methods for doing bulk conversion
  // of input vectors into AngularVectors and file I/O.  The idea here is to
  // provide a standard interface that users can adopt without having to
  // duplicate the same code over and over again.  Since most of the
  // astronomical community uses RA-DEC coordinates, that will be our default
  // as well.  Any failure to parse the input will return "false".  Like C++
  // arrays, the columns in an input file are zero-indexed, so the default
  // behavior is to read the first two columns of the input file as theta and
  // phi, respectively.
  static bool ToAngularVector(std::vector<double>& thetaVec,
			      std::vector<double>& phiVec,
			      AngularVector& ang,
			      Sphere sphere = Equatorial,
			      bool radians = false);
  static bool ToAngularVector(const std::string& input_file,
			      AngularVector& ang,
			      Sphere sphere = Equatorial,
			      bool radians = false,
			      uint8_t theta_column = 0,
			      uint8_t phi_column = 1);
  static bool FromAngularVector(AngularVector& ang,
				std::vector<double>& thetaVec,
				std::vector<double>& phiVec,
				Sphere sphere = Equatorial,
				bool radians = false);
  static bool FromAngularVector(AngularVector& ang,
				const std::string& output_file,
				Sphere sphere = Equatorial,
				bool radians = false);

 private:
  double us_x_, us_y_, us_z_;
};

class WeightedAngularCoordinate : public AngularCoordinate {
  // Sub-class of AngularCoordinate where we attach a weight value to that
  // angular position.

 public:
  WeightedAngularCoordinate();
  WeightedAngularCoordinate(double theta, double phi,
			    double weight, Sphere sphere = Survey,
			    bool radians = false);
  WeightedAngularCoordinate(double theta, double phi,
			    double weight, FieldDict& fields,
			    Sphere sphere = Survey,
			    bool radians = false);
  WeightedAngularCoordinate(double unit_sphere_x, double unit_sphere_y,
			    double unit_sphere_z, double weight);
  WeightedAngularCoordinate(double unit_sphere_x, double unit_sphere_y,
			    double unit_sphere_z, double weight,
			    FieldDict& fields);
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

  // Static methods for I/O to and from WAngularVectors.  For the input method
  // using a file, a negative value for the weight column results in unit
  // weights for all of the WAngularVector members.
  static bool ToWAngularVector(std::vector<double>& thetaVec,
			       std::vector<double>& phiVec,
			       std::vector<double>& weightVec,
			       WAngularVector& w_ang,
			       Sphere sphere = Equatorial,
			       bool radians = false);
  static bool ToWAngularVector(std::vector<double>& thetaVec,
			       std::vector<double>& phiVec,
			       double weight,
			       WAngularVector& w_ang,
			       Sphere sphere = Equatorial,
			       bool radians = false);
  static bool ToWAngularVector(const std::string& input_file,
			       WAngularVector& w_ang,
			       Sphere sphere = Equatorial,
			       bool radians = false,
			       uint8_t theta_column = 0,
			       uint8_t phi_column = 1,
			       int8_t weight_column = -1);
  static bool FromWAngularVector(WAngularVector& w_ang,
				 std::vector<double>& thetaVec,
				 std::vector<double>& phiVec,
				 std::vector<double>& weightVec,
				 Sphere sphere = Equatorial,
				 bool radians = false);
  static bool FromWAngularVector(WAngularVector& w_ang,
				 const std::string& output_file,
				 Sphere sphere = Equatorial,
				 bool radians = false);

  // For WeightedAngularCoordinate instances that have Field values, we have two
  // additional methods for reading from and writing to ascii files.  In both
  // cases the FieldColumnDict maps between Field name and column in the input
  // or output file.  Missing columns or Fields will result in 0.0 values in
  // the corresponding output.  If the columns assigned to field values
  // overlap those given for the angular coordinates or weight, then the
  // Field values will be output instead.
  static bool ToWAngularVector(const std::string& input_file,
			       WAngularVector& w_ang,
			       FieldColumnDict& field_columns,
			       Sphere sphere = Equatorial,
			       bool radians = false,
			       uint8_t theta_column = 0,
			       uint8_t phi_column = 1,
			       int8_t weight_column = -1);
  static bool FromWAngularVector(WAngularVector& w_ang,
				 FieldColumnDict& field_columns,
				 const std::string& output_file,
				 Sphere sphere = Equatorial,
				 bool radians = false,
				 uint8_t theta_column = 0,
				 uint8_t phi_column = 1,
				 uint8_t weight_column = 2);

  // We can also provide some wrapper methods for bulk adding of Field values
  // to WAngularVectors
  static bool AddField(WAngularVector& w_ang,
		       const std::string& field_name,
		       std::vector<double>& field_value);

 private:
  double weight_;
  FieldDict field_;
};

class CosmoCoordinate : public WeightedAngularCoordinate {
  // Sub-class of WeightedAngularCoordinate where we attach a redshift to the
  // angular position on the sphere.  This is instantiated with a redshift, so
  // the position is now a three dimensional coordinate.

 public:
  CosmoCoordinate();
  CosmoCoordinate(double theta, double phi, double redshift,
		  double weight, Sphere sphere = Survey, bool radians = false);
  CosmoCoordinate(double unit_sphere_x, double unit_sphere_y,
		  double unit_sphere_z, double redshift, double weight);
  ~CosmoCoordinate();

  // Since we have a redshift attached to the location, we can now attach some
  // functionality related to the 3-D coordinates.  Start with some conversions
  // between angular distance between an input point and our coordinate and
  // the projected radius at our coordinate's redshift.  Since we're getting
  // this data out of the Cosmology class, we follow that class's convention
  // of giving distance in comoving Mpc/h
  double ProjectedRadius(AngularCoordinate& ang);
  double ProjectedRadius(AngularCoordinate* ang);

  // We can aslo offer variations on the dot product that takes
  // CosmoCoordinates and do their calculations with the full 3-D vectors
  double DotProduct(CosmoCoordinate& ang);
  double DotProduct(CosmoCoordinate* ang);

  // Some methods for accessing the distance values implied by our redshift.
  double ComovingDistance();
  double AngularDiameterDistance();
  double LuminosityDistance();

  // Getter and setter for the redshift.
  double Redshift();
  void SetRedshift(double redshift);

  // Static methods for I/O to and from CosmoVectors.  For the input method
  // using a file, a negative value for the weight column results in unit
  // weights for all of the WAngularVector members.
  static bool ToCosmoVector(std::vector<double>& thetaVec,
			    std::vector<double>& phiVec,
			    std::vector<double>& redshiftVec,
			    std::vector<double>& weightVec,
			    CosmoVector& z_ang,
			    Sphere sphere = Equatorial,
			    bool radians = false);
  static bool ToCosmoVector(std::vector<double>& thetaVec,
			    std::vector<double>& phiVec,
			    std::vector<double>& redshiftVec,
			    double weight,
			    CosmoVector& z_ang,
			    Sphere sphere = Equatorial,
			    bool radians = false);
  static bool ToCosmoVector(const std::string& input_file,
			    CosmoVector& z_ang,
			    Sphere sphere = Equatorial,
			    bool radians = false,
			    uint8_t theta_column = 0,
			    uint8_t phi_column = 1,
			    uint8_t redshift_column = 2,
			    int8_t weight_column = -1);
  static bool FromCosmoVector(CosmoVector& z_ang,
			      std::vector<double>& thetaVec,
			      std::vector<double>& phiVec,
			      std::vector<double>& redshiftVec,
			      std::vector<double>& weightVec,
			      Sphere sphere = Equatorial,
			      bool radians = false);
  static bool FromCosmoVector(CosmoVector& z_ang,
			      const std::string& output_file,
			      Sphere sphere = Equatorial,
			      bool radians = false);

 private:
  double redshift_;
};

class IndexedAngularCoordinate : public AngularCoordinate {
  // Sub-class of AngularCoordinate where we attach an index value to the
  // position on the sphere.  The mathematics of what we would do with an
  // index value versus the floating point values in WeightedAngularCoordinate
  // necessitate a separate derived class for integers.

 public:
  IndexedAngularCoordinate();
  IndexedAngularCoordinate(double theta, double phi,
			   uint32_t index, Sphere sphere = Survey,
			   bool radians = false);
  IndexedAngularCoordinate(double unit_sphere_x, double unit_sphere_y,
			    double unit_sphere_z, uint32_t index);
  ~IndexedAngularCoordinate();

  void SetIndex(uint32_t index);
  uint32_t Index();

  // Static methods for I/O to and from IndexedVectors.  For the input method
  // using a file, a negative value for the weight column results in the line
  // number being used as the index.  Likewise for the vector-based methods.
  static bool ToIAngularVector(std::vector<double>& thetaVec,
			       std::vector<double>& phiVec,
			       std::vector<uint32_t>& indexVec,
			       IAngularVector& i_ang,
			       Sphere sphere = Equatorial,
			       bool radians = false);
  static bool ToIAngularVector(std::vector<double>& thetaVec,
			       std::vector<double>& phiVec,
			       IAngularVector& i_ang,
			       Sphere sphere = Equatorial,
			       bool radians = false);
  static bool ToIAngularVector(const std::string& input_file,
			       IAngularVector& i_ang,
			       Sphere sphere = Equatorial,
			       bool radians = false,
			       uint8_t theta_column = 0,
			       uint8_t phi_column = 1,
			       int8_t index_column = -1);
  static bool FromIAngularVector(IAngularVector& i_ang,
				 std::vector<double>& thetaVec,
				 std::vector<double>& phiVec,
				 std::vector<uint32_t>& indexVec,
				 Sphere sphere = Equatorial,
				 bool radians = false);
  static bool FromIAngularVector(IAngularVector& i_ang,
				 const std::string& output_file,
				 Sphere sphere = Equatorial,
				 bool radians = false);

 private:
  uint32_t index_;
};

} // end namespace Stomp

#endif
