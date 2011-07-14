%module stomp
%{
#include "../stomp/stomp_core.h"
#include "../stomp/stomp_angular_bin.h"
#include "../stomp/stomp_radial_bin.h"
#include "../stomp/stomp_angular_coordinate.h"
#include "../stomp/stomp_angular_correlation.h"
#include "../stomp/stomp_pixel.h"
#include "../stomp/stomp_scalar_pixel.h"
#include "../stomp/stomp_tree_pixel.h"
#include "../stomp/stomp_itree_pixel.h"
#include "../stomp/stomp_base_map.h"
#include "../stomp/stomp_map.h"
#include "../stomp/stomp_scalar_map.h"
#include "../stomp/stomp_tree_map.h"
#include "../stomp/stomp_itree_map.h"
#include "../stomp/stomp_geometry.h"
#include "../stomp/stomp_util.h"
%}

// catch these types of exceptions in python exceptions
%typemap(throws) const char * %{
    PyErr_SetString(PyExc_RuntimeError, $1);
    SWIG_fail;
%}

%feature("autodoc", "1");

%include stdint.i
%include stl.i
%include std_string.i
%include std_pair.i
%include std_map.i
%include std_vector.i
%include generators.i

%include typemaps.i
%apply double& INPUT { double& }
#%apply double *OUTPUT { double& }
%apply uint32_t& INPUT { uint32_t& }

// Need to add AngularCoordinate, WeightedAngularCoordinate, CosmoCoordinate
// TreePixel and TreeMap manually to avoid some over-loading issues.

namespace Stomp {

class Pixel;
class PixelOrdering;
class AngularBin;
class RadialBin;
class AngularCorrelation;
class AngularCoordinate;
class WeightedAngularCoordinate;
class CosmoCoordinate;
class IndexedAngularCoordinate;
class TreePixel;
class IndexedTreePixel;
class TreeNeighbor;
class IndexedTreeNeighbor;
class NearestNeighborPixel;
class NearestNeighborIndexedPixel;
class NearestNeighborPoint;
class NearestNeighborIndexedPoint;

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

typedef std::pair<double, TreePixel*> DistancePixelPair;
typedef std::priority_queue<DistancePixelPair,
  std::vector<DistancePixelPair>, NearestNeighborPixel> PixelQueue;

typedef std::pair<double, IndexedTreePixel*> DistanceIPixelPair;
typedef std::priority_queue<DistanceIPixelPair,
  std::vector<DistanceIPixelPair>, NearestNeighborIndexedPixel> IPixelQueue;

typedef std::pair<double, WeightedAngularCoordinate*> DistancePointPair;
typedef std::priority_queue<DistancePointPair,
  std::vector<DistancePointPair>, NearestNeighborPoint> PointQueue;

typedef std::pair<double, IndexedAngularCoordinate*> DistanceIPointPair;
typedef std::priority_queue<DistanceIPointPair,
  std::vector<DistanceIPointPair>, NearestNeighborIndexedPoint> IPointQueue;

typedef std::vector<uint32_t> IndexVector;
typedef IndexVector::iterator IndexIterator;

} # end namespace Stomp

%include "../stomp/stomp_core.h"
%include "../stomp/stomp_angular_bin.h"
%include "../stomp/stomp_radial_bin.h"
%include "../stomp/stomp_angular_correlation.h"
%include "../stomp/stomp_pixel.h"
%include "../stomp/stomp_scalar_pixel.h"
%include "../stomp/stomp_base_map.h"
%include "../stomp/stomp_map.h"
%include "../stomp/stomp_scalar_map.h"
%include "../stomp/stomp_geometry.h"
%include "../stomp/stomp_util.h"

namespace Stomp {

class AngularCoordinate {
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

  void SetSurveyCoordinates(double lambda, double eta, bool radians = false);
  void SetEquatorialCoordinates(double ra, double dec, bool radians = false);
  void SetGalacticCoordinates(double gal_lon, double gal_lat,
                              bool radians = false);
  void SetUnitSphereCoordinates(double unit_sphere_x, double unit_sphere_y,
                                double unit_sphere_z);
  void SetUnitSphereCoordinates(double unit_sphere_x, double unit_sphere_y,
                                double unit_sphere_z, Sphere sphere);
  void Set(double theta, double phi, Sphere sphere, bool radians = false);
  double Lambda();
  double Eta();
  double LambdaRadians();
  double EtaRadians();
  double RA();
  double DEC();
  double RARadians();
  double DECRadians();
  double GalLon();
  double GalLat();
  double GalLonRadians();
  double GalLatRadians();
  double UnitSphereX();
  double UnitSphereY();
  double UnitSphereZ();
  double UnitSphereX(Sphere sphere);
  double UnitSphereY(Sphere sphere);
  double UnitSphereZ(Sphere sphere);
  double AngularDistance(AngularCoordinate& ang);
  double DotProduct(AngularCoordinate& ang);
  AngularCoordinate CrossProduct(AngularCoordinate& ang,
                                 Sphere sphere = Equatorial);
  void GreatCircle(AngularCoordinate& ang, AngularCoordinate& great_circle,
                   Sphere sphere = Equatorial);
  double PositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double PositionAngle(Pixel& pix, Sphere sphere = Equatorial);
  double CosPositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double CosPositionAngle(Pixel& pix, Sphere sphere = Equatorial);
  double SinPositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double SinPositionAngle(Pixel& pix, Sphere sphere = Equatorial);
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle,
              Sphere sphere = Equatorial);
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle,
              AngularCoordinate& rotated_ang, Sphere sphere = Equatorial);
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle,
              Sphere sphere, double& unit_sphere_x, double& unit_sphere_y,
              double& unit_sphere_z);
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
  static double EtaMultiplier(double lam);
  static double RAMultiplier(double dec);
  static double GalLonMultiplier(double glat);
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
};

class WeightedAngularCoordinate : public AngularCoordinate {
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
  void SetWeight(double weight);
  double Weight();
  void SetField(const std::string& field_name, double weight);
  double Field(const std::string& field_name);
  uint16_t NFields();
  bool HasFields();
  void FieldNames(std::vector<std::string>& field_names);
  void CopyFields(WeightedAngularCoordinate& w_ang);
  void CopyFieldToWeight(const std::string& field_name);
  void RestoreOriginalWeight();
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
  static bool AddField(WAngularVector& w_ang,
                       const std::string& field_name,
                       std::vector<double>& field_value);
};

class CosmoCoordinate : public WeightedAngularCoordinate {
 public:
  CosmoCoordinate();
  CosmoCoordinate(double theta, double phi, double redshift,
                  double weight, Sphere sphere = Survey, bool radians = false);
  CosmoCoordinate(double unit_sphere_x, double unit_sphere_y,
                  double unit_sphere_z, double redshift, double weight);
  ~CosmoCoordinate();
  double ProjectedRadius(AngularCoordinate& ang);
  double DotProduct(CosmoCoordinate& ang);
  double ComovingDistance();
  double AngularDiameterDistance();
  double LuminosityDistance();
  double Redshift();
  void SetRedshift(double redshift);
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
};

class IndexedAngularCoordinate : public AngularCoordinate {
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
};

class TreePixel : public Pixel {
 public:
  friend class NearestNeighborPixel;
  TreePixel();
  TreePixel(const uint32_t resolution, const uint32_t pixnum,
	    const uint16_t maximum_points=200);
  TreePixel(AngularCoordinate& ang, const uint32_t resolution,
	    const uint16_t maximum_points=200);
  TreePixel(const uint32_t x, const uint32_t y, const uint32_t resolution,
	    const uint16_t maximum_points=200);
  virtual ~TreePixel();

  bool _InitializeSubPixels();
  uint32_t DirectPairCount(AngularCoordinate& ang, AngularBin& theta,
			   int16_t region = -1);
  uint32_t FindPairs(AngularCoordinate& ang, AngularBin& theta,
		     int16_t region = -1);
  uint32_t FindPairs(AngularCoordinate& ang,
		     double theta_min, double theta_max);
  uint32_t FindPairs(AngularCoordinate& ang, double theta_max);
  double DirectWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			     int16_t region = -1);
  double FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			   int16_t region = -1);
  double FindWeightedPairs(AngularCoordinate& ang,
			   double theta_min, double theta_max);
  double FindWeightedPairs(AngularCoordinate& ang, double theta_max);
  double DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
			     AngularBin& theta, int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   AngularBin& theta, int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_min, double theta_max);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_max);
  void FindPairs(AngularVector& ang, AngularBin& theta,
		 int16_t region = -1);
  void FindPairs(AngularVector& ang, AngularCorrelation& wtheta,
		 int16_t region = -1);
  void FindWeightedPairs(AngularVector& ang, AngularBin& theta,
			 int16_t region = -1);
  void FindWeightedPairs(AngularVector& ang, AngularCorrelation& wtheta,
			 int16_t region = -1);
  void FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta,
			 int16_t region = -1);
  void FindWeightedPairs(WAngularVector& w_ang,
			 AngularCorrelation& wtheta, int16_t region = -1);
  double DirectWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			     const std::string& field_name,
			     int16_t region = -1);
  double FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			   const std::string& field_name, int16_t region = -1);
  double FindWeightedPairs(AngularCoordinate& ang,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(AngularCoordinate& ang, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(AngularVector& ang, AngularBin& theta,
			 const std::string& field_name, int16_t region = -1);
  void FindWeightedPairs(AngularVector& ang, AngularCorrelation& wtheta,
			 const std::string& field_name, int16_t region = -1);
  double DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
			     AngularBin& theta, const std::string& field_name,
			     int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang, AngularBin& theta,
			   const std::string& field_name, int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta,
			 const std::string& field_name, int16_t region = -1);
  void FindWeightedPairs(WAngularVector& w_ang, AngularCorrelation& wtheta,
			 const std::string& field_name, int16_t region = -1);
  double DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
			     const std::string& ang_field_name,
			     AngularBin& theta, const std::string& field_name,
			     int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name, AngularBin& theta,
			   const std::string& field_name, int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang,
			 const std::string& ang_field_name,
			 AngularBin& theta, const std::string& field_name,
			 int16_t region = -1);
  void FindWeightedPairs(WAngularVector& w_ang,
			 const std::string& ang_field_name,
			 AngularCorrelation& wtheta,
			 const std::string& field_name, int16_t region = -1);
  uint16_t FindKNearestNeighbors(AngularCoordinate& ang, uint8_t n_neighbors,
				 WAngularVector& neighbors_ang);
  uint16_t FindNearestNeighbor(AngularCoordinate& ang,
			       WeightedAngularCoordinate& neighbor_ang);
  double KNearestNeighborDistance(AngularCoordinate& ang, uint8_t n_neighbors,
				  uint16_t& nodes_visited);
  double NearestNeighborDistance(AngularCoordinate& ang,
				 uint16_t& nodes_visited);
  bool ClosestMatch(AngularCoordinate& ang, double max_distance,
		    WeightedAngularCoordinate& match_ang);
  void InitializeCorners();
  bool AddPoint(WeightedAngularCoordinate& w_ang);
  bool AddPoint(AngularCoordinate& ang, double object_weight = 1.0);
  uint32_t NPoints();
  uint32_t NPoints(Pixel& pix);
  double PixelWeight(Pixel& pix);
  double Coverage();
  double Coverage(Pixel& pix);
  void Points(WAngularVector& w_ang);
  void Points(WAngularVector& w_ang, Pixel& pix);
  uint16_t Nodes();
  void AddToWeight(double weight);
  double FieldTotal(const std::string& field_name);
  double FieldTotal(const std::string& field_name, Pixel& pix);
  void AddToField(const std::string& field_name, double weight);
  uint16_t NField();
  bool HasFields();
  void FieldNames(std::vector<std::string>& field_names);
  void SetPixelCapacity(uint16_t maximum_points);
  uint16_t PixelCapacity();
  bool HasPoints();
  bool HasNodes();
  void Clear();
  virtual double UnitSphereX();
  virtual double UnitSphereY();
  virtual double UnitSphereZ();
  virtual double UnitSphereX_UL();
  virtual double UnitSphereY_UL();
  virtual double UnitSphereZ_UL();
  virtual double UnitSphereX_UR();
  virtual double UnitSphereY_UR();
  virtual double UnitSphereZ_UR();
  virtual double UnitSphereX_LL();
  virtual double UnitSphereY_LL();
  virtual double UnitSphereZ_LL();
  virtual double UnitSphereX_LR();
  virtual double UnitSphereY_LR();
  virtual double UnitSphereZ_LR();
  virtual void WithinAnnulus(AngularBin& theta, PixelVector& pix,
			     bool check_full_pixel);
};

class NearestNeighborPixel {
 public:
  int operator()(const DistancePixelPair& x, const DistancePixelPair& y) {
    return x.first > y.first;
  }
};

class NearestNeighborPoint {
 public:
  int operator()(const DistancePointPair& x, const DistancePointPair& y) {
    return x.first < y.first;
  }
};

class TreeNeighbor {
 public:
  friend class NearestNeighborPoint;
  TreeNeighbor(AngularCoordinate& reference_ang,
	       uint8_t n_neighbors = 1);
  ~TreeNeighbor();

  void NearestNeighbors(WAngularVector& w_ang, bool save_neighbors = true);
  uint8_t Neighbors();
  uint8_t MaxNeighbors();
  bool TestPoint(WeightedAngularCoordinate* test_ang);
  double MaxDistance();
  double MaxAngularDistance();
  uint16_t NodesVisited();
  void AddNode();
};

class IndexedTreePixel : public Pixel {
 public:
  friend class NearestNeighborIndexedPixel;
  IndexedTreePixel();
  IndexedTreePixel(const uint32_t resolution, const uint32_t pixnum,
                   const uint16_t maximum_points=200);
  IndexedTreePixel(AngularCoordinate& ang, const uint32_t resolution,
                   const uint16_t maximum_points=200);
  IndexedTreePixel(const uint32_t x, const uint32_t y,
                   const uint32_t resolution,
                   const uint16_t maximum_points=200);
  virtual ~IndexedTreePixel();

  void FindPairs(AngularCoordinate& ang, AngularBin& theta,
                 IAngularVector& i_angVec);
  void FindPairs(AngularCoordinate& ang, AngularBin& theta,
                 IndexVector& pair_indices);
  void FindPairs(AngularCoordinate& ang,
                 double theta_min, double theta_max,
                 IAngularVector& i_angVec);
  void FindPairs(AngularCoordinate& ang,
                 double theta_min, double theta_max,
                 IndexVector& pair_indices);
  void FindPairs(AngularCoordinate& ang, double theta_max,
                 IAngularVector& pair_indices);
  void FindPairs(AngularCoordinate& ang, double theta_max,
                 IndexVector& pair_indices);
  uint16_t FindKNearestNeighbors(AngularCoordinate& ang, uint8_t n_neighbors,
                                 IAngularVector& neighbors_ang);

  uint16_t FindNearestNeighbor(AngularCoordinate& ang,
                               IndexedAngularCoordinate& neighbor_ang);
  double KNearestNeighborDistance(AngularCoordinate& ang, uint8_t n_neighbors,
                                  uint16_t& nodes_visited);
  double NearestNeighborDistance(AngularCoordinate& ang,
                                 uint16_t& nodes_visited);
  bool ClosestMatch(AngularCoordinate& ang, double max_distance,
                    IndexedAngularCoordinate& match_ang);
  void InitializeCorners();
  bool AddPoint(IndexedAngularCoordinate& w_ang);
  bool AddPoint(AngularCoordinate& ang, uint32_t index);
  uint32_t NPoints();
  uint32_t NPoints(Pixel& pix);
  void Indices(Pixel& pix, IndexVector& indices);
  double Coverage();
  double Coverage(Pixel& pix);
  void Points(IAngularVector& i_ang);
  void Points(IAngularVector& i_ang, Pixel& pix);
  uint16_t Nodes();
  void SetPixelCapacity(uint16_t maximum_points);
  uint16_t PixelCapacity();
  bool HasPoints();
  bool HasNodes();
  void Clear();
  virtual double UnitSphereX();
  virtual double UnitSphereY();
  virtual double UnitSphereZ();
  virtual double UnitSphereX_UL();
  virtual double UnitSphereY_UL();
  virtual double UnitSphereZ_UL();
  virtual double UnitSphereX_UR();
  virtual double UnitSphereY_UR();
  virtual double UnitSphereZ_UR();
  virtual double UnitSphereX_LL();
  virtual double UnitSphereY_LL();
  virtual double UnitSphereZ_LL();
  virtual double UnitSphereX_LR();
  virtual double UnitSphereY_LR();
  virtual double UnitSphereZ_LR();
  virtual void WithinAnnulus(AngularBin& theta, PixelVector& pix,
                             bool check_full_pixel);
};

class NearestNeighborIndexedPixel {
 public:
  int operator()(const DistanceIPixelPair& x, const DistanceIPixelPair& y) {
    return x.first > y.first;
  }
};

class NearestNeighborIndexedPoint {
 public:
  int operator()(const DistanceIPointPair& x, const DistanceIPointPair& y) {
    return x.first < y.first;
  }
};

class IndexedTreeNeighbor {
 public:
  friend class NearestNeighborIndexedPoint;
  IndexedTreeNeighbor(AngularCoordinate& reference_ang,
               uint8_t n_neighbors = 1);
  IndexedTreeNeighbor(AngularCoordinate& reference_ang,
               uint8_t n_neighbors, double max_distance);
  ~IndexedTreeNeighbor();

  void NearestNeighbors(IAngularVector& i_ang, bool save_neighbors = true);
  uint8_t Neighbors();
  uint8_t MaxNeighbors();
  bool TestPoint(IndexedAngularCoordinate* test_ang);
  double MaxDistance();
  double MaxAngularDistance();
  uint16_t NodesVisited();
  void AddNode();
};

class TreeMap : public BaseMap {
 public:
  friend class NearestNeighborPixel;
  TreeMap(uint32_t resolution=HPixResolution,
	  uint16_t maximum_points=50);
  TreeMap(const std::string& input_file,
          uint32_t resolution=HPixResolution, uint16_t maximum_points=50,
          AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
          bool verbose = false, uint8_t theta_column = 0,
          uint8_t phi_column = 1, int8_t weight_column = -1);
  TreeMap(const std::string& input_file, FieldColumnDict& field_columns,
          uint32_t resolution=HPixResolution, uint16_t maximum_points=50,
          AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
          bool verbose = false, uint8_t theta_column = 0,
          uint8_t phi_column = 1, int8_t weight_column = -1);
  ~TreeMap();

  uint32_t FindPairs(AngularCoordinate& ang, AngularBin& theta);
  uint32_t FindPairs(AngularCoordinate& ang,
		     double theta_min, double theta_max);
  uint32_t FindPairs(AngularCoordinate& ang, double theta_max);
  void FindPairs(AngularVector& ang, AngularBin& theta);
  void FindPairs(AngularVector& ang, AngularCorrelation& wtheta);
  double FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta);
  double FindWeightedPairs(AngularCoordinate& ang,
			   double theta_min, double theta_max);
  double FindWeightedPairs(AngularCoordinate& ang, double theta_max);
  void FindWeightedPairs(AngularVector& ang, AngularBin& theta);
  void FindWeightedPairs(AngularVector& ang, AngularCorrelation& wtheta);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   AngularBin& theta);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_min, double theta_max);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_max);
  void FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta);
  void FindWeightedPairs(WAngularVector& w_ang,
			 AngularCorrelation& wtheta);
  double FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			   const std::string& field_name);
  double FindWeightedPairs(AngularCoordinate& ang,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(AngularCoordinate& ang, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(AngularVector& ang, AngularBin& theta,
			 const std::string& field_name);
  void FindWeightedPairs(AngularVector& ang, AngularCorrelation& wtheta,
			 const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang, AngularBin& theta,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta,
			 const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang, AngularCorrelation& wtheta,
			 const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name, AngularBin& theta,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang,
			 const std::string& ang_field_name,
			 AngularBin& theta, const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang,
			 const std::string& ang_field_name,
			 AngularCorrelation& wtheta,
			 const std::string& field_name);
  void FindPairsWithRegions(AngularVector& ang, AngularBin& theta);
  void FindPairsWithRegions(AngularVector& ang, AngularCorrelation& wtheta);
  void FindWeightedPairsWithRegions(AngularVector& ang, AngularBin& theta);
  void FindWeightedPairsWithRegions(AngularVector& ang,
                                    AngularCorrelation& wtheta);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang, AngularBin& theta);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang,
                                    AngularCorrelation& wtheta);
  void FindWeightedPairsWithRegions(AngularVector& ang, AngularBin& theta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(AngularVector& ang,
                                    AngularCorrelation& wtheta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang, AngularBin& theta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang,
                                    AngularCorrelation& wtheta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang,
                                    const std::string& ang_field_name,
                                    AngularBin& theta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang,
                                    const std::string& ang_field_name,
                                    AngularCorrelation& wtheta,
                                    const std::string& field_name);
  uint16_t FindKNearestNeighbors(AngularCoordinate& ang, uint8_t n_neighbors,
				 WAngularVector& neighbors_ang);
  uint16_t FindNearestNeighbor(AngularCoordinate& ang,
			       WeightedAngularCoordinate& neighbor_ang);
  double KNearestNeighborDistance(AngularCoordinate& ang, uint8_t n_neighbors,
				  uint16_t& nodes_visited);
  double NearestNeighborDistance(AngularCoordinate& ang,
				 uint16_t& nodes_visited);
  bool ClosestMatch(AngularCoordinate& ang, double max_distance,
		    WeightedAngularCoordinate& match_ang);
  bool AddPoint(WeightedAngularCoordinate& w_ang);
  bool AddPoint(AngularCoordinate& ang, double object_weight = 1.0);
  bool Read(const std::string& input_file,
	    AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
	    bool verbose = false, uint8_t theta_column = 0,
	    uint8_t phi_column = 1, int8_t weight_column = -1);
  bool Read(const std::string& input_file, FieldColumnDict& field_columns,
	    AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
	    bool verbose = false, uint8_t theta_column = 0,
	    uint8_t phi_column = 1, int8_t weight_column = -1);
  virtual void Coverage(PixelVector& superpix,
			uint32_t resolution = HPixResolution);
  bool Covering(Map& stomp_map, uint32_t maximum_pixels);
  virtual double FindUnmaskedFraction(Pixel& pix);
  virtual int8_t FindUnmaskedStatus(Pixel& pix);
  void NodeMap(Map& stomp_map);
  uint32_t Resolution();
  uint16_t PixelCapacity();
  void SetResolution(uint32_t resolution);
  void SetPixelCapacity(int pixel_capacity);
  uint32_t NPoints(uint32_t k = MaxPixnum);
  uint32_t NPoints(Pixel& pix);
  void Points(WAngularVector& w_ang);
  void Points(WAngularVector& w_ang, Pixel& pix);
  double Weight(uint32_t k = MaxPixnum);
  double Weight(Pixel& pix);
  double FieldTotal(const std::string& field_name,
			   uint32_t k = MaxPixnum);
  double FieldTotal(const std::string& field_name, Pixel& pix);
  uint16_t NField();
  bool HasFields();
  void FieldNames(std::vector<std::string>& field_names);
  uint16_t BaseNodes();
  uint16_t Nodes();
  virtual uint32_t Size();
  virtual double Area();
  void CalculateArea();
  virtual uint32_t MinResolution();
  virtual uint32_t MaxResolution();
  virtual uint8_t MinLevel();
  virtual uint8_t MaxLevel();
  virtual bool Empty();
  virtual void Clear();
};


class IndexedTreeMap : public BaseMap {
 public:
  friend class NearestNeighborPixel;
  IndexedTreeMap(uint32_t resolution=HPixResolution,
                 uint16_t maximum_points=50);
  IndexedTreeMap(
    const std::string& input_file,
    uint32_t resolution=HPixResolution, uint16_t maximum_points=50,
    AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
    bool verbose = false, uint8_t theta_column = 0,
    uint8_t phi_column = 1, int8_t index_column = -1);
  ~IndexedTreeMap();

  void FindPairs(AngularCoordinate& ang, AngularBin& theta,
                 IAngularVector& i_angVec);
  void FindPairs(AngularCoordinate& ang, AngularBin& theta,
                 IndexVector& pair_indices);
  void FindPairs(AngularCoordinate& ang,
                 double theta_min, double theta_max,
                 IAngularVector& i_angVec);
  void FindPairs(AngularCoordinate& ang,
                 double theta_min, double theta_max,
                 IndexVector& pair_indices);
  void FindPairs(AngularCoordinate& ang, double theta_max,
                 IAngularVector& pair_indices);
  void FindPairs(AngularCoordinate& ang, double theta_max,
                 IndexVector& pair_indices);
  uint16_t FindKNearestNeighbors(AngularCoordinate& ang, uint8_t n_neighbors,
                                 IAngularVector& neighbors_ang);
  uint16_t FindNearestNeighbor(AngularCoordinate& ang,
                               IndexedAngularCoordinate& neighbor_ang);
  double KNearestNeighborDistance(AngularCoordinate& ang, uint8_t n_neighbors,
                                  uint16_t& nodes_visited);
  double NearestNeighborDistance(AngularCoordinate& ang,
                                 uint16_t& nodes_visited);
  bool ClosestMatch(AngularCoordinate& ang, double max_distance,
                    IndexedAngularCoordinate& match_ang);
  bool AddPoint(IndexedAngularCoordinate& i_ang);
  bool AddPoint(AngularCoordinate& ang, uint32_t index);
  bool Read(const std::string& input_file,
            AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
            bool verbose = false, uint8_t theta_column = 0,
            uint8_t phi_column = 1, int8_t index_column = -1);
  virtual void Coverage(PixelVector& superpix,
                        uint32_t resolution = HPixResolution,
                        bool calculate_fraction = true);
  bool Covering(Map& stomp_map, uint32_t maximum_pixels);
  virtual double FindUnmaskedFraction(Pixel& pix);
  virtual int8_t FindUnmaskedStatus(Pixel& pix);
  void NodeMap(Map& stomp_map);
  uint32_t Resolution();
  uint16_t PixelCapacity();
  void SetResolution(uint32_t resolution);
  void SetPixelCapacity(int pixel_capacity);
  uint32_t NPoints(uint32_t k = MaxPixnum);
  uint32_t NPoints(Pixel& pix);
  void Points(IAngularVector& i_ang);
  void Points(IAngularVector& i_ang, Pixel& pix);
  void Indices(Pixel& pix, IndexVector& indices);
  uint16_t BaseNodes();
  uint16_t Nodes();
  virtual uint32_t Size();
  virtual double Area();
  void CalculateArea();
  virtual uint32_t MinResolution();
  virtual uint32_t MaxResolution();
  virtual uint8_t MinLevel();
  virtual uint8_t MaxLevel();
  virtual bool Empty();
  virtual void Clear();
};

} // end namespace Stomp

%template(AngularVector) std::vector<Stomp::AngularCoordinate>;
%template(ThetaVector) std::vector<Stomp::AngularBin>;
%template(RadialVector) std::vector<Stomp::RadialBin>;
%template(WAngularVector) std::vector<Stomp::WeightedAngularCoordinate>;
%template(CosmoVector) std::vector<Stomp::CosmoCoordinate>;
%template(IAngularVector) std::vector<Stomp::IndexedAngularCoordinate>;
%template(PixelVector) std::vector<Stomp::Pixel>;
%template(FieldDict) std::map<std::string, double>;
%template(FieldColumnDict) std::map<std::string, uint8_t>;
%template(DoubleVector) std::vector<double>;
%template(IndexVector) std::vector<uint32_t>;

SETUP_GENERATOR(std::vector<Stomp::AngularBin>::const_iterator)
ADD_GENERATOR(Stomp::AngularCorrelation, Bins,
std::vector<Stomp::AngularBin>::const_iterator, Stomp::AngularBin, Begin, End)

SETUP_GENERATOR(std::vector<Stomp::HistogramBin>::const_iterator)
ADD_GENERATOR(Stomp::Histogram, Bins,
std::vector<Stomp::HistogramBin>::const_iterator, Stomp::HistogramBin,
Begin, End)
