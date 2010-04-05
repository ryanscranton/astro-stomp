%module stomp
%{
#include "stomp_core.h"
#include "stomp_angular_bin.h"
#include "stomp_angular_coordinate.h"
#include "stomp_angular_correlation.h"
#include "stomp_pixel.h"
#include "stomp_scalar_pixel.h"
#include "stomp_tree_pixel.h"
#include "stomp_base_map.h"
#include "stomp_map.h"
#include "stomp_scalar_map.h"
#include "stomp_tree_map.h"
#include "stomp_geometry.h"
#include "stomp_util.h"
%}

%include stl.i
%include std_string.i
%include std_pair.i
%include std_map.i
%include std_vector.i

// Need to add AngularCoordinate, WeightedAngularCoordinate, CosmoCoordinate
// TreePixel and TreeMap manually to avoid some over-loading issues.

namespace Stomp {

class Pixel;                // class declaration in stomp_pixel.h
class AngularBin;           // class definition in stomp_angular_bin.h
class AngularCorrelation;   // class definition in stomp_angular_correlation.h
class AngularCoordinate;
class WeightedAngularCoordinate;
class CosmoCoordinate;
class TreePixel;
class TreeNeighbor;
class NearestNeighborPixel;
class NearestNeighborPoint;

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

typedef std::vector<CosmoCoordinate> CosmoVector;
typedef CosmoVector::iterator CosmoIterator;
typedef std::vector<CosmoCoordinate *> CosmoPtrVector;
typedef CosmoPtrVector::iterator CosmoPtrIterator;

typedef std::pair<double, TreePixel*> DistancePixelPair;
typedef std::priority_queue<DistancePixelPair,
  std::vector<DistancePixelPair>, NearestNeighborPixel> PixelQueue;

typedef std::pair<double, WeightedAngularCoordinate*> DistancePointPair;
typedef std::priority_queue<DistancePointPair,
  std::vector<DistancePointPair>, NearestNeighborPoint> PointQueue;
}

%include "stomp_core.h"
%include "stomp_angular_bin.h"
%include "stomp_angular_correlation.h"
%include "stomp_pixel.h"
%include "stomp_scalar_pixel.h"
%include "stomp_base_map.h"
%include "stomp_map.h"
%include "stomp_scalar_map.h"
%include "stomp_geometry.h"
%include "stomp_util.h"

namespace Stomp {

class AngularCoordinate {
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

  void SetSurveyCoordinates(double lambda, double eta);
  void SetEquatorialCoordinates(double ra, double dec);
  void SetGalacticCoordinates(double gal_lon, double gal_lat);
  void SetUnitSphereCoordinates(double unit_sphere_x, double unit_sphere_y,
				double unit_sphere_z);
  double Lambda();
  double Eta();
  double RA();
  double DEC();
  double GalLon();
  double GalLat();
  double UnitSphereX();
  double UnitSphereY();
  double UnitSphereZ();
  double AngularDistance(AngularCoordinate& ang);
  double DotProduct(AngularCoordinate& ang);
  AngularCoordinate CrossProduct(AngularCoordinate& ang);
  void GreatCircle(AngularCoordinate& ang, AngularCoordinate& great_circle);
  double PositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double PositionAngle(Pixel& pix, Sphere sphere = Equatorial);
  double CosPositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double CosPositionAngle(Pixel& pix, Sphere sphere = Equatorial);
  double SinPositionAngle(AngularCoordinate& ang, Sphere sphere = Equatorial);
  double SinPositionAngle(Pixel& pix, Sphere sphere = Equatorial);
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle);
  void Rotate(AngularCoordinate& fixed_ang, double rotation_angle,
	      AngularCoordinate& rotated_ang);
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
  static double EtaMultiplier(double lam);
  static double RAMultiplier(double dec);
  static double GalLonMultiplier(double glat);
};

class WeightedAngularCoordinate : public AngularCoordinate {
 public:
  WeightedAngularCoordinate();
  WeightedAngularCoordinate(double theta, double phi,
			    double weight, Sphere sphere = Survey);
  WeightedAngularCoordinate(double unit_sphere_x, double unit_sphere_y,
			    double unit_sphere_z, double weight);
  ~WeightedAngularCoordinate();

  void SetWeight(double weight);
  double Weight();
  void SetField(const std::string& field_name, double weight);
  double Field(const std::string& field_name);
  uint16_t NFields();
  bool HasFields();
  void FieldNames(std::vector<std::string>& field_names);
  FieldIterator FieldBegin();
  FieldIterator FieldEnd();
  void CopyFields(WeightedAngularCoordinate& w_ang);
  void CopyFieldToWeight(const std::string& field_name);
  void RestoreOriginalWeight();
};

class CosmoCoordinate : public WeightedAngularCoordinate {
 public:
  CosmoCoordinate();
  CosmoCoordinate(double theta, double phi, double weight,
		  double redshift, Sphere sphere = Survey);
  CosmoCoordinate(double unit_sphere_x, double unit_sphere_y,
		  double unit_sphere_z, double weight, double redshift);
  ~CosmoCoordinate();

  double ProjectedRadius(AngularCoordinate& ang);
  double DotProduct(CosmoCoordinate& ang);
  double ComovingDistance();
  double AngularDiameterDistance();
  double LuminosityDistance();
  double Redshift();
  void SetRedshift(double redshift);
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
  void _AddSubNodes(uint16_t& n_nodes);
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

class TreeMap : public BaseMap {
 public:
  friend class NearestNeighborPixel;
  TreeMap(uint32_t resolution=HPixResolution,
	  uint16_t maximum_points=50);
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
  bool AddPoint(WeightedAngularCoordinate& w_ang);
  bool AddPoint(AngularCoordinate& ang, double object_weight = 1.0);
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

} // end namespace Stomp

