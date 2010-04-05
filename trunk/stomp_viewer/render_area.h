#ifndef RENDERAREA_H
#define RENDERAREA_H

#include <QBrush>
#include <QPen>
#include <QPixmap>
#include <QWidget>
#include <QPolygonF>
#include <QMouseEvent>
#include <stomp.h>

class QColor;
class Palette;
class RenderArea;

class Palette {
 public:
  enum PaletteType {
    BlueTemperature,
    GreenTemperature,
    RedTemperature,
    GrayScale,
    InverseGrayScale,
    Rainbow,
    InverseRainbow
  };
  Palette();
  Palette(PaletteType palette);
  ~Palette();
  void Initialize(PaletteType palette);
  void InitializeBlueTemperature();
  void InitializeGreenTemperature();
  void InitializeRedTemperature();
  void InitializeGrayScale();
  void InitializeInverseGrayScale();
  void InitializeRainbow();
  void InitializeInverseRainbow();
  QColor Color(double weight);
  std::string CurrentPalette();
  PaletteType CurrentPaletteType();

 private:
  PaletteType palette_type_;
  std::vector<QColor> rgb_;
};

class RenderArea : public QWidget {
  Q_OBJECT

 public:
  RenderArea(QWidget *parent = 0);

  QSize minimumSizeHint() const;
  QSize sizeHint() const;
  double weightMin();
  double weightMax();
  double pointsWeightMin();
  double pointsWeightMax();
  double longitudeMin();
  double longitudeMax();
  double latitudeMin();
  double latitudeMax();
  uint32_t maxResolution();
  bool autoUpdating();
  bool fullSky();
  bool displayingCoordinates();
  bool displayingGrid();
  bool displayingPoints();
  bool filteringPoints();

 public slots:
  // We're storing the map image in a QPixmap so that paint events will be
  // more graceful.  However, if we do something that should change the map,
  // we need to update the underlying map image.
  void updatePixmap(bool send_update_call = true);

  // Likewise, if we alter the image, then we need to re-draw the QPixmap that
  // contains the coordinate grid and the one that contains any points that have
  // been loaded.
  void updatePoints(bool send_update_call = true);
  void updateGrid(bool send_update_call = true);

  // Read in a new map.  Return true if the operation was successful.
  bool readNewMap(QString& file_name);

  // Read in a new set of points.  Return true if the operation was successful.
  bool readNewPointsFile(QString& file_name,
			 Stomp::AngularCoordinate::Sphere sphere,
			 bool weighted_points);

  // Remove either the current Map or set of points.
  void clearMap();
  void clearPoints();

  // Write the current rendered image to a PNG file.  Return true on success.
  bool writeToPng(QString& file_name);

  // Auto-alias the drawing of our image.
  void setAntialiased(bool antialiased);

  // Make the display automatically update after any input fields are changed.
  void setAutoUpdate(bool auto_update);

  // In all cases, we're mapping from some scalar weight to an
  // RGB flux (as determined by the Palette object).  With this method, we
  // can manually set the allowed weight range for our color mapping.  Any
  // weights outside of that range will be treated as equal to the corresponding
  // bound values.
  void setWeightRange(double weight_min, double weight_max);
  void setPointsWeightRange(double weight_min, double weight_max);

  // Add a palette to our object to scale between weight and RGB color.
  void setPalette(Palette::PaletteType palette_type);
  void setPointsPalette(Palette::PaletteType palette_type);

  // In case we don't want to have our plotting area stretch all the way
  // to the edge of the image.
  void setVerticalBuffer(uint16_t buffer_pixels);
  void setHorizontalBuffer(uint16_t buffer_pixels);
  void setTopBuffer(uint16_t buffer_pixels);
  void setBottomBuffer(uint16_t buffer_pixels);
  void setLeftBuffer(uint16_t buffer_pixels);
  void setRightBuffer(uint16_t buffer_pixels);

  // We can plot our maps in any of the available coordinate systems and these
  // methods allow us to choose accordingly.  If the bounds given exceed the
  // natural bounds of the coordinates, then they will be silently fixed.
  void setImageBounds(double lonmin, double lonmax,
		      double latmin, double latmax);
  void setImageBounds(double lonmin, double lonmax,
		      double latmin, double latmax,
		      Stomp::AngularCoordinate::Sphere sphere);
  void setSurveyBounds(double lammin, double lammax,
		       double etamin, double etamax);
  void setEquatorialBounds(double ramin, double ramax,
			   double decmin, double decmax);
  void setGalacticBounds(double lonmin, double lonmax,
			 double latmin, double latmax);
  void setFullSky(bool full_sky);
  void setCoordinateSystem(Stomp::AngularCoordinate::Sphere sphere);
  void useSurveyCoordinates();
  void useEquatorialCoordinates();
  void useGalacticCoordinates();

  // In order to help navigate the map presented, we have zoom in and zoom out
  // functionality (default is x2 zoom).
  void zoomIn();
  void zoomOut();
  void zoomIn(double zoom_factor);
  void zoomOut(double zoom_factor);

  // Fill the polygons. Or don't.
  void fillPolygons(bool fill_polygons);

  // Control whether the image is displayed with in either a Cartesian or an
  // Aitoff-Hammer projection.
  void useAitoffProjection(bool aitoff_projection);

  // Controls for the display and appearance of the coordinate grid.
  void setDisplayCoordinates(bool display_coordinates);
  void setDisplayGrid(bool display_grid);

  // Controls for displaying the points.
  void setDisplayPoints(bool display_points);
  void setFilterPoints(bool filter_points);

  // In general, when we make a polygon out of a pixel, we need to account
  // for the curvature of the pixel boundary in our coordinate space.  In
  // practice, this means turning the pixel into a many sided polygon, where
  // the number of vertices is a function of the difference between the
  // resolution of the pixel in question and the highest resolution pixels we
  // choose to include in our image.  This parameter controls that latter
  // limit (by default this is set to 2048 in the constructor).  If an invalid
  // resolution value is entered the maximum resolution is unchanged and the
  // return value is false.
  void setMaxResolution(uint32_t resolution);

  // At the same time, there is a second resolution that we want to consider,
  // that used to actually render the image.  This method sets that value.
  void setRenderResolution(uint32_t resolution);


 signals:
  void newMapParameters();

 protected:
  // paintEvent is an over-ridden method from the QWidget class that we use to
  // render our displayed Map/Points/Grid image.  paintEvent is called every
  // time the window system needs to re-draw the window, so it's called very
  // frequently.  To keep from having to re-render each time paintEvent is
  // called, we put each of the Map/Points/Grid images into a QPixmap which
  // paintEvent re-draws each time it's called.
  void paintEvent(QPaintEvent *event);

  // Another over-ridden QWidget method.  By default, double-clicking on our
  // rendered area zooms in on the location clicked by a factor of 2.
  void mouseDoubleClickEvent(QMouseEvent *event);

  // These methods handle both the double click zoom and the zoom functionality
  // allowed by the zoom slot.
  void zoomIn(double new_longitude_center, double new_latitude_center,
	      double zoom_factor = 2.0);
  void zoomOut(double new_longitude_center, double new_latitude_center,
	       double zoom_factor = 2.0);

  // When we do need to change the Map/Points/Grid QPixmaps, these methods
  // are called.
  bool paintMap(QPaintDevice *device);
  bool paintPoints(QPaintDevice *device);
  bool paintGrid(QPaintDevice *device);

  // The remainder of these methods are for purely internal usage.  First up,
  // the methods for drawing the Grid QPixmap.
  void _drawLonTick(QPainter* painter, double longitude, bool major_tick);
  void _drawLonTickLabel(QPainter* painter, double longitude);
  void _drawLatTick(QPainter* painter, double latitude, bool major_tick);
  void _drawLatTickLabel(QPainter* painter, double latitude);
  void _drawLonLabel(QPainter* painter, QString& label);
  void _drawLatLabel(QPainter* painter, QString& label);
  void _drawLonGrid(QPainter* painter, double longitude);
  void _drawLatGrid(QPainter* painter, double latitude);

  // A set of methods for returning the weight and coordinate bounds based on
  // either an input Map or set of AngualarCoordinates.
  void _findNewWeightBounds(Stomp::Map* stomp_map);
  void _findNewWeightBounds(Stomp::WAngularVector& ang);
  void _findNewImageBounds(Stomp::Map* stomp_map);
  void _findNewImageBounds(Stomp::WAngularVector& ang);
  void _findNewMaxResolution(Stomp::Map* stomp_map);

  // For _checkBounds, we have no initial idea of what the coordinate bounds
  // are or should be.  Hence, it will handle cases under the assumption of
  // maximum ignorance.  For _enforceBounds, we know we have lon-lat values that
  // are valid, we just need to check against possible double precision errors
  // on the edges of our spherical bounds.
  void _checkBounds();
  void _enforceBounds(double& lon, double& lat);

  // When rendering the pixels, we have a few cases to consider.  If a pixel is
  // so small that rendering it would cover just a few pixels in our QPixmap,
  // we're better off just putting a point there than spending time drawing the
  // full polygon.  If the pixel is big enough to be rendered, but wraps past
  // the longitude disontinuity, then we have to draw both halves separately.
  // Finally, there's just the simple case of taking a pixel and turning it into
  // a polygon.
  void _pixelToPoint(uint32_t pixel_x, uint32_t pixel_y,
		     uint32_t resolution, QPointF& point);
  void _pixelToPolygon(uint32_t pixel_x, uint32_t pixel_y,
		       uint32_t resolution, QPolygonF& polygon);
  void _splitPixelToPolygons(uint32_t pixel_x, uint32_t pixel_y,
			     uint32_t resolution, QPolygonF& left_polygon,
			     QPolygonF& right_polygon);

  // AngularCoordinates are simply turned into QPointFs.
  bool _angToPoint(Stomp::WeightedAngularCoordinate& w_ang, QPointF& point);

  // Four methods for going back and forth between angular space and QPixmap
  // XY pixel-space.
  qreal _lonToX(double longitude);
  qreal _latToY(double latitude);
  double _xToLon(qreal pixel_x);
  double _yToLat(qreal pixel_y);

  // Thanks to the longitude discontinuity and the fact that the longitude spans
  // for our various coordinate systems aren't the same, it's easier to put
  // these values into functions rather than re-calculating them in various
  // methods.
  double _lonRange();
  double _latRange();
  double _lonCenter();
  double _latCenter();

  // Translate between Cartesian and Aitoff-Hammer projections.
  void _cartesianToAitoff(double& longitude, double& latitude);

  // Given an input weight, return a value between [0,1) that we can use to
  // select the appropriate Pallete color.
  double _normalizeWeight(double weight);
  double _normalizePointsWeight(double weight);

  // Based on the current display bounds, figure out the area of a QPixmap
  // pixel.  Any Pixels with areas smaller than that should be rendered with
  // _pixelToPoint instead of _pixelToPolygon
  double _renderPixelArea();

  // Based on the current display bounds, find the version of the input map
  // that contains a reasonable number of Pixels for rendering.
  void _findRenderLevel(uint32_t target_pixels = 100000);

  // Since our set of Maps is complicated (and may contain duplicate pointers),
  // we need to be careful when deleting them so we don't free a pointer twice.
  void _clearMaps();

 private:
  bool antialiased_;
  Palette palette_, points_palette_;
  QPixmap *pixmap_;
  QPixmap *points_pixmap_;
  QPixmap *grid_pixmap_;
  QColor background_color_;
  std::map<uint8_t, Stomp::Map*> stomp_map_;
  bool good_map_;
  Stomp::WAngularVector ang_;
  bool good_points_;
  uint16_t buffer_top_, buffer_bottom_, buffer_left_, buffer_right_;
  bool fill_;
  double lonmin_, lonmax_, latmin_, latmax_;
  bool full_sky_;
  Stomp::AngularCoordinate::Sphere sphere_;
  bool continuous_longitude_;
  bool aitoff_projection_;
  double weight_min_, weight_max_;
  double points_weight_min_, points_weight_max_;
  uint32_t max_resolution_;
  uint8_t render_level_;
  bool auto_update_;
  double *tick_divisions_;
  uint16_t major_tick_length_, minor_tick_length_;
  uint16_t tick_precision_, minor_ticks_per_major_, n_tick_divisions_;
  bool show_coordinates_;
  bool show_grid_;
  bool show_points_, filter_points_;
};

 #endif
