// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the RenderArea class.  This QWidget derivative acts
// as the canvas for displaying all of the Stomp::Map and Stomp::WAngularVector
// data.  It also handles communication between its element objects and the
// StompViewer class.

#ifndef RENDERAREA_H
#define RENDERAREA_H

#include <QBrush>
#include <QPen>
#include <QPixmap>
#include <QWidget>
#include <QPolygonF>
#include <QMouseEvent>
#include <stomp.h>
#include "palette.h"
#include "render_geometry.h"
#include "render_thread.h"
#include "reader_thread.h"

class Palette;
class RenderGeometry;
class RenderMapThread;
class ReadMapThread;
class RenderArea;

class RenderArea : public QWidget {
  Q_OBJECT

 public:
  RenderArea(QWidget *parent = 0);

  QSize minimumSizeHint() const;
  QSize sizeHint() const;
  double mapWeightMin();
  double mapWeightMax();
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
  double mouseLongitude();
  double mouseLatitude();

 public slots:
  // We're storing the map image in a QPixmap so that paint events will be
  // more graceful.  However, if we do something that should change the map,
  // we need to update the underlying map image.
  void updatePixmap(bool send_update_call = true);

  // Our new map image is created via a separate thread call.  This slot is
  // used by the thread to signal that the new image is ready to be displayed.
  void newMapImage(const QImage& new_map_image);

  // Same pair of methods to process changes to the points pixmap.
  void updatePoints(bool send_update_call = true);
  void newPointsImage(const QImage& new_points_image);

  // Likewise, if we alter the image, then we need to re-draw the QPixmap that
  // contains the coordinate grid and the one that contains any points that have
  // been loaded.
  void updateGrid(bool send_update_call = true);

  // Read in a new map.  This will launch the Map reader thread, contingent on
  // us verifying that the named file exists.  If so, we launch the thread and
  // return true.
  bool readNewMap(QString& file_name);

  // The Map reader thread reads in the file indicated by the argument to
  // readMapFile.  After doing so, it emits a signal indicating that the base
  // resolution Map is ready.  That signal is tied to this slot, which stores
  // the pointer, finds the Map bounds, updates the status and
  // starts the Map rendering thread to show the base Map.
  void newBaseMap(Stomp::Map* base_map);

  // After reading in the Map, the reader thread starts constructing the
  // softened Maps.  This takes some time (hence the separate thread), but
  // will make for quicker rendering once they're available.  This slot is tied
  // to this process.
  void newSoftenedMap(uint8_t level, Stomp::Map* softened_map);

  // While the softening is happening, the reader thread will alert us of the
  // progress via this slot.
  void readerProgress(int progress);

  // Read in a new set of points.  Return true if the operation was successful.
  bool readNewPointsFile(QString& file_name,
			 Stomp::AngularCoordinate::Sphere sphere,
			 bool weighted_points);

  // Remove either the current Map or set of points.
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
  void setMapWeightRange(double weight_min, double weight_max);
  void setPointsWeightRange(double weight_min, double weight_max);

  // Add a palette to our object to scale between weight and RGB color.
  void setMapPalette(Palette::PaletteType palette_type);
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
  void setFullSky(bool full_sky);
  void setCoordinateSystem(Stomp::AngularCoordinate::Sphere sphere);

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
  void newMousePosition();
  void newStatus(const QString& message, int timeout);
  void progressUpdate(int new_progress);

 protected:
  // paintEvent is an over-ridden method from the QWidget class that we use to
  // render our displayed Map/Points/Grid image.  paintEvent is called every
  // time the window system needs to re-draw the window, so it's called very
  // frequently.  To keep from having to re-render each time paintEvent is
  // called, we put each of the Map/Points/Grid images into a QPixmap which
  // paintEvent re-draws each time it's called.
  void paintEvent(QPaintEvent *event);

  // We enable mouse tracking in the constructor.  If the mouse is in the
  // RenderArea widget's area, then we convert the mouse position into
  // coordinates based on the Sphere we're plotting.
  void mouseMoveEvent(QMouseEvent *event);

  // Another over-ridden QWidget method.  By default, double-clicking on our
  // rendered area zooms in on the location clicked by a factor of 2.
  void mouseDoubleClickEvent(QMouseEvent *event);

  // We need to keep our RenderGeometry object sync'd up with the size of our
  // RenderArea object, so we need to catch this event.
  void resizeEvent(QResizeEvent *event);

  // These methods handle both the double click zoom and the zoom functionality
  // allowed by the zoom slot.
  void zoomIn(double new_longitude_center, double new_latitude_center,
	      double zoom_factor = 2.0);
  void zoomOut(double new_longitude_center, double new_latitude_center,
	       double zoom_factor = 2.0);

  // When we do need to change the Map/Points/Grid QPixmaps, these methods
  // are called.
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

  // Based on the current display bounds, find the version of the input map
  // that contains a reasonable number of Pixels for rendering.
  void _findRenderLevel(uint32_t target_pixels = 100000);

 private:
  bool antialiased_;
  Palette map_palette_, points_palette_;
  RenderGeometry geom_;
  RenderMapThread render_map_thread_;
  RenderPointsThread render_points_thread_;
  ReadMapThread read_map_thread_;
  QPixmap map_pixmap_;
  QPixmap points_pixmap_;
  QPixmap *grid_pixmap_;
  QColor background_color_;
  std::map<uint8_t, Stomp::Map*> stomp_map_;
  Stomp::Map* base_map_;
  bool good_map_;
  bool good_soften_;
  Stomp::WAngularVector ang_;
  bool good_points_;
  bool fill_;
  double mouse_lon_, mouse_lat_;
  uint32_t max_resolution_;
  uint8_t render_level_;
  bool auto_update_;
  double *tick_divisions_;
  uint32_t major_tick_length_, minor_tick_length_;
  uint32_t tick_precision_, minor_ticks_per_major_, n_tick_divisions_;
  bool show_coordinates_;
  bool show_grid_;
  bool show_points_, filter_points_;
};

#endif
