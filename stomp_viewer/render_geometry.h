// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the RenderGeometry class, which handles
// translations between angular coordinates for our Stomp objects and the pixel
// coordinates used in RenderArea.

#ifndef RENDERGEOMETRY_H
#define RENDERGEOMETRY_H

#include <QBrush>
#include <QPen>
#include <QPixmap>
#include <QWidget>
#include <QPolygonF>
#include <QPointF>
#include <stomp.h>

class RenderGeometry {
  // This class handles all of the translation between the spherical
  // coordinates of the Map or WAngularVector and the pixel coordinates of the
  // corresponding QPixmap.

 public:
  RenderGeometry();
  RenderGeometry(double longitude_min, double longitude_max,
                 double latitude_min, double latitude_max,
                 Stomp::AngularCoordinate::Sphere sphere,
                 int width, int height);

  // Simple getters for our basic parameters.
  double longitudeMin();
  double longitudeMax();
  double latitudeMin();
  double latitudeMax();
  Stomp::AngularCoordinate::Sphere sphere();
  int height();
  int width();
  bool aitoffProjection();
  bool fullSky();

  // In case we don't want to have our plotting area stretch all the way
  // to the edge of the image.
  void setVerticalBuffer(uint16_t buffer_pixels);
  void setHorizontalBuffer(uint16_t buffer_pixels);
  void setTopBuffer(uint16_t buffer_pixels);
  void setBottomBuffer(uint16_t buffer_pixels);
  void setLeftBuffer(uint16_t buffer_pixels);
  void setRightBuffer(uint16_t buffer_pixels);

  // We can plot our Maps in any of the available coordinate systems and these
  // methods allow us to choose accordingly.  If the bounds given exceed the
  // natural bounds of the coordinates, then they will be silently fixed.
  void setImageBounds(double lonmin, double lonmax,
                      double latmin, double latmax);
  void setImageBounds(double lonmin, double lonmax,
                      double latmin, double latmax,
                      Stomp::AngularCoordinate::Sphere sphere);
  void setCoordinateSystem(Stomp::AngularCoordinate::Sphere sphere);

  // Setters to keep our class in sync with the size of the RenderArea.
  void setHeight(int height);
  void setWidth(int width);
  void setSize(const QSize& new_size);

  // Change to or from an Aitoff-Hammer projection.
  void useAitoffProjection(bool use_aitoff_projection);

  // Increase current bounds to cover the entire sky for our current projection.
  void setFullSky(bool full_sky);

  // For _checkBounds, we have no initial idea of what the coordinate bounds
  // are or should be.  Hence, it will handle cases under the assumption of
  // maximum ignorance.  For _enforceBounds, we know we have lon-lat values that
  // are valid, we just need to check against possible double precision errors
  // on the edges of our spherical bounds.
  void checkBounds();
  void enforceBounds(double& lon, double& lat);

  // When rendering the pixels, we have a few cases to consider.  If a pixel is
  // so small that rendering it would cover just a few pixels in our QPixmap,
  // we're better off just putting a point there than spending time drawing the
  // full polygon.  If the pixel is big enough to be rendered, but wraps past
  // the longitude disontinuity, then we have to draw both halves separately.
  // Finally, there's just the simple case of taking a pixel and turning it into
  // a polygon.
  void pixelToPoint(uint32_t pixel_x, uint32_t pixel_y,
                    uint32_t resolution, QPointF& point);
  void pixelToPolygon(uint32_t pixel_x, uint32_t pixel_y,
                      uint32_t resolution, uint32_t max_resolution,
                      QPolygonF& polygon);
  void splitPixelToPolygons(uint32_t pixel_x, uint32_t pixel_y,
                            uint32_t resolution, uint32_t max_resolution,
                            QPolygonF& left_polygon, QPolygonF& right_polygon);

  // AngularCoordinates are simply turned into QPointFs.
  bool angToPoint(Stomp::WeightedAngularCoordinate& w_ang, QPointF& point);

  // Four methods for going back and forth between angular space and QPixmap
  // XY pixel-space.
  qreal lonToX(double longitude);
  qreal latToY(double latitude);
  double xToLon(qreal pixel_x);
  double yToLat(qreal pixel_y);

  // Thanks to the longitude discontinuity and the fact that the longitude spans
  // for our various coordinate systems aren't the same, it's easier to put
  // these values into functions rather than re-calculating them in various
  // methods.
  double lonRange();
  double latRange();
  double lonCenter();
  double latCenter();

  // Translate between Cartesian and Aitoff-Hammer projections.
  void cartesianToAitoff(double& longitude, double& latitude);

  // Based on the current display bounds, figure out the area of a QPixmap
  // pixel.  Any Stomp::Pixels with areas smaller than that should be rendered
  // with pixelToPoint instead of pixelToPolygon
  double renderPixelArea();

 private:
  int width_, height_;
  uint16_t buffer_top_, buffer_bottom_, buffer_left_, buffer_right_;
  double lonmin_, lonmax_, latmin_, latmax_;
  Stomp::AngularCoordinate::Sphere sphere_;
  bool continuous_longitude_;
  bool aitoff_projection_;
  bool full_sky_;
};

#endif
