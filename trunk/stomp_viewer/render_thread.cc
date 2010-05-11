#include <QtGui>
#include <QBrush>
#include <QPen>
#include <QWidget>
#include <QPolygonF>
#include "render_area.h"
#include "render_thread.h"

RenderMapThread::RenderMapThread(QObject *parent) : QThread(parent) {
  restart = false;
  abort = false;
}

RenderMapThread::~RenderMapThread() {
  mutex.lock();
  abort = true;
  condition.wakeOne();
  mutex.unlock();

  wait();
}

void RenderMapThread::renderMap(Stomp::Map* stomp_map, RenderGeometry& geometry,
				Palette& palette, uint32_t max_resolution,
				bool antialiased, bool fill) {
  QMutexLocker locker(&mutex);

  stomp_map_ = stomp_map;
  geom_ = geometry;
  palette_ = palette;
  palette_.initialize(palette.currentPaletteType()); // Need this for deep copy
  max_resolution_ = max_resolution;
  antialiased_ = antialiased;
  fill_ = fill;

  if (!isRunning()) {
    start(LowPriority);
  } else {
    restart = true;
    condition.wakeOne();
  }
}

void RenderMapThread::run() {
  forever {
    mutex.lock();
    Stomp::Map* stomp_map = stomp_map_;
    RenderGeometry geom = geom_;
    Palette palette = palette_;
    palette.initialize(palette_.currentPaletteType());
    uint32_t max_resolution = max_resolution_;
    bool antialiased = antialiased_;
    bool fill = fill_;
    mutex.unlock();

    QImage image(geom.width(), geom.height(),
		 QImage::Format_ARGB32_Premultiplied);
    image.fill(QColor(Qt::white).rgb());

    QPainter painter;

    bool accessed_device = painter.begin(&image);

    if (accessed_device) {
      if (antialiased) {
	painter.setRenderHint(QPainter::Antialiasing, true);
	painter.translate(+0.5, +0.5);
      }

      if (!stomp_map->Empty()) {

	double render_pixel_area = geom.renderPixelArea();

	for (Stomp::MapIterator iter=stomp_map->Begin();
	     iter!=stomp_map->End();stomp_map->Iterate(&iter)) {
	  if (restart) break;
	  if (abort) return;
	  Stomp::PixelIterator pix = iter.second;
	  if (pix->Resolution() <= max_resolution &&
	      (geom.fullSky() || pix->IntersectsBounds(geom.longitudeMin(),
						       geom.longitudeMax(),
						       geom.latitudeMin(),
						       geom.latitudeMax(),
						       geom.sphere()))) {
	    QBrush pixel_brush = QBrush(palette.color(pix->Weight()));
	    painter.setBrush(pixel_brush);
	    painter.setPen(QPen(pixel_brush, 0));

	    if (pix->Area() < render_pixel_area) {
	      // If the pixel is so small that displaying it wouldn't cover at
	      // least a single pixel, we just use a single QPointF.
	      QPointF point;
	      geom.pixelToPoint(pix->PixelX(), pix->PixelY(),
				pix->Resolution(), point);
	      painter.drawPoint(point);
	    } else {
	      if (pix->ContinuousBounds(geom.sphere())) {
		QPolygonF polygon;
		geom.pixelToPolygon(pix->PixelX(), pix->PixelY(),
				    pix->Resolution(), max_resolution, polygon);
		if (fill) {
		  painter.drawConvexPolygon(polygon);
		} else {
		  painter.drawPolyline(polygon);
		}
	      } else {
		QPolygonF left_polygon;
		QPolygonF right_polygon;
		geom.splitPixelToPolygons(pix->PixelX(), pix->PixelY(),
					  pix->Resolution(), max_resolution,
					  left_polygon, right_polygon);
		if (fill) {
		  painter.drawConvexPolygon(left_polygon);
		  painter.drawConvexPolygon(right_polygon);
		} else {
		  painter.drawPolyline(left_polygon);
		  painter.drawPolyline(right_polygon);
		}
	      }
	    }
	  }
	}
      }
      accessed_device = painter.end();
    }
    if (!restart && accessed_device) emit renderedImage(image);

    mutex.lock();
    if (!restart) condition.wait(&mutex);
    restart = false;
    mutex.unlock();
  }
}

RenderPointsThread::RenderPointsThread(QObject *parent) : QThread(parent) {
  restart = false;
  abort = false;
}

RenderPointsThread::~RenderPointsThread() {
  mutex.lock();
  abort = true;
  condition.wakeOne();
  mutex.unlock();

  wait();
}

void RenderPointsThread::renderPoints(Stomp::WAngularVector* ang,
				      RenderGeometry& geometry,
				      Palette& palette, bool antialiased) {
  QMutexLocker locker(&mutex);

  ang_ = ang;
  geom_ = geometry;
  palette_ = palette;
  palette_.initialize(palette.currentPaletteType()); // Need this for deep copy
  antialiased_ = antialiased;

  if (!isRunning()) {
    start(LowPriority);
  } else {
    restart = true;
    condition.wakeOne();
  }
}

void RenderPointsThread::run() {
  forever {
    mutex.lock();
    Stomp::WAngularVector* ang = ang_;
    RenderGeometry geom = geom_;
    Palette palette = palette_;
    palette.initialize(palette_.currentPaletteType());
    bool antialiased = antialiased_;
    mutex.unlock();

    QImage image(geom.width(), geom.height(),
		 QImage::Format_ARGB32_Premultiplied);
    image.fill(QColor(Qt::transparent).rgb());

    QPainter painter;

    bool accessed_device = painter.begin(&image);

    if (accessed_device) {
      if (antialiased) {
	painter.setRenderHint(QPainter::Antialiasing, true);
	painter.translate(+0.5, +0.5);
      }

      if (!ang->empty()) {
	for (Stomp::WAngularIterator iter=ang->begin();
	     iter!=ang->end();++iter) {
	  if (restart) break;
	  if (abort) return;
	  QPointF point;
	  if (geom_.angToPoint(*iter, point)) {
	    QBrush pixel_brush = QBrush(palette.color(iter->Weight()));
	    painter.setBrush(pixel_brush);
	    painter.setPen(QPen(pixel_brush, 0));
	    painter.drawPoint(point);
	  }
	}
      }
      accessed_device = painter.end();
      if (!restart && accessed_device) emit renderedImage(image);
    }

    mutex.lock();
    if (!restart) condition.wait(&mutex);
    restart = false;
    mutex.unlock();
  }
}
