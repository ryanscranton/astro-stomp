#include <QtGui>
#include "reader_thread.h"

ReadMapThread::ReadMapThread(QObject *parent) : QThread(parent) {
  restart = false;
  abort = false;
  good_base_map_ = false;
  softened_level_ = Stomp::MaxPixelLevel + 1;
}

ReadMapThread::~ReadMapThread() {
  mutex.lock();
  abort = true;
  condition.wakeOne();
  mutex.unlock();

  wait();
}

void ReadMapThread::readMap(const std::string& map_file_name) {
  QMutexLocker locker(&mutex);

  map_file_name_ = map_file_name;

  if (!isRunning()) {
    start(LowPriority);
  } else {
    restart = true;
    condition.wakeOne();
  }
}

void ReadMapThread::clearMaps() {
  if (softened_level_ < Stomp::MaxPixelLevel) {
    delete softened_maps_[Stomp::MaxPixelLevel];
    for (uint8_t level=Stomp::MaxPixelLevel-1;level>=softened_level_;level--) {
      if (softened_maps_[level] != softened_maps_[level+1])
	delete softened_maps_[level];
    }
  }

  softened_maps_.clear();
  good_base_map_ = false;
  softened_level_ = Stomp::MaxPixelLevel + 1;
}

void ReadMapThread::run() {
  forever {
    mutex.lock();
    std::string map_file_name = map_file_name_;
    clearMaps();
    mutex.unlock();

    std::cout << "Reading from " << map_file_name << "...\n";
    Stomp::Map* base_map = new Stomp::Map(map_file_name);
    std::cout << "Done.  Softening...\n";
    if (restart) break;
    if (abort) return;

    bool good_base_map = !base_map->Empty();

    if (good_base_map) {
      mutex.lock();
      base_map_ = base_map;
      good_base_map_ = good_base_map_;
      mutex.unlock();
      emit newBaseMap(base_map);
      emit readerProgress(1);

      // Now, we iterate over our possible resolution values for rendering
      // the Map.  If the resolution limit is higher than the maximum resolution
      // of the basic Map, then we assign a copy of the pointer to
      // that resolution.
      for (uint8_t level=Stomp::MaxPixelLevel;level>=Stomp::HPixLevel;level--) {
	if (restart) break;
	if (abort) return;
	mutex.lock();
	uint32_t resolution = 1 << level;
	std::cout << "\t" << static_cast<int>(level) << "\n";
	Stomp::Map* soft_map;
	if (resolution >= base_map->MaxResolution()) {
	  soft_map = base_map;
	} else {
	  // For the lower resolution versions of the map, we create softened
	  // versions of the map.  We can do this progressively, so that the
	  // computation is faster than resampling from the basic Map every
	  // time.
	  soft_map = new Stomp::Map();

	  softened_maps_[level+1]->Soften(*soft_map, resolution, true);
	}
	softened_maps_[level] = soft_map;
	softened_level_ = level;
	emit newSoftenedMap(level, soft_map);
	emit readerProgress(static_cast<int>(17 - level));
	mutex.unlock();
      }
      std::cout << "Done\n";
    }

    mutex.lock();
    if (!restart) condition.wait(&mutex);
    restart = false;
    mutex.unlock();
  }
}

