#include <QtGui>
#include "reader_thread.h"

ReadMapThread::ReadMapThread(QObject *parent) : QThread(parent) {
  restart = false;
  abort = false;
  good_base_map_ = false;
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

void ReadMapThread::run() {
  forever {
    mutex.lock();
    std::string map_file_name = map_file_name_;
    mutex.unlock();

    if (restart) break;
    if (abort) return;
    std::cout << "Reading from " << map_file_name << "...\n";
    Stomp::Map* base_map = new Stomp::Map(map_file_name);
    std::cout << "Done.\n";
    if (restart) break;
    if (abort) return;

    bool good_base_map = !base_map->Empty();

    if (good_base_map) {
      mutex.lock();
      base_map_ = base_map;
      good_base_map_ = good_base_map_;
      mutex.unlock();
      emit newBaseMap(base_map);
    }

    mutex.lock();
    if (!restart) condition.wait(&mutex);
    restart = false;
    mutex.unlock();
  }
}

