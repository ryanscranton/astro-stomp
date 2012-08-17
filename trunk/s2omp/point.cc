/*
 * point.cc
 *
 *  Created on: Aug 16, 2012
 *      Author: morrison
 */




#include "point.h"

pixel* point::to_pixel() const {
  return new pixel(id());
}

pixel* point::to_pixel(int level) const {
  return to_pixel()->parent(level);
}
