/*
 * coverer.cc
 *
 *  Created on: Jul 25, 2012
 *      Author: cbmorrison
 */


#include "coverer.h"

namespace s2omp {

coverer::coverer() {
  coverer(0, MAX_LEVEL);
}

coverer::coverer(int min_level, int max_level) {
  min_level_ = min_level;
  max_level_ = max_level;
  max_pixels_ = 8;
}

uint32_t coverer::get_covering(const bound_interface& bound,
    pixel_vector* pixels) {
  return covering(8, bound, pixels);
}

uint32_t coverer::get_covering(uint32_t max_pixels, const bound_interface& bound,
    pixel_vector* pixels) {
  // We first clear our input pixel vector and our priority queue.
  if (!pixels->empty()) pixels->clear();
  if (!pix_q_.empty()) pix_q_.clear();

  // Next we grab a list of initial candidates that cover the circle bound for
  // our target bound.
  pixel_vector initial_candidates;
  get_initial_covering(bound.get_bound(), &initial_candidates);

  // Need something here to test if initial_candidates is empty and if so start
  // with the face pixels and test may_intersect

  // If we have already met or exceeded the number of pixels with the initial
  // candidates then there is no need to continue
  if (initial_candidates.size() >= max_pixels) {
    while (!initial_candidates.empty()) {
      pixels->push_back(initial_candidates.pop_back());
    }
  }

  // We score these initial pixels and store them in the the priority queue
  // based on score from low to high (low being the best score).
  for (pixel_iterator iter = initial_candidates.begin();
      iter != initial_candidates.end(); ++iter) {
    pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
  }


  // While we still have candidates left in the queue, iterate through the queue
  // and test those pixels.
  while (!pix_q_.empty()) {
    // Pop the current highest priority pixel out of the queue.
    pixel candidate = pix_q_.top().second;
    pix_q_.pop();
    // If the pixel we are considering is at lower level than the min_level
    // specified. Break up the pixel into it's children, add them to the queue
    // and continue.
    // Note this could cause more pixels to be returned than requested.
    if (candidate.level() < min_level_) {
      for (pixel_iterator iter = candidate.child_begin(min_level_);
          iter != candidate.child_end(min_level_); ++iter) {
        if (bound.may_intersect(*iter)) {
          pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
        }
      }
      continue;
    }

    // We test to see if the number of pixels left in the queue plus those we
    // have already stored is enough meet the max number of pixels requested.
    // If it is we flush the pixel queue and store the pixels.
    if (pix_q_.size() + pixels->size() + 1 >= max_pixels) {
      pixels->push_back(candidate);
      flush_queue(pixels);
    } else {
      // We check to see if this pixel is already at the maximum level or if it
      // is fully contained within the bound. In either case we add it to the
      // output.
      if (candidate.level() + 1 > max_level_ || bound.contains(candidate)) {
        pixels->push_back(candidate);
      } else {
        // If we can't add the pixel outright we get all of this pixels children
        // and test for may_intersect.
        pixel_vector child_candidates;
        for (pixel_iterator iter = candidate.child_begin();
            iter != candidate.child_end(); ++iter) {
          if (bound.may_intersect(*iter)) {
            child_candidates.push_back(*iter);
          }
        }
        // If the number of intersecting children in this pixel, plus the those
        // already stored, plus those left in the queue exceed the number of
        // pixels requested, add this candidate pixel to the output.
        if (child_candidates.size() + pix_q_.size() + pixels->size() >
          max_pixels) {
          pixels->push_back(candidate);
        } else {
          // If we won't exceed max_pixels with these children then we score
          // them and add them to the queue.
          for (pixel_iterator iter = child_candidates.begin();
              iter != child_candidates.end(); ++iter) {
            pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
          }
        }
      }
    }
  }
  return pixels->size();
}

bool coverer::get_covering(double fractional_area_tolerace,
  const bound_interface& bound, pixel_vector* pixels) {
  // We first clear our input pixel vector and our priority queue.
  if (!pixels->empty()) pixels->clear();
  if (!pix_q_.empty()) pix_q_.clear();
  // In order to compute the current fractional area we need to store both the
  // area of the bound as well as a running tally of both the pixels already
  // stored in the output and those within the queue.
  double bound_area = bound.area();
  double covered_area = 0.0;

  // If the user has specified a fractional tolerance that is not achievable
  // given the mass pixels we exit without returning a covering.
  if (pixel.get_level_from_area(fractional_area_tolerace*bound_area) >
        max_level_) {
    std::cout << "coverer::covering - Fractional area tolerance not "
              << "achievable with specified max_level.\n\tExiting";
    exit(2);
  }

  // Next we grab a list of initial candidates that cover the circle bound for
  // our target bound.
  pixel_vector initial_candidates;
  get_initial_covering(bound.get_bound(), &initial_candidates);

  // Need something here to test if initial_candidates is empty and if so start
  // with the face pixels and test may_intersect

  // We start by adding up the total area covered so far given the initial list
  // of candidates.
  for (pixel_iterator iter = initial_candidates.begin();
      initial_candidates.end(); ++iter) {
    covered_area += iter->exact_area();
  }
  // If we have already met or exceeded the tolerance with the initial
  // candidates then there is no need to continue. We simply push those to the
  // output. Note that for a covering of this type (as in not interior) the
  // covered area will always exceed the bound area.
  if ((covered_area - bound_area)/bound_area <=
      fractional_area_tolerace) {
    while (!initial_candidates.empty()) {
      pixels->push_back(initial_candidates.pop_back());
    }
    // The initial candidate list is now empty. So the for loop and while loops
    // that are next will not run.
  }

  // We score these initial pixels and store them in the the priority queue
  // based on score from low to high (low being the best score).
  for (pixel_iterator iter = initial_candidates.begin();
      iter != initial_candidates.end(); ++iter) {
    pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
  }

  // While we still have candidates left in the queue, iterate through the queue
  // and test those pixels.
  while (!pix_q_.empty()) {
    // If the pixel we are considering is at lower level than the min_level
    // specified. Break up the pixel into it's children, add them to the queue
    // and continue.
    // Note this could cause the actual tolerance to exceed the fractional
    // tolerance.
    if (pix_q_.top().second.level() < min_level_) {
      pixel candidate = pix_q_.top().second;
      pix_q_.pop();
      covered_area -= candidate.exact_area();
      for (pixel_iterator iter = candidate.child_begin(min_level_);
          iter != candidate.child_end(min_level_); ++iter) {
        if (bound.may_intersect(*iter)) {
          pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
          covered_area += iter->exact_area();
        }
      }
      continue;
    }

    // If we have reached the desired tolerance there is no need to continue
    // and hence we flush the pixels from the queue and add them to the output.
    // Note that for this covering the covered_area will always be greater than
    // the bound area.
    if ((covered_area - bound_area)/bound_area <=
        fractional_area_tolerace) {
      flush_queue(pixels);
    } else {
      // Now we pop the highest priority candidate from the queue.
      pixel candidate = pix_q_.top().second;
      pix_q_.pop();
      // Since we've taken the pixel candidate out of the queue, we need to
      // remove it's area from the the running tally. If the pixel is at the max
      // level or it is fully contained within the bound we add the area back
      // in. If these conditions are not met then we will add back in the
      // the area of the children of this candidate that intersect the area.
      double candidate_area = candidate.exact_area();
      covered_area -= candidate_area();
      // If our candidate is at the max level or fully contained by the bound
      // we add it to the output and add it's area back to the running tally.
      if (candidate.level() + 1 > max_level_ || bound.contains(candidate)) {
        pixels->push_back(candidate);
        covered_area += candidate_area;
      } else {
        // If we can't add the pixel outright we get all of this pixels children
        // and test for may_intersect. If the child intersects we score it and
        // add it to the priority queue. We also add it's area into the covered
        // area tally.
        for (pixel_iterator iter = candidate.child_begin();
            iter != candidate.child_end(); ++iter) {
          if (bound.may_intersect(*iter)) {
            pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
            covered_area += iter->exact_area();
          }
        }
      }
    }
  }
  // Since we don't have a goal number of pixels to reach we instead return a
  // boolean specifying if the desired tolerance was reached or not.
  return (covered_area - bound_area)/bound_area <
      fractional_area_tolerace;
}

uint32_t coverer::get_interior_covering(const bound_interface& bound,
    pixel_vector* pixels) {
  // We first clear our input pixel vector and our priority queue.
  if (!pixels->empty()) pixels->clear();
  if (!pix_q_.empty()) pix_q_.clear();

  // Next we grab a list of initial candidates that cover the circle bound for
  // our target bound.
  pixel_vector initial_candidates;
  get_initial_covering(bound.get_bound(), &initial_candidates);

  // Need something here to test if initial_candidates is empty and if so start
  // with the face pixels and test may_intersect

  // We score these initial pixels and store them in the the priority queue
  // based on score from low to high (low being the best score).
  for (pixel_iterator iter = initial_candidates.begin();
      iter != initial_candidates.end(); ++iter) {
    pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
  }

  // While we still have candidates left in the queue, iterate through the queue
  // and test those pixels.
  while (!pix_q_.empty()) {
    pixel candidate = pix_q_.top().second;
    pix_q_.pop();
    // If the pixel we are considering is at lower level than the min_level
    // specified. Break up the pixel into it's children, add them to the queue
    // and continue. Since we are just adding it the children into the queue, we
    // don't need to test for contains just yet.
    if (candidate.level() < min_level_) {
      for (pixel_iterator iter = candidate.child_begin(min_level_);
          iter != candidate.child_end(min_level_); ++iter) {
        if (bound.may_intersect(*iter)) {
          pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
        }
      }
      continue;
    }
    // For an interior covering, we only add pixels to the output if they are
    // fully contained within the bound.
    if (bound.contains(candidate)) {
      pixels->push_back(candidate);
    } else if (candidate.level() < max_level_) {
      // If the candidate pixel is not contained and is less than the max level
      // we immediately refine it into it's children and add those children that
      // may intersect into the queue.
      for (pixel_iterator iter = candidate.child_begin();
          iter != candidate.child_end(); ++iter) {
        if (bound.may_intersect(*iter)) {
          pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
        }
      }
    }
  }
  // Return the number of pixels used in the covering.
  return pixels->size();
}

bool coverer::get_interior_covering(double fractional_area_tolerace,
    const bound_interface& bound, pixel_vector* pixels) {
  // We first clear our input pixel vector and our priority queue.
  if (!pixels->empty()) pixels->clear();
  if (!pix_q_.empty()) pix_q_.clear();
  double bound_area = bound.area();
  double covered_area = 0.0;

  if (pixel.get_level_from_area(fractional_area_tolerace*bound_area) >
      max_level_) {
    std::cout << "coverer::interior_covering - Fractional area tolerance not "
              << "achievable with specified max_level.\n\tExiting";
    exit(2);
  }

  // Next we grab a list of initial candidates that cover the circle bound for
  // our target bound.
  pixel_vector initial_candidates;
  get_initial_covering(bound.get_bound(), &initial_candidates);

  // We score these initial pixels and store them in the the priority queue
  // based on score from low to high (low being the best score).
  for (pixel_iterator iter = initial_candidates.begin();
      iter != initial_candidates.end(); ++iter) {
    pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
  }

  // While we still have candidates left in the queue and the fractional
  // tolerance has not been met, iterate through the queue and test the pixels
  // in the queue. Note that for an interior covering the covered area will
  // always be less than the true area of the bound.
  while (!pix_q_.empty() && (bound_area - covered_area)/bound_area >
      fractional_area_tolerace) {
    pixel candidate = pix_q_.top().second;
    pix_q_.pop();
    // If the pixel we are considering is at lower level than the min_level
    // specified. Break up the pixel into it's children, add them to the queue
    // and continue. Since we are just adding it the children into the queue, we
    // don't need to test for contains just yet.
    if (candidate.level() < min_level_) {
      for (pixel_iterator iter = candidate.child_begin(min_level_);
          iter != candidate.child_end(min_level_); ++iter) {
        if (bound.may_intersect(*iter)) {
          pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
        }
      }
      continue;
    }
    // If the candidate pixel is fully contained within the bound we add it to
    // the output and add the area to the running tally.
    if (bound.contains(candidate)) {
      pixels->push_back(candidate);
      covered_area += candidate->exact_area();
    } else if (candidate.level() < max_level_) {
      // If we can't add the pixel outright we get all of this pixel's children
      // and test for may_intersect. If the child intersects the bound we score
      // it and add it to the queue.
      for (pixel_iterator iter = candidate.child_begin();
          iter != candidate.child_end(); ++iter) {
        if (bound.may_intersect(*iter)) {
          pix_q_.push(pixel_pair(score_pixel(bound, *iter), *iter));
        }
      }
    }
  }
  // This is the only method where we are not guaranteed to clear the pixel
  // queue before exiting the while loop. As such, if we want to free up some
  // memory and keep this coverer around we need to clear the queue by hand.
  if (!pix_q_.empty()) pix_q_.clear();

  // Again, instead of returning the number of pixels we instead return a
  // boolean specifying if the fractional tolerance was met before reaching max
  // level.
  return (bound_area - covered_area)/bound_area >
    fractional_area_tolerace;
}

void coverer::get_simple_covering(const bound_interface& bound, int level,
    pixel_vector* pixels) {
  // Clear the input pixel vector.
  if (!pixels->empty()) pixels->clear();

  // The simple covering starts with a single point which we convert to a pixel
  // and test it as well as it's neighbors (and neighbors neighbors) until the
  // bound is covered. The starting pixel we choose is the central point of the
  // circle bound that covers the input bound.
  pixel start = bound.get_bound().axis().to_pixel(level);
  while (!bound.may_intersect(start)) {
    start = start.next_wrap();
  }
  pixel_map kept_map;
  pixel_vector canidates;
  canidates.push_back(start);

  while (!canidates.empty()) {
    pixel pix = (*canidates.back());
    canidates.pop_back();
    if (!bound.may_itersect(pix)) continue;
    pixels->push_back(pix);

    pixel_vector neighbors;
    pix.neighbors(&neighbors);
    kept_map.insert(std::pair<uint64, pixel>(pix.id(), pix));
    for (pixel_iterator iter = neighbors.begin();
        iter != neighbors.end(); ++iter) {
      if (kept_map.find(iter->id())->empty()) {
        canidates.push_back(*iter);
      }
    }
  }
}

void coverer::get_initial_covering(const circle_bound& bound,
    pixel_vector* pixels) {
  if (!pixels.empty()) pixels.clear();

  int level = pixel.get_level_from_area(bound.area());
  if (level < 0) {
    level = 0;
  } else if (level > max_level_) {
    level = max_level_;
  }

  get_simple_covering(bound, level, pixels);
}

int coverer::score_pixel(const bound_interface& bound, const pixel& pix) {
  // We want to sort the pixels in our priority queue first by their size,
  // next by the number of children that may_intersect the bound, and then
  // by the number of children that are terminal (for non-interior coverings
  /// this means child.level() == max_level or bound.contains(child))
  int n_children, n_terminals = 0;
  for (pixel_iterator iter = pix.child_begin();
      iter != pix.child_end(); ++iter) {
    if (bound.may_intersect(*iter)) {
      n_children++;
      if (iter->level() + 1 > max_level_ || bound.contains(*iter)) {
        n_terminals++;
      }
    }
  }
  // This is the value that is used within the priority queue. First we
  // compute the level of the pixel and shift by 2 bits (since there at most 4
  // children) and store the children. We then shift 2 bits and store the
  // number of terminal children. The negative sign flips all the bits so that
  // a pixel with a low level, most children, and most terminals will be
  // checked first in the queue. The pixel with the highest possible priority
  // would be one at level 0 with only one child that may_intersect the bound.
  return -(((pix.level() << 2) + n_children << 2) + n_terminals);
}

void coverer::flush_queue(pixel_vector* pixels) {
  while (!pix_q_.empty()) {
    pixels->push_back(pix_q_.top().second);
    pix_q_.pop();
  }
}

} // end namespace s2omp

