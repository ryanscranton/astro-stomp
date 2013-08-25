/*
 * coverer.cc
 *
 *  Created on: Jul 25, 2012
 *      Author: cbmorrison
 */

#include <s2/s2.h>
#if defined(OS_MACOSX)
#include <ext/hash_set>
#else
#include <hash_set>
#endif
using __gnu_cxx::hash_set;

#include "bound_interface.h"
#include "coverer.h"
#include "pixel.h"
#include "point.h"

namespace s2omp {

coverer::coverer() {
  min_level_ = 0;
  max_level_ = MAX_LEVEL;
}

coverer::coverer(int min_level, int max_level) {
  min_level_ = min_level;
  max_level_ = max_level;
}

coverer::~coverer() {
}

bool coverer::get_covering(const bound_interface& bound, pixel_vector* pixels) {
  return generate_covering(bound, DEFAULT_COVERING_PIXELS, false, -1.0, pixels);
}

bool coverer::get_covering(
    long max_pixels, const bound_interface& bound, pixel_vector* pixels) {
  return generate_covering(bound, max_pixels, false, -1.0, pixels);
}

bool coverer::get_covering(
    double fractional_area_tolerace, const bound_interface& bound,
    pixel_vector* pixels) {
  return generate_covering(bound, 0, false, fractional_area_tolerace,
      pixels);
}

bool coverer::get_interior_covering(
    const bound_interface& bound, pixel_vector* pixels) {
  return generate_covering(bound, 0, true, -1.0, pixels);
}

bool coverer::get_interior_covering(
    double fractional_area_tolerace, const bound_interface& bound,
    pixel_vector* pixels) {
  return generate_covering(bound, 0, true, fractional_area_tolerace, pixels);
}

bool coverer::generate_covering(
    const bound_interface& bound, long max_pixels, bool interior,
    double fraction, pixel_vector* pixels) {
  // We first clear our input pixel vector and our priority queue.
  if (!pixels->empty()) {
    pixels->clear();
  }
  while (!pix_q_.empty()) {
    pix_q_.pop();
  }

  // Define a few convenience variables for storing the bound area and the
  // current area of the stored pixels. We will use these if fractional is
  // requested.
  double bound_area = bound.area();
  double covered_area = 0.0;

  // Next we grab a list of initial candidates that cover the circle bound for
  // our target bound.
  pixel_vector initial_candidates;
  get_initial_covering(bound, &initial_candidates);

  // Need something here to test if initial_candidates is empty and if so start
  // with the face pixels and test may_intersect

  // We score these initial pixels and store them in the the priority queue
  // based on score from low to high (low being the best score).
  for (pixel_iterator iter = initial_candidates.begin(); iter
      != initial_candidates.end(); ++iter) {
    new_candidate(bound, *iter);
    if (!interior && fraction > 0.0)
      covered_area += iter->exact_area();
  }

  // If we have already reached our convergence criteria for a non-interior
  // covering, we flush the queue and output the initial candidates. We
  // test if the max pixels have been reached or if a fractional area has been
  // requested, we test if this tolerance has been reached.
  if (!interior) {
    if (fraction < 0.0 && initial_candidates.size() >= max_pixels) {
      flush_queue(pixels);
    } else if (fraction >= 0.0 &&
        (covered_area - bound_area)/bound_area <= fraction) {
      flush_queue(pixels);
    }
  }

  // While we still have candidates left in the queue, iterate through the queue
  // and test those pixels.
  while (!pix_q_.empty()) {

    // If we have an interior covering and we have already reached the precision
    // set in fraction, we quit and clear the queue.
    if (interior && fraction > 0.0 &&
        (bound_area - covered_area)/bound_area < fraction) {
      while (!pix_q_.empty()) {
        pix_q_.pop();
      }
      break;
    }

    // Pop the current highest priority pixel out of the queue and if our
    // covering is not interior we subtract it's area from the total in case we
    // refine the pixel later.
    pixel_candidate candidate = pix_q_.top().second;
    pix_q_.pop();
    if (!interior && fraction > 0.0)
      covered_area -= candidate.pix.exact_area();

    // If the pixel we are considering is at lower level than the min_level
    // specified. Break up the pixel into it's children, add them to the queue,
    // and continue.
    // Note this could cause more pixels to be returned or a smaller fraction
    // than initially requested if those are the modes.
    if (candidate.pix.level() < min_level_) {
      for (pixel child = candidate.pix.child_begin(min_level_); child
               != candidate.pix.child_end(min_level_); child = child.next()) {
        if (bound.may_intersect(child)) {
          new_candidate(bound, child);
          if (!interior && fraction > 0.0)
            covered_area += child.exact_area();
        }
      }
      continue;
    }

    // The cases for non-interior coverings are more complicated, though less
    // computationally intensive, hence we tackle those first.
    if (!interior) {

      // We first see if with the current pixels both in the output and in the
      // queue we have reached our goal tolerance in max_pixels or area. If
      // we have we add the current and all remaining candidates to the output,
      // clearing the queue.
      if ((fraction <= 0.0 &&
            pix_q_.size() + pixels->size() + 1 >= max_pixels) ||
          (fraction > 0.0 &&
            (covered_area + candidate.pix.exact_area() - bound_area)/
            bound_area < fraction)) {
        pixels->push_back(candidate.pix);
        flush_queue(pixels);
        break;
      }

      // If we haven't reached our goal yet we test if this pixel can be added.
      // We test is_terminal which says that all of the child pixels of this
      // candidate are intersecting and are either contained or at max_level_ or
      // the candidate is at max_level_. For a non-interior covering we then
      // add this pixel.
      if (candidate.is_terminal) {
        pixels->push_back(candidate.pix);
        if (fraction > 0.0)
          covered_area += candidate.pix.exact_area();
      } else {

        // Since we can't add this pixel outright, we first test to see, for a
        // non-fractional covering, if adding this pixel's children to the queue
        // will exceed our max points. If it does we add this candidate, if it
        // doesn't we resolve it to the next level.
        if (fraction <= 0.0 && pixels->size() + pix_q_.size() +
            candidate.n_children >= max_pixels) {
          pixels->push_back(candidate.pix);
        } else {
          for (pixel child = candidate.pix.child_begin();
               child != candidate.pix.child_end(); child = child.next()) {
            if (bound.may_intersect(child)) {
              new_candidate(bound, child);
              if (fraction > 0.0)
                covered_area += child.exact_area();
            }
          }
        }
      }
    } else {

      //The interior covering is a bit simpler if more expensive
      // computationally. We first test if the pixel is fully contained. If it
      // is we add it to the output and, if a fractional area is specified, we
      // add this pixel's area to the total.
      if (bound.contains(candidate.pix)) {
        pixels->push_back(candidate.pix);
        if (fraction > 0.0)
          covered_area += candidate.pix.exact_area();
      } else {

        // If we can't add this pixel, we resolve it's children and add them to
        // the queue.
        for (pixel child = candidate.pix.child_begin(); child
                 != candidate.pix.child_end(); child = child.next()) {
          if (bound.may_intersect(child)) {
            new_candidate(bound, child);
          }
        }
      }
    }
  }

  // We define "success" if a) we return fewer than or equal to the number of
  // pixels requested b) we have reached the fractional area tolerance or
  // c), for interior coverings only, we return a non-empty vector.
  bool success = false;
  if (!interior && fraction <= 0.0) {
    if (!pixels->empty() && pixels->size() <= max_pixels)
      success = true;
  } else if (fraction > 0.0) {
    if (!pixels->empty() &&
        fabs(bound_area - covered_area)/bound_area < fraction)
      success = true;
  } else {
    success = !pixels->empty();
  }
  return success;
}

void coverer::get_simple_covering(
    const bound_interface& bound, int level, pixel_vector* pixels) {
  // Clear the input pixel vector.
  if (!pixels->empty()) pixels->clear();

  // The simple covering starts with a single point which we convert to a pixel
  // and test it as well as it's neighbors (and neighbors neighbors) until the
  // bound is covered. The starting pixel we choose is the central point of the
  // circle bound that covers the input bound.
  pixel starting_pixel = bound.get_center().to_pixel(level);
  while (!bound.may_intersect(starting_pixel)) {
    starting_pixel = starting_pixel.next_wrap();
  }

  pixel_vector candidates;
  candidates.push_back(starting_pixel);

  hash_set<uint64> candidate_ids;
  candidate_ids.insert(starting_pixel.id());

  while (!candidates.empty()) {
    pixel pix = candidates.back();
    candidates.pop_back();
    if (!bound.may_intersect(pix)) continue;
    pixels->push_back(pix);

    pixel_vector neighbors;
    pix.neighbors(&neighbors);
    for (pixel_iterator iter = neighbors.begin();
         iter != neighbors.end(); ++iter) {
      if (candidate_ids.insert(iter->id()).second) {
        candidates.push_back(*iter);
      }
    }
  }

  sort(pixels->begin(), pixels->end());
}

void coverer::get_center_covering(
    const bound_interface& bound, int level, pixel_vector* pixels) {
  // Clear the input pixel vector.
  if (!pixels->empty()) pixels->clear();

  // Follow the same logic here as get_simple_covering, but check that
  // the center of the candidate pixel is contained rather than if the
  // pixel intersects.
  pixel starting_pixel = bound.get_center().to_pixel(level);
  while (!bound.contains(starting_pixel.center_point())) {
    starting_pixel = starting_pixel.next_wrap();
  }

  pixel_vector candidates;
  candidates.push_back(starting_pixel);

  hash_set<uint64> candidate_ids;
  candidate_ids.insert(starting_pixel.id());

  while (!candidates.empty()) {
    pixel pix = candidates.back();
    candidates.pop_back();
    if (!bound.contains(pix.center_point())) continue;

    pixels->push_back(pix);

    pixel_vector neighbors;
    pix.neighbors(&neighbors);
    for (pixel_iterator iter = neighbors.begin();
         iter != neighbors.end(); ++iter) {
      if (candidate_ids.insert(iter->id()).second) {
        candidates.push_back(*iter);
      }
    }
  }

  sort(pixels->begin(), pixels->end());
}

bool coverer::set_min_max_level(int min, int max) {
  if (min <= max && min >= 0 && min <= MAX_LEVEL &&
      max >= 0 && max < MAX_LEVEL) {
    min_level_ = min;
    max_level_ = max;
    return true;
  }
  return false;
}

void coverer::get_initial_covering(
    const bound_interface& bound, pixel_vector* pixels) {
  if (!pixels->empty())
    pixels->clear();

  int level = pixel::get_level_from_area(bound.area());
  if (level < 0) {
    level = 0;
  } else if (level > max_level_) {
    level = max_level_;
  }

  get_simple_covering(bound, level, pixels);
}

void coverer::flush_queue(pixel_vector* pixels) {
  while (!pix_q_.empty()) {
    pixels->push_back(pix_q_.top().second.pix);
    pix_q_.pop();
  }
}

void coverer::new_candidate(const bound_interface& bound,
    const pixel& pix) {

  // Create the pixel_candidate object and initialize default values.
  pixel_candidate pix_cand;
  pix_cand.pix = pix;
  pix_cand.n_children = 0;
  pix_cand.is_terminal = false;

  // We use score pixel in a two fold sense. First it tells us where to place
  // this candidate in the priority queue. Second, it computes both the total
  // number of children that may intersect the bound as well as computing if
  // this pixel is terminal that is pix.level > max_leve_ or pix.contains().
  int score = score_pixel(bound, &pix_cand);

  // Once we have our candidate scored and initialized we add it too the queue.
  pix_q_.push(pixel_entry(score, pix_cand));
}

int coverer::score_pixel(const bound_interface& bound,
    pixel_candidate* pix_cand) {

  // We want to sort the pixels in our priority queue first by their size,
  // next by the number of children that may_intersect the bound, and then
  // by the number of children that are terminal (for non-interior coverings
  /// this means child.level() == max_level or bound.contains(child))
  int n_children = 0, n_terminals = 0;
  for (pixel child = pix_cand->pix.child_begin();
       child != pix_cand->pix.child_end(); child = child.next()) {
    if (bound.may_intersect(child)) {
      n_children++;
      if (child.level() + 1 > max_level_ || bound.contains(child)) {
        n_terminals++;
      }
    }
  }
  pix_cand->n_children = n_children;
  if (pix_cand->pix.level() + 1 > max_level_ || n_terminals == 4)
    pix_cand->is_terminal = true;


  // This is the value that is used within the priority queue. First we
  // compute the level of the pixel and shift by 2 bits (since there at most 4
  // children) and store the children. We then shift 2 bits and store the
  // number of terminal children. The negative sign flips all the bits so that
  // a pixel with a low level, most children, and most terminals will be
  // checked first in the queue. The pixel with the highest possible priority
  // would be one at level 0 with only one child that may_intersect the bound.
  return -(((pix_cand->pix.level() << 2) + n_children << 2) + n_terminals);
}

} // end namespace s2omp

