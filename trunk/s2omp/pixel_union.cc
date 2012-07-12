/*
 * pixel_union.cc
 *
 *  Created on: Jul 11, 2012
 *      Author: cbmorrison
 */




pixel_union::pixel_union() {
	min_level_ = 64; // don't know if these numbers are right, should look this up
	max_level_ = 0;
	double area_ = 0.0;
	bool initialized_ = false;
}

void pixel_union::init(const pixel_vector& pixels) {
	// uuhhh I'm going to look harder at Stomp and S2 before doing this
	// can't think of a good way to do this right now that isn't slow.

	// okay, here's an idea. take the pixels and separate them into vectors by
	// level, ordered by id. If we see a pixel that is first in it's parent, we
	// look to see if the next 3 pixels form the parent. If they do we create the
	// parent and insert it into the vector for the above level. In this we can
	// create new vector levels if we can create parents at the highest level.
	// we then report back the concatenated pixel vector order by low level(
	// large pixel) first and then id.

  if (!pixels_.empty()) pixels_.clear();
  pixels_.reserve(pixels.size());

  max_level_ = -1;
  min_level_ = 64;
  while (!pixels.empty()) {
    int level = 64;
    uint64 id = 2^64;
    unit32_t idx -1;
    unit32_t idx_run = 0;

    for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
      if (iter->level() <= level) {
        if (iter->id() < id) {
          level = iter->level();
          id = iter->id();
          idx = idx_run;
        }
      }
      ++idx_run;
    }
    pixels_.push_back(pixesl.pop(idx));
  }
}

void soften() {
	// Since init has already given us an ordered set of pixels all we need to do
	// is loop through the pixels until we hit the resolution we are intrested in
	// softening to. Starting from the end of the vector and looping backwards
	// sounds like the best bet here.
}

void combine(const pixel_union& u) {
	// This method is easy enough. We can simply concatenate the two pixel vectors
	// in each union and re-initialize the map. This may not be the best plan as
	// since both maps are already unions they should conform to the ordering in
	// init. Using this assumption could simplify things.

}

void pixel_union::intersect(const pixel_union& u) {
	// loop through of input union testing first for may_interect with each this
	// union. The next test is contains and if it is then, we keep the pixel. If
	// the test is only may_intersect then we break the pixel up and perform the
	// same test on it's children until contains is true or may_intersect is
	// false.

	pixel_vector kept_pixels;
	for (pixel_iterator iter = u.begin(); iter != u.end(); ++iter) {
	  pixel_vector tmp_pixels = pixel_intersection(*iter);
	  if (!tmp_pixels.empty()) {
	    for (pixel_iterator t_iter = tmp_pixels.begin();
	        t_iter != tmp_pixel.end(); ++t_iter) kept_pixels.push_back(*t_iter);
	  }
	}
	init(kept_pixels);
}

void pixel_union::exclude(const pixel_union& u) {
	// same as above but we want to flip the tests so that if a pixel in our
	// current union is within the input union we want to drop that pixel.

  pixel_vector kept_pixels;
  for (pixel_iterator iter = u.begin(); iter != u.end(); ++iter) {
    pixel_vector tmp_pixels = pixel_exclusion(*iter);
    if (!tmp_pixels.empty()) {
      for (pixel_iterator t_iter = tmp_pixels.begin();
          t_iter != tmp_pixel.end(); ++t_iter) kept_pixels.push_back(*t_iter);
    }
  }
  init(kept_pixels);
}

void pixel_union::init_from_combination(const pixel_union& a,
		const pixel_union& b) {
	// blah blah same as above.
}

void pixel_union::init_from_intersection(const pixel_union& a,
		const pixel_union& b) {

  pixel_vector kept_pixels;
  for (pixel_iterator iter = b.begin(); iter != b.end(); ++iter) {
    pixel_vector tmp_pixels = a.pixel_intersection(*iter);
    if (!tmp_pixels.empty()) {
      for (pixel_iterator t_iter = tmp_pixels.begin();
          t_iter != tmp_pixel.end(); ++t_iter) kept_pixels.push_back(*t_iter);
    }
  }
  init(kept_pixels);
}

void pixel_union::init_from_exclusion(const pixel_union& a,
		const pixel_union& b) {

  pixel_vector kept_pixels;
  for (pixel_iterator iter = b.begin(); iter != b.end(); ++iter) {
    pixel_vector tmp_pixels = a.pixel_exclusion(*iter);
    if (!tmp_pixels.empty()) {
      for (pixel_iterator t_iter = tmp_pixels.begin();
          t_iter != tmp_pixel.end(); ++t_iter) kept_pixels.push_back(*t_iter);
    }
  }
  init(kept_pixels);
}

pixel_vector pixel_union::pixel_intersection(pixel& pix) {
  pixel_vector pixels;
  if (may_intersect(pix)) {
    if (contains(pix)) pix_vect.push_back(pix);
    else if (!pix.is_leaf()) {
      pixel_vector tmp_pixels;
      for (pixel_iterator child_iter = pix.child_begin();
          child_iter != pix.child_end(); ++child_iter) {
        tmp_pixels = pixel_intersection(*child_iter);
        if (!tmp_pixels.empty()) {
          for (pixle_iterator iter = tmp_pixels.begin();
              iter != tmp_pixels.end(); ++iter)
            pixels.push_back(*iter);
        }
      }
    }
  }
  return pixels;
}

pixel_vector pixel_union::pixel_exclusion(pixel& pix) {
  pixel_vector pixels;
  if (!contains(pix)) {
    if (!may_intersect(pix)) pix_vect.push_back(pix);
    else if (!pix.is_leaf()) {
      pixel_vector tmp_pixels;
      for (pixel_iterator child_iter = pix.child_begin();
          child_iter != pix.child_end(); ++child_iter) {
        tmp_pixels = pixel_intersection(*child_iter);
        if (!tmp_pixels.empty()) {
          for (pixle_iterator iter = tmp_pixels.begin();
              iter != tmp_pixels.end(); ++iter)
            pixels.push_back(*iter);
        }
      }
    }
  }
  return pixels;
}




