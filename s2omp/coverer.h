/*
 * coverer.h
 *
 *  Created on: Jul 24, 2012
 *      Author: cbmorrison
 */

#ifndef COVERER_H_
#define COVERER_H_

namespace s2omp {

class coverer {

public:
  coverer();
  virtual ~coverer();

  static void get_simple_covering(const bound_interface& bound, int level,
      pixel_vector* pixels);

private:
  int min_level_, max_level_;
  uint32_t max_pixels_;
};

} // end namespace s2omp


#endif /* COVERER_H_ */
