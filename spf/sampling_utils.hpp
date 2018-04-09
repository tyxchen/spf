/*
 * sampling_util.hpp
 *
 *  Created on: 30 Mar 2018
 *      Author: seonghwanjun
 */

#ifndef SRC_SAMPLING_UTILS_HPP_
#define SRC_SAMPLING_UTILS_HPP_

#include <vector>

#include <gsl/gsl_rng.h>

using namespace std;

gsl_rng* get_random(long seed);

// draw one sample and return the index of that sample {1, ..., normalized_probs.size()}
int multinomial(const gsl_rng *random, vector<double> normalized_probs);

// ASSERT: the length of result is same as normalized_probs
void multinomial(const gsl_rng *random, unsigned int N, vector<double> normalized_probs, unsigned int *result);

// ASSERT: the length of indices is N
// each element of indices taking on values {0, ..., normalized_probs.size()}
void multinomial_sample_indices(const gsl_rng *random, unsigned int N, vector<double> normalized_probs, unsigned int *indices);

#endif /* SRC_SAMPLING_UTILS_HPP_ */
