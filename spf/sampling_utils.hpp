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

int multinomial(const gsl_rng *random, vector<double> normalized_probs);
void multinomial(const gsl_rng *random, unsigned int N, vector<double> normalized_probs, unsigned int *result); // ASSERT: the length of result is same as normalized_probs
void sample_indices(const gsl_rng *random, unsigned int N, vector<double> normalized_probs, unsigned int *indices); // ASSERT: the length of indices is N
gsl_rng* get_random(long seed);


#endif /* SRC_SAMPLING_UTILS_HPP_ */
