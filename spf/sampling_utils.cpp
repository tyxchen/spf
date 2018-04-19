/*
 * sampling_util.cpp
 *
 *  Created on: 30 Mar 2018
 *      Author: seonghwanjun
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "sampling_utils.hpp"

gsl_rng* generate_random_object(long seed)
{
	// initialize random
    const gsl_rng_type * random_type = gsl_rng_default;
	gsl_rng *random;

	/* create a generator chosen by the
	     environment variable GSL_RNG_TYPE */
	gsl_rng_env_setup();

	random = gsl_rng_alloc(random_type);
	gsl_rng_set(random, seed);
	return random;
}

int multinomial(const gsl_rng *random, vector<double> normalized_probs)
{
    double u = gsl_rng_uniform(random);
    double sum = 0.0;
    for (int i = 0; i < normalized_probs.size(); i++)
    {
        sum += normalized_probs[i];
        if (u <= sum) {
            return i;
        }
    }
    
    // probs does not sum to 1: check this and handle the error
    return -1;
}

void uniform(const gsl_rng *random, unsigned int N, double *ret)
{
    for (int n = 0; n < N; n++)
    {
        ret[n] = gsl_rng_uniform(random);
    }
}

void multinomial(const gsl_rng *random, unsigned int N, vector<double> normalized_probs, unsigned int *result)
{
    double probs[normalized_probs.size()];
    std::copy(normalized_probs.begin(), normalized_probs.end(), probs);
    gsl_ran_multinomial(random, N, N, probs, result);
}

void multinomial_sample_indices(const gsl_rng *random, unsigned int N, vector<double> normalized_probs, unsigned int *indices)
{
    // draw N uniform values
    double *uvec = new double[N];
    uniform(random, N, uvec); // N
    sort(uvec, uvec + N); // N log N
    double sum = 0.0;
    int idx = 0;
    for (int n = 0; n < N; n++)
    {
        while (uvec[n] > sum + normalized_probs[idx]) {
            sum += normalized_probs[idx];
            idx++;
        }
        indices[n] = idx;
    }

    delete [] uvec;
}
