//
//  resampling.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-04-07.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef resampling_h
#define resampling_h

#include <math.h>

#include "particle_population.hpp"
#include "sampling_utils.hpp"
#include "smc_options.hpp"

template <typename P>
void multinomial_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);
void stratified_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);
void residual_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);
void systematic_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);

void multinomial_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices)
{
    sample_indices(random, N, *probs, indices);
}

void stratified_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices)
{
    
}

void systematic_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices)
{
    
}

void residual_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices)
{
    
}

void systematic_resampling()
{
    
}

#endif /* resampling_h */
