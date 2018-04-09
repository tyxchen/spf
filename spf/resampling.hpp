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
void systematic_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);
void residual_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);

void multinomial_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices)
{
    multinomial_sample_indices(random, N, *probs, indices);
}

void stratified_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices)
{
    // draw N uniform numbers from (n/N, (n+1)/N]
    double num_samples = (double)N;
    double *u = new double[N];
    for (int n = 0; n < N; n++) {
        double a = n/num_samples;
        double b = (n+1)/num_samples;
        u[n] = gsl_ran_flat(random, a, b);
    }
    
    // u's are already sorted -- determine the indices
    double sum = 0.0;
    int idx = 0;
    for (unsigned int n = 0; n < N; n++) {
        while (u[n] > (sum + (*probs)[idx])) {
            sum += (*probs)[idx];
            idx++;
        }
        indices[n] = idx;
    }
    
    delete [] u;
}

void systematic_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices)
{
    // draw a single uniform from (0, 1/N]
    double u = gsl_ran_flat(random, 0, 1./N);
    double u_i = u;
    double inc = 1./N;
    double sum = 0.0;
    int idx = 0;
    for (unsigned int n = 0; n < N; n++) {
        // find the index
        while (u_i > (sum + (*probs)[idx])) {
            sum += (*probs)[idx];
            idx++;
        }
        indices[n] = idx;
        u_i += inc;
    }
}

// NOTE: Implementation not complete -- seems to be introducing bias
 void residual_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices)
 {
 unsigned long population_size = probs->size();
 unsigned int N_i = 0;
 unsigned int R = 0;
 unsigned int idx = 0;
 vector<double> *residual_probs = new vector<double>();
 
 // compute the residuals
 for (int i = 0; i < population_size; i++) {
 double Nw = N * (*probs)[i];
 N_i = (int)floor(Nw);
 residual_probs->push_back(Nw - N_i);
 R += N_i;
 for (int j = 0; j < N_i; j++) {
 indices[idx] = i;
 idx++;
 }
 }
 normalize_destructively(*residual_probs);
 multinomial_sample_indices(random, N-R, *residual_probs, indices + R);
 
 // shuffle the indices?
 gsl_ran_shuffle(random, indices, N, sizeof(unsigned int));
 
 delete residual_probs;
 }

#endif /* resampling_h */
