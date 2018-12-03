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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

#include "particle_population.hpp"
#include "sampling_utils.hpp"

template <typename P>
void multinomial_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);
void stratified_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);
void systematic_resampling(const gsl_rng *random, const vector<double> *probs, int N, unsigned int *indices);

void multinomial_resampling_sorted_uniform(const gsl_rng *random, int N, double *u);
void stratified_resampling_sorted_uniform(const gsl_rng *random, int N, double *u);
void systematic_resampling_sorted_uniform(const gsl_rng *random, int N, double *u);

void multinomial_resampling(const gsl_rng *random, const vector<double> *probs, unsigned int N, unsigned int *indices)
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

void multinomial_resampling_sorted_uniform(const gsl_rng *random, int N, double *uvec)
{
    uniform(random, N, uvec); // N
    sort(uvec, uvec + N); // N log N
}

void stratified_resampling_sorted_uniform(const gsl_rng *random, int N, double *u)
{
    // draw N uniform numbers from (n/N, (n+1)/N]
    double num_samples = (double)N;
    for (int n = 0; n < N; n++) {
        double a = n/num_samples;
        double b = (n+1)/num_samples;
        u[n] = gsl_ran_flat(random, a, b);
    }
}

void systematic_resampling_sorted_uniform(const gsl_rng *random, int N, double *u)
{
    double one_over_N = 1./N;
    double initial_unif = gsl_ran_flat(random, 0, one_over_N);
    double double_N = (double)N;
    
    for (int n = 0; n < N; n++) {
        u[n] = (n - 1)/double_N + initial_unif;
    }
}

#endif /* resampling_h */
