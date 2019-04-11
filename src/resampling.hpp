//
//  resampling.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-04-07.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef resampling_h
#define resampling_h

#include <algorithm>
#include <cmath>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

using namespace std;

void multinomial_resampling(const gsl_rng *random, const vector<double> *probs, unsigned int N, unsigned int *indices);
void stratified_resampling(const gsl_rng *random, const vector<double> *probs, unsigned int N, unsigned int *indices);
void systematic_resampling(const gsl_rng *random, const vector<double> *probs, unsigned int N, unsigned int *indices);

void multinomial_resampling_sorted_uniform(const gsl_rng *random, unsigned int N, double *u);
void stratified_resampling_sorted_uniform(const gsl_rng *random, unsigned int N, double *u);
void systematic_resampling_sorted_uniform(const gsl_rng *random, unsigned int N, double *u);

#endif /* resampling_h */
