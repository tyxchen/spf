//
//  linear_gaussian_model.hpp
//  spf
//
//  simplest possible model (y,x are both scalars)
//  y_t | x_t ~ N(Bx_t + b, \sigma^2)
//  x_t | x_{t-1} ~ N(Ax_t + a, \tau^2 I)
//  y_t, a, b, \sigma^2 are scalars, A,x_t are vectors of dim d
//
//  Created by Seong-Hwan Jun on 2018-07-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef linear_gaussian_model_hpp
#define linear_gaussian_model_hpp

#include <stdio.h>
#include <vector>

#include "linear_gaussian_model_params.hpp"
#include "smc_model.hpp"

using namespace std;

class LinearGaussianModel : public ProblemSpecification<double, LinearGaussianModelParams>
{
    vector<double> *obs;
public:
    LinearGaussianModel(vector<double> *obs);
    unsigned long num_iterations();
    // propose new sample and return it with the log weight
    pair<double, double> propose_initial(gsl_rng *random, LinearGaussianModelParams &params);
    pair<double, double> propose_next(gsl_rng *random, int t, double curr, LinearGaussianModelParams &params);
};

#endif /* linear_gaussian_model_hpp */
