//
//  linear_gaussian_model.cpp
//  spf
//
//  x_1 ~ Normal(0, \nu^2)
//  x_t | x_{t-1} ~ Normal(Ax_{t-1} + a, \sigma^2)
//  y_t | x_t ~ Normal(Bx_t + b, \tau^2)
//
//  Created by Seong-Hwan Jun on 2018-07-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_randist.h>

#include "linear_gaussian_model.hpp"

LinearGaussianModel::LinearGaussianModel(vector<double> *obs)
{
    this->obs = obs;
}

unsigned long LinearGaussianModel::num_iterations()
{
    return obs->size();
}

pair<double, double> LinearGaussianModel::propose_initial(gsl_rng *random, LinearGaussianModelParams &params)
{
    double x_new = gsl_ran_gaussian(random, params.nu);
    double logw = log(gsl_ran_gaussian_pdf((*obs)[0] - (params.B*x_new + params.b), params.tau));
    return make_pair(x_new, logw);
}

pair<double, double> LinearGaussianModel::propose_next(gsl_rng *random, int t, double curr, LinearGaussianModelParams &params)
{
    double x_new = gsl_ran_gaussian(random, params.sigma) + (params.A * curr + params.a);
    double logw = log(gsl_ran_gaussian_pdf((*obs)[t] - (params.B*x_new + params.b), params.tau));
    return make_pair(x_new, logw);
}

