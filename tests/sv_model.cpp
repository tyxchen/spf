//
//  sv_model.cpp
//  sv
//
//  Created by Seong-Hwan Jun on 2018-11-25.
//

#include "sv_model.hpp"

#include <iostream>
#include <math.h>
#include <gsl/gsl_randist.h>

SVModel::SVModel(vector<double> &obs)
{
    this->obs = obs;
}

unsigned long SVModel::num_iterations()
{
    return obs.size();
}

Particle<double> *SVModel::propose_initial(gsl_rng *random, SVModelParams &params)
{
    double x1 = gsl_ran_gaussian(random, params.sigma);
    double obs_log_lik = log(gsl_ran_gaussian_pdf(obs[0], sqrt(exp(x1)) * params.beta));
    return make_particle(&x1, obs_log_lik);
}

Particle<double> *SVModel::propose_next(gsl_rng *random, int t, double &curr, SVModelParams &params)
{
    double xt = params.phi * curr + gsl_ran_gaussian(random, params.sigma);
    double obs_log_lik = log(gsl_ran_gaussian_pdf(obs[t], sqrt(exp(xt)) * params.beta));
    return make_particle(&xt, obs_log_lik);
}

void SVModel::generate_data(gsl_rng *random, size_t T, SVModelParams &params, vector<double> &latent, vector<double> &obs)
{
    for (size_t t = 0; t < T; t++)
    {
        double x, y;
        if (t == 0) {
            x = gsl_ran_gaussian(random, params.sigma);
        } else {
            x = params.phi * latent[t-1] + gsl_ran_gaussian(random, params.sigma);
        }
        y = gsl_ran_gaussian(random, sqrt(exp(x)) * params.beta);
        latent.push_back(x);
        obs.push_back(y);
    }
}

SVModel::~SVModel()
{
    
}
