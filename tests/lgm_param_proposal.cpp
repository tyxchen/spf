//
//  lgm_param_proposal.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-07-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <gsl/gsl_randist.h>

#include "lgm_param_proposal.hpp"

LGMParamProposal::LGMParamProposal(double mu_A, double sd_A,
                                   double mu_B, double sd_B,
                                   double mu_a, double sd_a,
                                   double mu_b, double sd_b,
                                   double alpha_x0, double beta_x0,
                                   double alpha_x, double beta_x,
                                   double alpha_y, double beta_y)
{
    // Normal dist
    this->mu_A = mu_A;
    this->sd_A = sd_A;

    this->mu_B = mu_B;
    this->sd_B = sd_B;

    this->mu_a = mu_a;
    this->sd_a = sd_a;

    this->mu_b = mu_b;
    this->sd_b = sd_b;

    // Gamma (or IG) dist
    this->alpha_x0 = alpha_x0;
    this->beta_x0 = beta_x0;

    this->alpha_x = alpha_x;
    this->beta_x = beta_x;

    this->alpha_y = alpha_y;
    this->beta_y = beta_y;
}

LinearGaussianModelParams *LGMParamProposal::sample_from_prior(gsl_rng *random)
{
    LinearGaussianModelParams *param = new LinearGaussianModelParams();
    param->A = gsl_ran_gaussian(random, sd_A) + mu_A;
    param->B = gsl_ran_gaussian(random, sd_B) + mu_B;

    param->a = gsl_ran_gaussian(random, sd_a) + mu_a;
    param->b = gsl_ran_gaussian(random, sd_b) + mu_b;

    param->nu = gsl_ran_gamma(random, alpha_x0, beta_x0);
    param->sigma = gsl_ran_gamma(random, alpha_x, beta_x);
    param->tau = gsl_ran_gamma(random, alpha_y, beta_y);

    return param;
}

// independent sampler (i.e., does not use the current value)
LinearGaussianModelParams *LGMParamProposal::propose(gsl_rng *random, LinearGaussianModelParams *curr)
{
    LinearGaussianModelParams *param = new LinearGaussianModelParams();
    param->A = gsl_ran_gaussian(random, sd_A) + mu_A;
    param->B = gsl_ran_gaussian(random, sd_B) + mu_B;
    
    param->a = gsl_ran_gaussian(random, sd_a) + mu_a;
    param->b = gsl_ran_gaussian(random, sd_b) + mu_b;
    
    param->nu = gsl_ran_gamma(random, alpha_x0, beta_x0);
    param->sigma = gsl_ran_gamma(random, alpha_x, beta_x);
    param->tau = gsl_ran_gamma(random, alpha_y, beta_y);
    
    return param;
}

double LGMParamProposal::log_proposal(LinearGaussianModelParams *curr, LinearGaussianModelParams *prev)
{
    // for independent sampler, it is just the log of the prior
    double log_density = gsl_ran_gaussian_pdf(curr->A - mu_A, sd_A);
    log_density += gsl_ran_gaussian_pdf(curr->B - mu_B, sd_B);
    log_density += gsl_ran_gaussian_pdf(curr->a - mu_a, sd_a);
    log_density += gsl_ran_gaussian_pdf(curr->b - mu_b, sd_b);
    
    log_density += gsl_ran_gamma_pdf(curr->nu, alpha_x0, beta_x0);
    log_density += gsl_ran_gamma_pdf(curr->sigma, alpha_x, beta_x);
    log_density += gsl_ran_gamma_pdf(curr->tau, alpha_y, beta_y);
    
    return log_density;
}

