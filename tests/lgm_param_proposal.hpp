//
//  lgm_param_proposal.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-07-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef lgm_param_proposal_hpp
#define lgm_param_proposal_hpp

#include <stdio.h>

#include "linear_gaussian_model_params.hpp"
#include "param_proposal.hpp"

class LGMParamProposal : public ParamProposal<LinearGaussianModelParams>
{
    double mu_A; double sd_A;
    double mu_B; double sd_B;
    double mu_a; double sd_a;
    double mu_b; double sd_b;
    double alpha_x0; double beta_x0;
    double alpha_x; double beta_x;
    double alpha_y; double beta_y;
public:
    LGMParamProposal(double mu_A, double sd_A,
                     double mu_B, double sd_B,
                     double mu_a, double sd_a,
                     double mu_b, double sd_b,
                     double alpha_x0, double beta_x0,
                     double alpha_x, double beta_x,
                     double alpha_y, double beta_y);
    LinearGaussianModelParams *sample_from_prior(gsl_rng *random);
    LinearGaussianModelParams *propose(gsl_rng *random, LinearGaussianModelParams *curr);
    double log_proposal(LinearGaussianModelParams *curr, LinearGaussianModelParams *prev);
};

#endif /* lgm_param_proposal_hpp */
