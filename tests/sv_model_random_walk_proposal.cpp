//
//  sv_model_param_proposal.cpp
//  sv
//
//  Created by Seong-Hwan Jun on 2018-11-25.
//

#include <math.h>
#include <gsl/gsl_randist.h>
#include "numerical_utils.hpp"
#include "sv_model_random_walk_proposal.hpp"

SVModelParams *SVModelRandomWalkProposal::sample_from_prior(gsl_rng *random)
{
    SVModelParams *param = new SVModelParams(1.0, 0.16, gsl_ran_flat(random, 0.0, 2.0));
    return param;
}

SVModelParams *SVModelRandomWalkProposal::propose(gsl_rng *random, SVModelParams *curr)
{
    double beta = curr->beta + gsl_ran_gaussian(random, 0.02);
    SVModelParams *param = new SVModelParams(1.0, 0.16, beta);
    return param;
}

double SVModelRandomWalkProposal::log_proposal(SVModelParams *curr, SVModelParams *prev)
{
    return 0.0;
}

double SVModelRandomWalkProposal::log_prior(SVModelParams *curr)
{
    if (curr->beta < 0 || curr->beta > 2) {
        return DOUBLE_NEG_INF;
    }
    return log(gsl_ran_flat_pdf(curr->beta, 0, 2));
}

void SVModelRandomWalkProposal::adapt(size_t num_accepts, size_t curr_iter)
{
    
}
