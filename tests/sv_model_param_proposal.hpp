//
//  sv_model_param_proposal.hpp
//  sv
//
//  Created by Seong-Hwan Jun on 2018-11-25.
//

#ifndef sv_model_param_proposal_hpp
#define sv_model_param_proposal_hpp

#include <stdio.h>

#include "param_proposal.hpp"
#include "sv_model_params.hpp"

class SVModelParamProposal : public ParamProposal<SVModelParams>
{
public:
    SVModelParams *sample_from_prior(gsl_rng *random);
    SVModelParams *propose(gsl_rng *random, SVModelParams *curr);
    double log_proposal(SVModelParams *curr, SVModelParams *prev);
    double log_prior(SVModelParams *curr); // log p(curr)
    void adapt(size_t num_accepts, size_t curr_iter);
};

#endif /* sv_model_param_proposal_hpp */
