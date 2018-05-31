//
//  tssb_problem_spec.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-23.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef tssb_problem_spec_h
#define tssb_problem_spec_h

#include "smc_model.hpp"
#include "tssb_state.hpp"

class TSSBProblemSpecification : public ProblemSpecification<PartialCancerPhylogenyState *>
{
    unsigned long num_ssms;
    PartialCancerPhylogenyState *initial_state;
    CancerPhyloParameters &params;
public:
    TSSBProblemSpecification(vector<SomaticMutation> *data_points, CancerPhyloParameters &params);
    unsigned long num_iterations();
    std::pair<PartialCancerPhylogenyState *, double> propose_initial(gsl_rng *random);
    std::pair<PartialCancerPhylogenyState *, double> propose_next(gsl_rng *random, int t, PartialCancerPhylogenyState *curr);
    ~TSSBProblemSpecification();
};

#endif /* tssb_problem_spec_h */
