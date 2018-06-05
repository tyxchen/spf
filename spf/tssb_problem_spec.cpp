//
//  tssb_problem_spec.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-23.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <iostream>
#include "tssb_problem_spec.hpp"

TSSBProblemSpecification::TSSBProblemSpecification(vector<SomaticMutation> *data_points, CancerPhyloParameters &params) :
params{params}
{
    this->num_ssms = data_points->size();
    this->initial_state = new PartialCancerPhylogenyState(data_points);
}

unsigned long TSSBProblemSpecification::num_iterations() {
    return num_ssms;
}

std::pair<PartialCancerPhylogenyState*, double> TSSBProblemSpecification::propose_initial(gsl_rng *random) {
    return propose_next(random, 0, initial_state);
}

std::pair<PartialCancerPhylogenyState*, double> TSSBProblemSpecification::propose_next(gsl_rng *random, int t, PartialCancerPhylogenyState *curr) {
    // make a copy of the current state
    PartialCancerPhylogenyState *new_state = new PartialCancerPhylogenyState(*curr);
    double ret = new_state->sample_next_state(random, params);
    //cout << ret << ", " << new_state->print() << endl;
    //return make_pair(new_state, ret - curr->get_log_likelihood());
    return make_pair(new_state, ret);
}

TSSBProblemSpecification::~TSSBProblemSpecification() { }

