//
//  tssb_problem_spec.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-23.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include "tssb_problem_spec.hpp"

TSSBProblemSpecification::TSSBProblemSpecification(vector<SomaticMutation> data_points, CancerPhyloParameters &params) :
params{params}
{
    this->num_ssms = data_points.size();
    this->initial_state = new PartialCancerPhylogenyState(data_points);
}

unsigned long TSSBProblemSpecification::num_iterations() {
    return num_ssms;
}

std::pair<PartialCancerPhylogenyState*, double> TSSBProblemSpecification::propose_initial(gsl_rng *random) {
    return make_pair(initial_state, 0.0);
}

std::pair<PartialCancerPhylogenyState*, double> TSSBProblemSpecification::propose_next(gsl_rng *random, int t, PartialCancerPhylogenyState curr) {
    return curr.assign_data_point(random, params);
}

TSSBProblemSpecification::~TSSBProblemSpecification() { }

