//
//  discrete_hmm_params.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-06-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include "discrete_hmm_params.hpp"

DiscreteHMMParams::DiscreteHMMParams(vector<double> initial_distn,
                                     vector<vector<double> > transition_probs,
                                     vector<vector<double> > emission_probs)
{
    this->initial_distn = initial_distn;
    this->transition_probs = transition_probs;
    this->emission_probs = emission_probs;
}
