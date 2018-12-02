//
//  discrete_hmm_params.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-06-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include "discrete_hmm_params.hpp"

DiscreteHMMParams::DiscreteHMMParams(vector<double> &initial,
                                     vector<vector<double> > &transition,
                                     vector<vector<double> > &emission)
: initial_distn(initial), transition_probs(transition), emission_probs(emission)
{
}
