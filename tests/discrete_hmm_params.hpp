//
//  discrete_hmm_params.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-06-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef discrete_hmm_params_hpp
#define discrete_hmm_params_hpp

#include <vector>

using namespace std;

class DiscreteHMMParams
{
public:
    vector<double> &initial_distn;
    vector<vector<double> > &transition_probs;
    vector<vector<double> > &emission_probs;
    DiscreteHMMParams(vector<double> &initial_distn,
                      vector<vector<double> > &transition_probs,
                      vector<vector<double> > &emission_probs);
};

#endif /* discrete_hmm_params_hpp */
