//
//  cancer_phylo_param.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-25.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef cancer_phylo_param_h
#define cancer_phylo_param_h

#include <string>

using namespace std;

class CancerPhyloParameters
{
public:
    double gamma;
    double lambda;
    double alpha_0;
    double birth_rate;
    double death_rate;
    double weibull_scale;
    double weibull_shape;
    double sequencing_error_prob;
    double alpha(string node_str);
};

#endif /* cancer_phylo_param_h */
