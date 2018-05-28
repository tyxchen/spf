//
//  cancer_phylo_param.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-28.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <cmath>
#include "cancer_phylo_param.hpp"

double CancerPhyloParameters::alpha(string node_str)
{
    size_t j = node_str.length();
    return pow(this->lambda, j) * this->alpha_0;
}
