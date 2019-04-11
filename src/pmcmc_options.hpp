//
//  pmcmc_options.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-06-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef pmcmc_options_h
#define pmcmc_options_h

#include "sampling_utils.hpp"

class PMCMCOptions
{
public:
    bool verbose = false;
    size_t burn_in = 0;
    gsl_rng *random;
    size_t num_iterations;
    PMCMCOptions(long seed, size_t num_iter);
    
    ~PMCMCOptions();
};

#endif /* pmcmc_options_h */

