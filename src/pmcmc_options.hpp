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
    size_t num_iterations;
    size_t burn_in = 2000;
    gsl_rng *random;
    PMCMCOptions(long seed, size_t num_iter);
};

PMCMCOptions::PMCMCOptions(long seed, size_t num_iter)
{
    this->random = generate_random_object(seed);
    this->num_iterations = num_iter;
}

#endif /* pmcmc_options_h */

