//
//  pmcmc_options.cpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-04-10.
//

#include "pmcmc_options.hpp"

PMCMCOptions::PMCMCOptions(long seed, size_t num_iter)
{
    this->random = generate_random_object(seed);
    this->num_iterations = num_iter;
}

PMCMCOptions::~PMCMCOptions()
{
    delete random;
}
