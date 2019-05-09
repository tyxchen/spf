//
//  smc_options.cpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-04-10.
//

#include "smc_options.hpp"

#include "sampling_utils.hpp"

void SMCOptions::init()
{
    if (main_random == 0) {
        main_random = generate_random_object(main_seed);
    }
    if (resampling_random == 0) {
        resampling_random = generate_random_object(resampling_seed);
    }
    
    for (size_t i = 0; i < num_particles; i++) {
        proposal_randoms.push_back(generate_random_object(gsl_rng_get(main_random)));
    }
}

SMCOptions::~SMCOptions()
{
    //delete main_random;
    //delete resampling_random;
    gsl_rng_free(main_random);
    gsl_rng_free(resampling_random);
    for (size_t i = 0; i < num_particles; i++) {
        gsl_rng_free(proposal_randoms[i]);
        //delete proposal_randoms[i];
    }
}
