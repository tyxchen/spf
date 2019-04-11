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
}

SMCOptions::~SMCOptions()
{
    delete main_random;
    delete resampling_random;
}
