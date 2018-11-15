//
//  smc_options.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-04-05.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef smc_options_h
#define smc_options_h

class SMCOptions
{
public:
    unsigned int num_particles = 1000;
    unsigned int max_virtual_particles = 1000000;

    enum ResamplingScheme {
        MULTINOMIAL = 0,
        STRATIFIED = 1,
        SYSTEMATIC = 2
    };
    ResamplingScheme resampling_scheme = MULTINOMIAL;
    double essThreshold = 0.5;
    bool track_population = false;
    bool resample_last_round = false;
    bool debug = true;

    gsl_rng *resampling_random;
    gsl_rng *main_random;

    // at some point, look into implementing multi-threading
    // unsigned int num_threads = 4;
};

#endif /* smc_options_h */