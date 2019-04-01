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
    double ess_threshold = 0.5;
    bool track_population = false;
    bool resample_last_round = false;
    bool use_SPF = false;
    bool debug = false;

    long main_seed = 1;
    long resampling_seed = 2;
    
    gsl_rng *main_random = 0;
    gsl_rng *resampling_random = 0;

    unsigned int num_threads = 1;
    
    void init();
};

void SMCOptions::init()
{
    if (main_random == 0) {
        main_random = generate_random_object(main_seed);
    }
    if (resampling_random == 0) {
        resampling_random = generate_random_object(resampling_seed);
    }
}

#endif /* smc_options_h */
