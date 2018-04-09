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
    enum ResamplingScheme {
        MULTINOMIAL = 0,
        STRATIFIED = 1,
        SYSTEMATIC = 2,
        RESIDUAL = 3
    };
    ResamplingScheme resampling = MULTINOMIAL;
    double essThreshold = 0.5;
    bool track_population = false;
    
    // at some point, look into implementing multi-threading
    // unsigned int num_threads = 4;
};

#endif /* smc_options_h */
