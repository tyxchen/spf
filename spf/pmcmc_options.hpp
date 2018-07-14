//
//  pmcmc_options.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-06-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef pmcmc_options_h
#define pmcmc_options_h

class PMCMCOptions
{
public:
    size_t num_iterations;
    gsl_rng *random;
};

#endif /* pmcmc_options_h */
