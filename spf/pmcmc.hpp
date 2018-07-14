//
//  pmcmc.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-06-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef pmcmc_hpp
#define pmcmc_hpp

#include <math.h>
#include <vector>
#include <gsl/gsl_rng.h>

#include "numerical_utils.hpp"
#include "particle_population.hpp"
#include "param_proposal.hpp"
#include "pmcmc_options.hpp"
#include "resampling.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"
#include "sampling_utils.hpp"

using namespace std;

template <class S, class P> class ParticleMMH
{
    PMCMCOptions &options;
    SMC<S, P> &smc;
    ParamProposal<P> param_proposal;
public:
    ParticleMMH(PMCMCOptions &options, SMC<S, P> &smc, ParamProposal<P> param_specification);
    run();
};

template <class S, class P>
ParticleMMH<S,P> ParticleMMH(PMCMCOptions &options, SMC<S, P> &smc, ParamProposal<P> param_proposal)
{
    this->options = options;
    this->smc = smc;
    this->param_proposal = param_proposal;
}

template <class S, class P>
ParticleMMH<S,P> run()
{
    // initialize the parameters
    P *param_curr = param_proposal->sample_from_prior(options->random);
    P *param_new = 0;
    smc.run_smc(param_curr);
    double logZold = smc.get_log_marginal_likelihood();
    double logZcurr = 0;
    double log_accept_ratio = 0.0;
    double unif = 0.0;
    for (size_t i = 0; i < options->num_iterations; i++)
    {
        // propose new param
        param_new = param_proposal->propose(options->random, param_curr);

        // run SMC
        smc.run_smc(param_new);
        logZcurr = smc.get_log_marginal_likelihood();

        // compute the acceptance probability
        log_accept_ratio = (logZcurr - logZold) + (log_proposal(param_prev, param_curr) - log_proposal(param_curr, param_prev));
        unif = log(uniform(options->random));
        if (unif < accept_ratio) {
            // accept
            param_curr = param_new;
            logZcurr = logZold;
        }
    }
}

#endif /* pmcmc_hpp */
