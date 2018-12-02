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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "numerical_utils.hpp"
#include "particle_population.hpp"
#include "param_proposal.hpp"
#include "pmcmc_options.hpp"
#include "resampling.hpp"
#include "smc.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"
#include "sampling_utils.hpp"

using namespace std;

template <class S, class P> class ParticleMMH
{
    PMCMCOptions *options;
    SMC<S, P> *smc;
    ParamProposal<P> *param_proposal;
    vector<P*> *parameters;
    vector<S*> *states;
public:
    ParticleMMH(PMCMCOptions *options, SMC<S, P> *smc, ParamProposal<P> *param_proposal);
    void run();
    vector<P*> *get_parameters();
};

template <class S, class P>
ParticleMMH<S,P>::ParticleMMH(PMCMCOptions *options, SMC<S, P> *smc, ParamProposal<P> *param_proposal)
{
    this->options = options;
    this->smc = smc;
    this->param_proposal = param_proposal;
    this->parameters = new vector<P*>();
    this->states = new vector<S*>();
}

template <class S, class P>
void ParticleMMH<S,P>::run()
{
    // initialize the parameters
    P *param_curr = param_proposal->sample_from_prior(options->random);
    P *param_new = 0;
    
    S *curr_state = 0;
    S *new_state = 0;
    
    // initialize
    smc->run_smc(*param_curr);
    curr_state = smc->sample(options->random);
    
    double log_Z_old = smc->get_log_marginal_likelihood();
    double log_Z_curr = 0;
    double log_accept_ratio = 0.0;
    double unif = 0.0;
    double max_log_z = DOUBLE_NEG_INF;
    size_t num_accepts = 0;
    size_t num_iters = 1;
    for (size_t i = 0; i < options->num_iterations; i++)
    {
        if (num_iters % 100 == 0 && i < options->burn_in) {
            param_proposal->adapt(num_accepts, num_iters-1);
            num_accepts = 0;
            num_iters = 1;
        }
        parameters->push_back(param_curr);
        states->push_back(curr_state);

        // propose new param: proposal can be quite general
        // as it can combine Gibbs sampling and form MH within MH
        param_new = param_proposal->propose(options->random, param_curr);

        // check that the param_new is in the support
        if (param_proposal->log_prior(param_new) == DOUBLE_NEG_INF) {
            continue;
        }

        // run SMC
        smc->run_smc(*param_new);
        new_state = smc->sample(options->random);
        log_Z_curr = smc->get_log_marginal_likelihood();

        // compute the acceptance probability
        log_accept_ratio = (log_Z_curr - log_Z_old) + (param_proposal->log_proposal(param_curr, param_new) - param_proposal->log_proposal(param_new, param_curr)) + (param_proposal->log_prior(param_new) - param_proposal->log_prior(param_curr));

        unif = log(uniform(options->random));
        if (unif < log_accept_ratio) {
            num_accepts++;
            // accept
            cout << "logZ*: " << log_Z_curr << endl;
            param_curr = param_new;
            curr_state = new_state;
            log_Z_old = log_Z_curr;
            if (log_Z_curr > max_log_z) {
                cout << "new maximum log likelihood is found!" << endl;
                max_log_z = log_Z_curr;
            }
        }
        num_iters++;
    }
}

template <class S, class P>
vector<P*> *ParticleMMH<S,P>::get_parameters()
{
    return parameters;
}

#endif /* pmcmc_hpp */
