//
//  pg.hpp
//  spf-lib
//
//  Created by Seong-Hwan Jun on 2018-12-04.
//

#ifndef pg_h
#define pg_h

#include <math.h>
#include <vector>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "numerical_utils.hpp"
#include "particle_population.hpp"
#include "pg_proposal.hpp"
#include "pmcmc_options.hpp"
#include "resampling.hpp"
#include "csmc.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"
#include "sampling_utils.hpp"

using namespace std;

template <class S, class P> class ParticleGibbs
{
    PMCMCOptions *options;
    ConditionalSMC<S, P> *csmc;
    PGProposal<S, P> *param_proposal;
    vector<P*> *parameters;
    vector<ParticleGenealogy<S> *> *states; // store the genealogy along with the log weights

public:
    ParticleGibbs(PMCMCOptions *options, ConditionalSMC<S, P> *smc, PGProposal<S, P> *param_proposal);
    void run();
    vector<P*> *get_parameters();
    vector<vector<pair<S, double>> *> *get_states();
};

template <class S, class P>
ParticleGibbs<S,P>::ParticleGibbs(PMCMCOptions *options, ConditionalSMC<S, P> *csmc, PGProposal<S, P> *param_proposal)
{
    this->options = options;
    this->csmc = csmc;
    this->param_proposal = param_proposal;
    this->parameters = new vector<P*>();
    this->states = new vector<ParticleGenealogy<S> *>();
}

template <class S, class P>
void ParticleGibbs<S,P>::run()
{
    // initialize the parameters
    P *param = param_proposal->sample_from_prior(options->random);
    // initialize state
    ParticleGenealogy<S> *genealogy = csmc->initialize(*param);

    double log_Z = csmc->get_log_marginal_likelihood();
    double max_log_z = DOUBLE_NEG_INF;
    for (size_t i = 0; i < options->num_iterations; i++)
    {
        cout << "logZ: " << log_Z << endl;
        if (log_Z > max_log_z) {
            cout << "new maximum log likelihood is found!" << endl;
            max_log_z = log_Z;
        }

        parameters->push_back(param);
        states->push_back(genealogy);
        
        // propose new param: proposal can be quite general
        // as it can combine Gibbs sampling and form MH within MH
        // but it is always accepted
        param = param_proposal->propose(options->random, param, genealogy);
        // run cSMC
        genealogy = csmc->run_csmc(*param, genealogy);
        
        log_Z = csmc->get_log_marginal_likelihood();
    }
}

template <class S, class P>
vector<P*> *ParticleGibbs<S,P>::get_parameters()
{
    return parameters;
}

template <class S, class P>
vector<vector<pair<S, double>> *> *ParticleGibbs<S,P>::get_states()
{
    return states;
}

#endif /* pg_h */
