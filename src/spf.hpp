//
//  spf.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-04-17.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef SRC_SPF_HPP_
#define SRC_SPF_HPP_

#include <vector>
#include <gsl/gsl_rng.h>

#include "compact_particle_population.hpp"
#include "numerical_utils.hpp"
#include "particle_population.hpp"
#include "permutation_stream.hpp"
#include "resampling.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"
#include "sampling_utils.hpp"

using namespace std;

template <class S, class P> class SPF
{
    SMCOptions &options;
    ProblemSpecification<S, P> &proposal;
    double log_marginal_likelihood = 0;
    CompactParticlePopulation<S> propose_compact_population(PermutationStream &stream, ParticlePopulation<S> *prev_pop, int iter, P &params);
    ParticlePopulation<S> *contraction(PermutationStream &stream, ParticlePopulation<S> *prev_pop, CompactParticlePopulation<S> &compact_pop, int N, int iter, P &params);
    vector<ParticlePopulation<S> *> *populations = 0;

    gsl_rng *main_random;
    gsl_rng *resampling_random;

public:
    SPF(ProblemSpecification<S, P> &roposal, SMCOptions &options);
    void run_spf(P &params);
    
    double get_log_marginal_likelihood();
    inline ParticlePopulation<S>* get_curr_population() { return populations->back(); }
    ~SPF();
};

template <class S, class P>
SPF<S,P>::SPF(ProblemSpecification<S, P> &proposal, SMCOptions &options) :
options(options), proposal(proposal)
{
    this->populations = new vector<ParticlePopulation<S> *>();
}

template <class S, class P>
void SPF<S,P>::run_spf(P &params)
{
    unsigned long R = proposal.num_iterations();
    options.init();

    main_random = options.main_random;
    resampling_random = options.resampling_random;

    ParticlePopulation<S> *curr_pop = 0;

    log_marginal_likelihood = 0.0;
    long seed;
    PermutationStream stream(options.num_particles);
    for (int r = 0; r < R; r++)
    {
        seed = gsl_rng_get(main_random);
        stream.set_seed(seed);

        // expansion
        CompactParticlePopulation<S> compact_pop = propose_compact_population(stream, curr_pop, r, params);
        log_marginal_likelihood += compact_pop.logZ();
        if (options.debug) {
            cout << "ESS: " << compact_pop.ess() << ", sum_of_sq_weights: " << compact_pop.get_log_sum_of_square_weights() << ", nParticles: " << compact_pop.get_num_particles() << ", logZ: " << log_marginal_likelihood << endl;
        }

        // contraction
        stream.reset();
        curr_pop = contraction(stream, curr_pop, compact_pop, options.num_particles, r, params);
        populations->push_back(curr_pop);
    }
}

template <class S, class P>
CompactParticlePopulation<S> SPF<S,P>::propose_compact_population(PermutationStream &stream, ParticlePopulation<S> *pop, int iter, P &params)
{
    unsigned int idx = 0;
    CompactParticlePopulation<S> compact_pop;

    vector<shared_ptr<S> > *curr_particles = 0;
    if (pop != 0) {
        curr_particles = pop->get_particles();
    }
    gsl_rng *random = stream.get_random();
    while (compact_pop.get_num_particles() < options.num_particles ||
            (compact_pop.get_num_particles() < options.max_virtual_particles &&
             compact_pop.ess()/options.num_particles < options.ess_threshold))
    {
        double log_w = 0.0;
        idx = stream.pop();
        if (pop == 0) {
            proposal.propose_initial(random, log_w, params);
        } else {
            S &parent = *curr_particles->at(idx).get();
            proposal.propose_next(random, iter, parent, log_w, params);
        }
        compact_pop.add_weight(log_w);
    }

    return compact_pop;
}

template <class S, class P>
ParticlePopulation<S> *SPF<S,P>::contraction(PermutationStream &stream, ParticlePopulation<S> *pop, CompactParticlePopulation<S> &compact_pop, int N, int iter, P &params)
{
    // generate N uniform numbers and sort it
    // propagate particles one at a time, using the log norm to determine, which particle is to be stored
    // draw N uniform values
    unsigned int idx = 0;
    double normalized_partial_sum = 0.0;
    double *sorted_uniform = new double[N];
    double log_norm = compact_pop.get_log_sum_weights();

    vector<shared_ptr<S> > *new_particles = new vector<shared_ptr<S> >(N);
    vector<shared_ptr<S> > *curr_particles = 0;
    if (pop != 0) {
        curr_particles = pop->get_particles();
    }

    switch (options.resampling_scheme)
    {
        case SMCOptions::ResamplingScheme::MULTINOMIAL:
            multinomial_resampling_sorted_uniform(resampling_random, N, sorted_uniform);
            break;
        case SMCOptions::ResamplingScheme::STRATIFIED:
            stratified_resampling_sorted_uniform(resampling_random, N, sorted_uniform);
            break;
        case SMCOptions::ResamplingScheme::SYSTEMATIC:
            systematic_resampling_sorted_uniform(resampling_random, N, sorted_uniform);
            break;
    }

    gsl_rng *random = stream.get_random();
    CompactParticlePopulation<S> sanity_check;
    shared_ptr<S> s;
    for (int n = 0; n < N; n++)
    {
        while (normalized_partial_sum < sorted_uniform[n]) {
            idx = stream.pop();
            double log_w = 0.0;
            if (pop == 0) {
                s = proposal.propose_initial(random, log_w, params);
            } else {
                S &parent = *curr_particles->at(idx).get();
                s = proposal.propose_next(random, iter, parent, log_w, params);
            }
            normalized_partial_sum += exp(log_w - log_norm);
            sanity_check.add_weight(log_w);
        }
        (*new_particles)[n] = s;
    }
    
    int num_iter_left = compact_pop.get_num_particles() - sanity_check.get_num_particles();
    for (int j = 0; j < num_iter_left; j++) {
        idx = stream.pop();
        double log_w = 0.0;
        if (pop == 0) {
            proposal.propose_initial(random, log_w, params);
        } else {
            S &parent = *curr_particles->at(idx).get();
            proposal.propose_next(random, iter, parent, log_w, params);
        }
        sanity_check.add_weight(log_w);
    }

    cout << "|" << log_norm << " - " << sanity_check.get_log_sum_weights() << "| = " << abs(log_norm - sanity_check.get_log_sum_weights()) << endl;

    //ParticlePopulation<S> *new_pop = new ParticlePopulation<S>(new_particles, new_log_weights);
    ParticlePopulation<S> *new_pop = new ParticlePopulation<S>(new_particles);
    delete [] sorted_uniform;
    return new_pop;

}

template <class S, class P>
double SPF<S,P>::get_log_marginal_likelihood()
{
    return log_marginal_likelihood;
}

template <class S, class P>
SPF<S,P>::~SPF()
{
}


#endif /* spf_h */

