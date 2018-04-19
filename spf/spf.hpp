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

template <class P> class SPF
{
    SMCOptions *options;
    ProblemSpecification<P> *proposal;
    double log_marginal_likelihood = 0;
    CompactParticlePopulation<P> propose_compact_population(PermutationStream &stream, ParticlePopulation<P> *prev_pop, int iter);
    ParticlePopulation<P> *contraction(PermutationStream &stream, ParticlePopulation<P> *prev_pop, CompactParticlePopulation<P> &compact_pop, int N, int iter);
    vector<ParticlePopulation<P> *> *populations = 0;

public:
    SPF(ProblemSpecification<P> *proposal, SMCOptions *options);
    void run_spf();
    
    double get_log_marginal_likelihood();
    inline ParticlePopulation<P>* get_curr_population() { return populations->back(); }
    ~SPF();
};

template <class P>
SPF<P>::SPF(ProblemSpecification<P> *proposal, SMCOptions *options)
{
    this->proposal = proposal;
    this->options = options;
    this->populations = new vector<ParticlePopulation<P> *>();
}

template <class P>
void SPF<P>::run_spf()
{
    unsigned long R = proposal->num_iterations();

    ParticlePopulation<P> *curr_pop = 0;

    log_marginal_likelihood = 0.0;
    for (int r = 0; r < R; r++)
    {
        PermutationStream stream(options->num_particles, gsl_rng_get(options->main_random));

        // expansion
        CompactParticlePopulation<P> compact_pop = propose_compact_population(stream, curr_pop, r);
        log_marginal_likelihood += compact_pop.logZ();
        cout << "ESS: " << compact_pop.ess() << ", sum_of_sq_weights: " << compact_pop.get_log_sum_of_square_weights() << ", nParticles: " << compact_pop.get_num_particles() << ", logZ: " << log_marginal_likelihood << endl;

        // contraction
        stream.reset();
        curr_pop = contraction(stream, curr_pop, compact_pop, options->num_particles, r);
        populations->push_back(curr_pop);
    }
}

template<class P>
CompactParticlePopulation<P> SPF<P>::propose_compact_population(PermutationStream &stream, ParticlePopulation<P> *pop, int iter)
{
    unsigned int idx = 0;
    CompactParticlePopulation<P> compact_pop;
    pair<int, double> ret;

    vector<P> *curr_particles = 0;
    if (pop != 0) {
        curr_particles = pop->get_particles();
    }
    gsl_rng *random = stream.get_random();
    while (compact_pop.get_num_particles() < options->num_particles ||
            (compact_pop.get_num_particles() < options->max_virtual_particles &&
             compact_pop.ess()/options->num_particles < options->essThreshold))
    {
        idx = stream.pop();
        if (pop == 0) {
            ret = proposal->propose_initial(random);
        } else {
            ret = proposal->propose_next(random, iter, (*curr_particles)[idx]);
        }
        compact_pop.add_weight(ret.second);
    }
    
    return compact_pop;
}

template <class P>
ParticlePopulation<P> *SPF<P>::contraction(PermutationStream &stream, ParticlePopulation<P> *pop, CompactParticlePopulation<P> &compact_pop, int N, int iter)
{
    // generate N uniform numbers and sort it
    // propagate particles one at a time, using the log norm to determine, which particle is to be stored
    // draw N uniform values
    unsigned int idx = 0;
    double normalized_partial_sum = 0.0;
    double *sorted_uniform = 0;
    pair<int, double> ret;
    double log_norm = compact_pop.get_log_sum_weights();
    
    vector<P> *new_particles = new vector<P>(N);
    vector<double> *new_log_weights = new vector<double>(N);

    vector<P> *curr_particles = 0;
    if (pop != 0) {
        curr_particles = pop->get_particles();
    }

    switch (options->resampling)
    {
        case SMCOptions::ResamplingScheme::MULTINOMIAL:
            sorted_uniform = multinomial_resampling(options->resampling_random, N);
            break;
        case SMCOptions::ResamplingScheme::STRATIFIED:
            sorted_uniform = stratified_resampling(options->resampling_random, N);
            break;
        case SMCOptions::ResamplingScheme::SYSTEMATIC:
            sorted_uniform = systematic_resampling(options->resampling_random, N);
            break;
    }

    gsl_rng *random = stream.get_random();
    CompactParticlePopulation<P> sanity_check;
    double log_weight = log(1./N);
    for (int n = 0; n < N; n++)
    {
        while (normalized_partial_sum < sorted_uniform[n]) {
            idx = stream.pop();
            if (pop == 0) {
                ret = proposal->propose_initial(random);
            } else {
                ret = proposal->propose_next(random, iter, (*curr_particles)[idx]);
            }
            normalized_partial_sum += exp(ret.second - log_norm);
            sanity_check.add_weight(ret.second);
        }
        (*new_particles)[n] = ret.first;
        (*new_log_weights)[n] = log_weight;
    }
    
    for (int j = 0; j < (compact_pop.get_num_particles() - sanity_check.get_num_particles()); j++) {
        idx = stream.pop();
        if (pop == 0) {
            ret = proposal->propose_initial(random);
        } else {
            ret = proposal->propose_next(random, iter, (*curr_particles)[idx]);
        }
        sanity_check.add_weight(ret.second);
    }
    
    cout << log_norm << " vs " << sanity_check.get_log_sum_weights() << endl;
    
    ParticlePopulation<P> *new_pop = new ParticlePopulation<P>(new_particles, new_log_weights);
    return new_pop;

    delete [] sorted_uniform;
}

template <class P>
double SPF<P>::get_log_marginal_likelihood()
{
    return log_marginal_likelihood;
}

template <class P>
SPF<P>::~SPF()
{
}


#endif /* spf_h */

