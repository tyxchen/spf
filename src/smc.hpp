/*
 * smc.hpp
 *
 * T: the latent variable type
 *
 *  Created on: 6 Mar 2018
 *      Author: seonghwanjun
 */

#ifndef SRC_SMC_HPP_
#define SRC_SMC_HPP_

#include <vector>
#include <gsl/gsl_rng.h>

#include "numerical_utils.hpp"
#include "particle_population.hpp"
#include "resampling.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"
#include "sampling_utils.hpp"

using namespace std;

// S: state
// P: parameter
template <class S, class P> class SMC
{
    SMCOptions *options;
	ProblemSpecification<S, P> *proposal; // pointer to proposal object to be passed into SMC constructor

    double log_marginal_likelihood = 0;
    pair<double, ParticlePopulation<S>*> propose(gsl_rng *random, ParticlePopulation<S> *pop, int iter, int num_proposals, P &params);
    ParticlePopulation<S> *resample(const gsl_rng *random, SMCOptions::ResamplingScheme resampling_scheme, ParticlePopulation<S> *pop, int N);
    vector<ParticlePopulation<S> *> *populations = 0;

public:
	SMC(ProblemSpecification<S, P> *proposal, SMCOptions *options);
	void run_smc(P &params);
    S *sample(gsl_rng *random);
    double get_log_marginal_likelihood();
    inline ParticlePopulation<S>* get_curr_population() { return populations->back(); }
    inline ParticlePopulation<S> *get_population(size_t i) { return populations->at(i); }
    ~SMC();
};

template <class S, class P>
SMC<S,P>::SMC(ProblemSpecification<S, P> *proposal, SMCOptions *options)
{
    this->proposal = proposal;
    this->options = options;
    this->options->init();
    this->populations = new vector<ParticlePopulation<S> *>();
}

template <class S, class P>
void SMC<S,P>::run_smc(P &params)
{
    unsigned long R = proposal->num_iterations();

    pair<double, ParticlePopulation<S>*> ret_val;
    ParticlePopulation<S> *curr_pop = 0;
    gsl_rng *random = generate_random_object(options->main_seed);

    log_marginal_likelihood = 0.0;
    //double log_num_particles = log(options->num_particles);
    for (int r = 0; r < R; r++)
    {
        if (options->debug)
            cout << "iter: " << r << endl;
        ret_val = propose(random, curr_pop, r, options->num_particles, params);
        curr_pop = ret_val.second;
        populations->push_back(curr_pop);
        if (!options->track_population && r > 0) {
            if (!(*populations)[r-1])
                delete (*populations)[r-1]; // delete the ParticlePopulation (particles, log_weights, and normalized_weights)
        }
        log_marginal_likelihood += ret_val.first;
        if (r == R - 1 && options->resample_last_round == false) {
            break;
        }

        if (options->ess_threshold == 0.0) // no resampling (essneitally SIS)
            continue;
        if (options->ess_threshold >= 1.0 || curr_pop->get_ess() <= options->ess_threshold) {
            // resample
            curr_pop = resample(random, options->resampling_scheme, curr_pop, options->num_particles);
            delete (*populations)[r]; // delete pre-resampling population
            (*populations)[r] = curr_pop;
        }
    }
}

template <class S, class P>
pair<double, ParticlePopulation<S>*> SMC<S,P>::propose(gsl_rng *random, ParticlePopulation<S> *pop, int iter, int num_proposals, P &params)
{
    double log_Z_ratio = DOUBLE_NEG_INF;
    vector<S> *curr_particles = 0;
    vector<double> *curr_log_weights = 0;
    bool is_resampled = true;
    if (pop != 0) {
        curr_particles = pop->get_particles();
        curr_log_weights = pop->get_log_weights();
        is_resampled = pop->is_resampled();
    }

    vector<S> *new_particles = new vector<S>(num_proposals);
    vector<double> *new_log_weights = new vector<double>(num_proposals);

    double log_norm = DOUBLE_NEG_INF;
    Particle<S> *ret;
    double log_N = log(num_proposals);
    for (int n = 0; n < num_proposals; n++)
    {
        if (pop == 0) {
            ret = proposal->propose_initial(random, params);
        } else {
            ret = proposal->propose_next(random, iter, (*curr_particles)[n], params);
        }
        (*new_particles)[n] = *ret->get_state();
        double log_alpha = ret->get_log_alpha();
        if (is_resampled) {
            (*new_log_weights)[n] = log_alpha - log_N;
        } else {
            (*new_log_weights)[n] = log_alpha + (*curr_log_weights)[n];
        }
        log_norm = log_add(log_norm, log_alpha);
        log_Z_ratio = log_add(log_Z_ratio, log_alpha);
    }
    log_Z_ratio -= log(num_proposals);
    ParticlePopulation<S> *new_pop = new ParticlePopulation<S>(new_particles, new_log_weights, log_norm, log_Z_ratio);
    return make_pair(log_Z_ratio, new_pop);
}

template <class S, class P>
ParticlePopulation<S>* SMC<S,P>::resample(const gsl_rng *random, SMCOptions::ResamplingScheme resampling_scheme, ParticlePopulation<S> *pop, int N)
{
    unsigned int *indices = new unsigned int[N];
    switch (resampling_scheme)
    {
        case SMCOptions::ResamplingScheme::MULTINOMIAL:
            multinomial_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
        case SMCOptions::ResamplingScheme::STRATIFIED:
            stratified_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
        case SMCOptions::ResamplingScheme::SYSTEMATIC:
            systematic_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
    }

    vector<S> *curr_particles = pop->get_particles();

    vector<S> *particles = new vector<S>();
    for (int n = 0; n < N; n++)
    {
        particles->push_back((*curr_particles)[indices[n]]);
    }

    ParticlePopulation<S> *new_pop = new ParticlePopulation<S>(particles);

    delete [] indices;
    return new_pop;
}

template <class S, class P>
double SMC<S,P>::get_log_marginal_likelihood()
{
    return log_marginal_likelihood;
}

template <class S, class P>
S *SMC<S,P>::sample(gsl_rng *random)
{
    ParticlePopulation<S> *pop = get_curr_population(); // it could be either normalized already or not (doesn't matter either way)
    int N = 1;
    unsigned int *indices = new unsigned int[N];
    switch (options->resampling_scheme)
    {
        case SMCOptions::ResamplingScheme::MULTINOMIAL:
            multinomial_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
        case SMCOptions::ResamplingScheme::STRATIFIED:
            stratified_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
        case SMCOptions::ResamplingScheme::SYSTEMATIC:
            systematic_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
    }
    
    vector<S> *curr_particles = pop->get_particles();
    S *sample = new S((*curr_particles)[indices[0]]); // use default copy constructor to create a pointer object
    return sample;
}

template <class S, class P>
SMC<S,P>::~SMC()
{
}
#endif /* SRC_SMC_HPP_ */
