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

template <class P> class SMC
{
    SMCOptions *options;
	ProblemSpecification<P> *proposal; // pointer to proposal object to be passed into SMC constructor
    double log_marginal_likelihood = 0;
    ParticlePopulation<P> *propose(gsl_rng *random, ParticlePopulation<P> *pop, int iter, int num_proposals);
    ParticlePopulation<P> *resample(const gsl_rng *random, SMCOptions::ResamplingScheme resampling_scheme, ParticlePopulation<P> *pop, int N);
    vector<ParticlePopulation<P> *> *populations = 0;

public:
	SMC(ProblemSpecification<P> *proposal, SMCOptions *options);
	void run_smc(gsl_rng *random, unsigned int num_particles);
    double get_log_marginal_likelihood();
    inline ParticlePopulation<P>* get_curr_population() { return populations->back(); }
    ~SMC();
};

template <class P>
SMC<P>::SMC(ProblemSpecification<P> *proposal, SMCOptions *options)
{
    this->proposal = proposal;
    this->options = options;
    this->populations = new vector<ParticlePopulation<P> *>();
}

template <class P>
void SMC<P>::run_smc(gsl_rng *random, unsigned int num_particles)
{
    unsigned long R = proposal->num_iterations();

    ParticlePopulation<P> *curr_pop = 0;

    log_marginal_likelihood = 0.0;
    double log_num_particles = log(num_particles);
    for (int r = 0; r < R; r++)
    {
        curr_pop = propose(random, curr_pop, r, num_particles);
        populations->push_back(curr_pop);
        if (!options->track_population && r > 0) {
            delete (*populations)[r-1]; // delete the ParticlePopulation (particles, log_weights, and normalized_weights)
        }
        log_marginal_likelihood += (curr_pop->get_log_norm() - log_num_particles);
        if (curr_pop->get_ess() <= options->essThreshold) {
            // resample
            curr_pop = resample(random, options->resampling, curr_pop, num_particles);
            delete (*populations)[r]; // delete pre-resampling population
            (*populations)[r] = curr_pop;
        }
    }
}

template<class P>
ParticlePopulation<P>* SMC<P>::propose(gsl_rng *random, ParticlePopulation<P> *pop, int iter, int num_proposals)
{
    vector<P> *curr_particles = 0;
    if (pop != 0) {
         curr_particles = pop->get_particles();
        //vector<double> *curr_log_weights = pop->get_log_weights();
    }

    vector<P> *new_particles = new vector<P>(num_proposals);
    vector<double> *new_log_weights = new vector<double>(num_proposals);

    std::pair<int, double> ret;
    for (int n = 0; n < num_proposals; n++)
    {
        if (pop == 0) {
            ret = proposal->propose_initial(random);
        } else {
            ret = proposal->propose_next(random, iter, (*curr_particles)[n]);
        }
        (*new_particles)[n] = ret.first;
        (*new_log_weights)[n] = ret.second;
        //(*new_log_weights)[n] = (*curr_log_weights)[n] + ret.second;
    }
    ParticlePopulation<P> *new_pop = new ParticlePopulation<P>(new_particles, new_log_weights);
    return new_pop;
}

template <typename P>
ParticlePopulation<P>* SMC<P>::resample(const gsl_rng *random, SMCOptions::ResamplingScheme resampling_scheme, ParticlePopulation<P> *pop, int N)
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
        case SMCOptions::ResamplingScheme::RESIDUAL:
            residual_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
    }
    
    vector<P> *curr_particles = pop->get_particles();
    
    vector<P> *particles = new vector<P>();
    vector<double> *log_weights = new vector<double>();
    double log_w = log(N);
    for (int n = 0; n < N; n++)
    {
        particles->push_back((*curr_particles)[indices[n]]);
        log_weights->push_back(log_w);
    }
    
    ParticlePopulation<P> *new_pop = new ParticlePopulation<P>(particles, log_weights);
    
    delete [] indices;
    return new_pop;
}

template <class P>
double SMC<P>::get_log_marginal_likelihood()
{
    return log_marginal_likelihood;
}

template <class P>
SMC<P>::~SMC()
{
}
#endif /* SRC_SMC_HPP_ */
