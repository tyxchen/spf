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

#include "smc_model.hpp"
#include "sampling_utils.hpp"

using namespace std;

template <class T> class SMC
{
	vector<T> particles;
	vector<double> log_weights;
	ProblemSpecification<T> *proposal; // pointer to proposal object to be passed into SMC constructor

public:
	SMC(ProblemSpecification<T> *proposal);
	void run_smc(gsl_rng *random, unsigned int num_particles, bool final_resampling);
    inline vector<T> get_particles() { return particles; }
    inline vector<double> get_log_weights() { return log_weights; }
};

template <class P>
SMC<P>::SMC(ProblemSpecification<P> *proposal)
{
    this->proposal = proposal;
}

template <class P>
void SMC<P>::run_smc(gsl_rng *random, unsigned int num_particles, bool final_resampling)
{
    // begin SMC code
    unsigned long R = proposal->num_iterations();
    unsigned int indices[num_particles];
    
    for (int n = 0; n < num_particles; n++)
    {
        pair<P, double> ret = proposal->propose_initial(random);
        particles.push_back(ret.first);
        log_weights.push_back(ret.second);
    }
    
    for (int r = 1; r < R; r++)
    {
        normalize_destructively(log_weights);
        multinomial(random, num_particles, log_weights, indices);
        
        for (int n = 0; n < num_particles; n++)
        {
            // resample a particle
            P curr = particles[indices[n]];
            pair<P, double> ret = proposal->propose_next(random, r, curr);
            particles[n] = ret.first;
            log_weights[n] = ret.second;
        }
    }
    
    // resample last round
    if (final_resampling)
    {
        normalize_destructively(log_weights);
    }
}

#endif /* SRC_SMC_HPP_ */
