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
#include "smc_model.hpp"
#include "sampling_utils.hpp"

using namespace std;

template <class P> class SMC
{
    double log_marginal_likelihood = 0;
	vector<P> particles;
	vector<double> log_weights;
    vector<double> *normalized_weights = 0;
	ProblemSpecification<P> *proposal; // pointer to proposal object to be passed into SMC constructor

public:
	SMC(ProblemSpecification<P> *proposal);
	void run_smc(gsl_rng *random, unsigned int num_particles);
    inline vector<P> get_particles() { return particles; }
    inline vector<double> get_log_weights() { return log_weights; }
    vector<double> get_normalized_weights();
    double get_log_marginal_likelihood();
    ~SMC();
};

template <class P>
SMC<P>::SMC(ProblemSpecification<P> *proposal)
{
    this->proposal = proposal;
}

template <class P>
void SMC<P>::run_smc(gsl_rng *random, unsigned int num_particles)
{
    // declare local pointers to be used
    unsigned int *indices = new unsigned int[num_particles];
    vector<P> *temp_particles = new vector<P>(num_particles);
    vector<double> *temp_log_weights = new vector<double>(num_particles);
    normalized_weights = new vector<double>(num_particles);

    unsigned long R = proposal->num_iterations();

    for (int n = 0; n < num_particles; n++)
    {
        pair<P, double> ret = proposal->propose_initial(random);
        particles.push_back(ret.first);
        log_weights.push_back(ret.second);
    }
    log_marginal_likelihood = normalize(log_weights, *normalized_weights) - log(num_particles);

    for (int r = 1; r < R; r++)
    {
        sample_indices(random, num_particles, *normalized_weights, indices);
        
        for (int n = 0; n < num_particles; n++)
        {
            // resample a particle
            P curr = particles[indices[n]];
            pair<P, double> ret = proposal->propose_next(random, r, curr);
            (*temp_particles)[n] = ret.first;
            //(*temp_log_weights)[n] = log_weights[indices[n]] + ret.second;
            (*temp_log_weights)[n] = ret.second;
        }
        // assigns the vector pointed by temp_particles to particles, this copies the values due to operator overloading on =
        particles = *temp_particles;
        log_weights = *temp_log_weights;

        // update the log_marginal_likelihood thus far -- log_weights will be normalized in preparation for the next iteration
        log_marginal_likelihood += (normalize(log_weights, *normalized_weights) - log(num_particles));
    }
        
    // delete local pointers to be used
    delete temp_particles;
    delete temp_log_weights;
    delete [] indices;
}

template <class P>
vector<double> SMC<P>::get_normalized_weights()
{
    return *normalized_weights;
}

template <class P>
double SMC<P>::get_log_marginal_likelihood()
{
    return log_marginal_likelihood;
}

template <class P>
SMC<P>::~SMC()
{
    delete normalized_weights;
}
#endif /* SRC_SMC_HPP_ */
