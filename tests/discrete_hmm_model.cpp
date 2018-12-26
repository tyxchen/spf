/*
 * kitagawa_model.cpp
 *
 *  Created on: 30 Mar 2018
 *      Author: seonghwanjun
 */

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "discrete_hmm_model.hpp"

#include "sampling_utils.hpp"
#include "smc_model.hpp"

DiscreteHMM::DiscreteHMM(vector<int> &obs)
{
    this->obs = obs;
	// TODO: check dimensionality of mu and P
}

unsigned long DiscreteHMM::num_iterations()
{
	return obs.size();
}

std::pair<int, double> DiscreteHMM::propose_initial(gsl_rng *random, DiscreteHMMParams &params)
{
	// sample from multinomial distribution parameterized by mu
	int state = multinomial(random, params.initial_distn);
    double logw = log(params.emission_probs[state][obs[0]]);
    return make_pair(state, logw);
}

std::pair<int, double> DiscreteHMM::propose_next(gsl_rng *random, int t, int &curr, DiscreteHMMParams &params)
{
	// sample from multinomial distribution parameterized by mu
	int state = multinomial(random, params.transition_probs[curr]);
    double logw = log(params.emission_probs[state][obs[t]]);
	return make_pair(state, logw);
}

int DiscreteHMM::initial(gsl_rng *random, DiscreteHMMParams &params)
{
    int state = multinomial(random, params.initial_distn);
    return state;
}

int DiscreteHMM::forward(gsl_rng *random, int curr, DiscreteHMMParams &params)
{
    int state = multinomial(random, params.transition_probs[curr]);
    return state;
}

int DiscreteHMM::emission(gsl_rng *random, int curr, DiscreteHMMParams &params)
{
    int state = multinomial(random, params.emission_probs[curr]);
    return state;
}

void DiscreteHMM::generate_data(gsl_rng *random, size_t T, DiscreteHMMParams &params, vector<int> &latent, vector<int> &obs)
{
    // initial proposal
    latent.push_back(initial(random, params));
    obs.push_back(emission(random, latent[0], params));
    for (size_t t = 1; t < T; t++) {
        latent.push_back(forward(random, latent[t-1], params));
        obs.push_back(emission(random, latent[t], params));
    }
}

DiscreteHMM::~DiscreteHMM()
{

}


