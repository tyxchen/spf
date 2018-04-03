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

DiscreteHMM::DiscreteHMM(unsigned long num_states, vector<double> mu, vector<vector<double> > transition_probs, vector<vector<double> > emission_probs, vector<int> obs)
{
    this->obs = obs;
	this->num_states = num_states;
	this->initial_distn = mu;
	this->transition_probs = transition_probs;
    this->emission_probs = emission_probs;

	// TODO: check dimensionality of mu and P
}

unsigned long DiscreteHMM::num_iterations()
{
	return obs.size();
}

std::pair<int, double> DiscreteHMM::propose_initial(gsl_rng *random)
{
	// sample from multinomial distribution parameterized by mu
	int state = multinomial(random, initial_distn);
    return make_pair(state, log(emission_probs[state][obs[0]]));
}

std::pair<int, double> DiscreteHMM::propose_next(gsl_rng *random, int t, int curr)
{
	// sample from multinomial distribution parameterized by mu
	int state = multinomial(random, transition_probs[curr]);
    double logw = log(emission_probs[state][obs[t]]);
	return make_pair(state, logw);
}

DiscreteHMM::~DiscreteHMM()
{

}


