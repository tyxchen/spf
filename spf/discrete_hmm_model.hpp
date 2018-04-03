/*
 * kitagawa_model.hpp
 *
 * Defines header file for HMM with discrete latent and observation variables
 *
 *  Created on: 30 Mar 2018
 *      Author: seonghwanjun
 */

#ifndef SRC_DISCRETE_HMM_MODEL_HPP_
#define SRC_DISCRETE_HMM_MODEL_HPP_

#include <utility>
#include <vector>

#include "smc_model.hpp"

using namespace std;

class DiscreteHMM : public ProblemSpecification<int>
{
	unsigned long num_states;
	vector<double> initial_distn;
	vector<vector<double> > transition_probs;
    vector<vector<double> > emission_probs;
    vector<int> obs;
public:
	DiscreteHMM(unsigned long num_states, vector<double> mu, vector<vector<double> > transition_probs, vector<vector<double> > emission_probs, vector<int> obs);
	unsigned long num_iterations();
	std::pair<int, double> propose_initial(gsl_rng *random);
	std::pair<int, double> propose_next(gsl_rng *random, int t, int curr);
    //double compute_log_weight(int t, int x_t, int y_t);
    ~DiscreteHMM();
};



#endif /* SRC_DISCRETE_HMM_MODEL_HPP_ */
