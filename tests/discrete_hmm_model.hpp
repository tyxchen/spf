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

#include "discrete_hmm_params.hpp"
#include "smc_model.hpp"

using namespace std;

class DiscreteHMM : public ProblemSpecification<int, DiscreteHMMParams>
{
	unsigned long num_states;
	//vector<double> initial_distn;
	//vector<vector<double> > transition_probs;
    //vector<vector<double> > emission_probs;
    vector<int> obs;
public:
	DiscreteHMM(unsigned long num_states, vector<int> obs);
    unsigned long num_iterations();
    std::pair<int, double> propose_initial(gsl_rng *random, DiscreteHMMParams &params);
    std::pair<int, double> propose_next(gsl_rng *random, int t, int curr, DiscreteHMMParams &params);
    ~DiscreteHMM();
};



#endif /* SRC_DISCRETE_HMM_MODEL_HPP_ */
