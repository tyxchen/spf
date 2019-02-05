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
#include "particle.hpp"
#include "smc_model.hpp"

using namespace std;

class DiscreteHMM : public ProblemSpecification<int, DiscreteHMMParams>
{
    vector<int> obs;
public:
	DiscreteHMM(vector<int> &obs);
    unsigned long num_iterations();
    int *propose_initial(gsl_rng *random, double &log_w, DiscreteHMMParams &params);
    int *propose_next(gsl_rng *random, int t, const int &curr, double &log_w, DiscreteHMMParams &params);

    static int initial(gsl_rng *random, DiscreteHMMParams &params);
    static int forward(gsl_rng *random, int curr, DiscreteHMMParams &params);
    static int emission(gsl_rng *random, int curr, DiscreteHMMParams &params);
    static void generate_data(gsl_rng *random, size_t T, DiscreteHMMParams &params, vector<int> &latent, vector<int> &obs);
    ~DiscreteHMM();
};



#endif /* SRC_DISCRETE_HMM_MODEL_HPP_ */
