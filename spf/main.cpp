/*
 * main.cpp
 *
 *  Created on: 6 Mar 2018
 *      Author: seonghwanjun
 */

#include <iostream>
#include <iterator>
#include <math.h>
#include <vector>
#include <gsl/gsl_randist.h>

#include "discrete_hmm_model.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"
#include "smc.hpp"
#include "smc_model.hpp"

using namespace std;

void test_numerical_utils()
{
	// test numerical utils
	double log_w[3];
	log_w[0] = log(2);
	log_w[1] = log(1);
	log_w[2] = log(3);
	cout << "actual: " << log_add(log_w[0], log_w[1]) << " expected: " << log(3) << endl;
	double log_norm = log_add(log_w, 3);
	cout << "actual: " << log_norm << " expected: " << log(6) << endl;
}

int main()
{
	int length = 3;
	unsigned int num_latent_states = 3;
    int num_obs_states = 3;
	vector<double> mu;
	vector<vector<double> > P;
    vector<vector<double> > Q;
    vector<int> latent;
    vector<int> obs;
    
    gsl_rng* random = get_random(121);

	double alpha[num_latent_states];
    alpha[0] = 1; alpha[1] = 1; alpha[2] = 1;
    double beta[num_obs_states];
    beta[0] = 1; beta[1] = 1; beta[2] = 1;

	mu.push_back(1); mu.push_back(0);

    double pi[num_latent_states];
    double eta[num_obs_states];
	for (int i = 0; i < num_latent_states; i++)
	{
		// sample a distribution from Dirichlet distribution
		gsl_ran_dirichlet(random, num_latent_states, alpha, pi);
		vector<double> row1(pi, pi+num_latent_states);
		P.push_back(row1);
        cout << "P_" + to_string(i) << ": \n";
        print_vector(row1);
        cout << endl;
        
        gsl_ran_dirichlet(random, num_obs_states, beta, eta);
        vector<double> row2(eta, eta+num_obs_states);
        Q.push_back(row2);
        cout << "Q_" + to_string(i) << ": \n";
        print_vector(row2);
        cout << endl;
	}
    
    int x_t = 0, y_t = 0;
    for (int t = 0; t < length; t++)
    {
        if (t == 0) {
            x_t = multinomial(random, mu);
        } else {
            x_t = multinomial(random, P[latent[t-1]]);
        }
        y_t = multinomial(random, Q[x_t]);
        latent.push_back(x_t);
        obs.push_back(y_t);
        cout << x_t << "->" << y_t << endl;
    }

	SMC<int> smcAlgorithm(new DiscreteHMM(num_latent_states, mu, P, Q, obs));
    smcAlgorithm.run_smc(random, 1000, true);
    vector<int> particles = smcAlgorithm.get_particles();
    vector<double> log_weights = smcAlgorithm.get_log_weights();
    cout << latent[num_latent_states-1] << endl;
    double sum = 0.0;
    for (unsigned int i = 0; i < particles.size(); i++)
    {
        if (particles[i] == latent[num_latent_states-1])
            sum += log_weights[i];
    }
    cout << "P(x_T = " << latent[num_latent_states - 1] << ") = " << sum << endl;
    
	return 0;
}



