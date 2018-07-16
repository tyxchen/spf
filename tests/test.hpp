//
//  test.cpp
//  code containing tests for spf
//
//  Created by Seong-Hwan Jun on 2018-04-02.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <vector>
#include <gsl/gsl_randist.h>

#include "discrete_hmm_model.hpp"
#include "discrete_hmm_params.hpp"
#include "lgm_param_proposal.hpp"
#include "linear_gaussian_model.hpp"
#include "linear_gaussian_model_params.hpp"
#include "numerical_utils.hpp"
#include "pmcmc.hpp"
#include "sampling_utils.hpp"
#include "smc.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"
#include "spf.hpp"

void test_numerical_utils();
double *compare_resampling_schemes(long seed, int num_particles, vector<SMCOptions::ResamplingScheme> resamplingSchemes);
void test_smc();
void test_resampling_schemes();

void test_numerical_utils()
{
    // test log_add
    double log_w[3];
    log_w[0] = log(2);
    log_w[1] = log(1);
    log_w[2] = log(3);
    cout << "actual: " << log_add(log_w[0], log_w[1]) << " expected: " << log(3) << endl;
    double log_norm = log_add(log_w, 3);
    cout << "actual: " << log_norm << " expected: " << log(6) << endl;
}

void test_smc()
{
    long seed = 20180402;
    int length = 3;
    vector<double> mu({0.1, 0.4, 0.5});
    vector<vector<double> > P({{0.2, 0.4, 0.4}, {0.1, 0.8, 0.1}, {0.9, 0.05, 0.05}});
    vector<vector<double> > Q({{0.4, 0.1, 0.5}, {0.1, 0.1, 0.8}, {0.3, 0.2, 0.5}});
    
    unsigned long num_latent_states = mu.size();
    //unsigned long num_obs_states = Q[0].size();

    vector<int> latent;
    vector<int> obs;

    gsl_rng* random = generate_random_object(seed);
    
    /* uncomment to sample the transition and emission probs
    
    double alpha[num_latent_states];
    alpha[0] = 1; alpha[1] = 1; alpha[2] = 1;
    double beta[num_obs_states];
    beta[0] = 1; beta[1] = 1; beta[2] = 1;
    
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
     */
    
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

    SMCOptions *options = new SMCOptions();
    options->essThreshold = 1;
    options->resampling_scheme = SMCOptions::ResamplingScheme::STRATIFIED;
    options->num_particles = 10000;
    options->resampling_random = generate_random_object(1);
    options->main_random = generate_random_object(1);
    DiscreteHMMParams params(mu, P, Q);
    SMC<int, DiscreteHMMParams> smc(new DiscreteHMM(num_latent_states, obs), options);
    smc.run_smc(params);
    ParticlePopulation<int> *pop = smc.get_curr_population();
    vector<int> *particles = pop->get_particles();
    //vector<double> *log_weights = pop->get_log_weights();
    vector<double> *normalized_weights = pop->get_normalized_weights();

    // get the estimate of the marginal likelihood: p(y_{1:T})
    double true_log_marginal_lik = -4.481802772037807;
    double log_marginal_lik = smc.get_log_marginal_likelihood();
    cout << "Estimated log p(y_{1:T}) = " << log_marginal_lik << endl;
    cout << "True log p(y_{1:T}) = " << true_log_marginal_lik << endl;
    cout << "Diff: " << abs(true_log_marginal_lik - log_marginal_lik) << endl;
    
    vector<double> probs_T(num_latent_states);
    for (int n = 0; n < options->num_particles; n++)
    {
        probs_T[(*particles)[n]] += (*normalized_weights)[n];
    }
    
    double truth_T[]{0.41014761778484926, 0.16648987890038008, 0.4233625033147706};
    double total_error = 0.0;
    for (int i = 0; i < num_latent_states; i++)
    {
        double err = abs(probs_T[i] - truth_T[i]);
        total_error += err;
        cout << "Estimated P(x_T = " << i << ") = " << probs_T[i] << endl;
        cout << "True P(x_T = " << i << ") = " << truth_T[i] << endl;
        cout << "Diff: " << err << endl;
    }
    cout << "total error: " << total_error << endl;
    
    // TODO: put this into a unit test framework
}

void test_resampling_schemes()
{
    vector<SMCOptions::ResamplingScheme> resampling_schemes;
    resampling_schemes.push_back(SMCOptions::ResamplingScheme::MULTINOMIAL);
    resampling_schemes.push_back(SMCOptions::ResamplingScheme::SYSTEMATIC);
    resampling_schemes.push_back(SMCOptions::ResamplingScheme::STRATIFIED);
    double vars[3];
    int num_runs = 20;
    for (int i = 0; i < num_runs; i++) {
        int seed = rand();
        double *ret = compare_resampling_schemes(seed, 10000, resampling_schemes);
        vars[0] += ret[0];
        vars[1] += ret[1];
        vars[2] += ret[2];
    }
    cout << vars[0] << ", " << vars[1] << ", " << vars[2] << endl;
}

void test_spf()
{
    long seed = 20180402;
    int length = 3;
    vector<double> mu({0.1, 0.4, 0.5});
    vector<vector<double> > P({{0.2, 0.4, 0.4}, {0.1, 0.8, 0.1}, {0.9, 0.05, 0.05}});
    vector<vector<double> > Q({{0.4, 0.1, 0.5}, {0.1, 0.1, 0.8}, {0.3, 0.2, 0.5}});

    unsigned long num_latent_states = mu.size();
    //unsigned long num_obs_states = Q[0].size();

    vector<int> latent;
    vector<int> obs;
    
    gsl_rng* random = generate_random_object(seed);
    
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

    SMCOptions *options = new SMCOptions();
    options->num_particles = 10000;
    options->max_virtual_particles = 20000;
    options->essThreshold = DOUBLE_INF;
    options->resampling_random = generate_random_object(seed);
    options->main_random = generate_random_object(1);
    options->resampling_scheme = SMCOptions::ResamplingScheme::STRATIFIED;
    DiscreteHMMParams params(mu, P, Q);
    SPF<int, DiscreteHMMParams> spf(new DiscreteHMM(num_latent_states, obs), options);
    spf.run_spf(params);
    ParticlePopulation<int> *pop = spf.get_curr_population();
    vector<int> *particles = pop->get_particles();

    // get the estimate of the marginal likelihood: p(y_{1:T})
    double true_log_marginal_lik = -4.481802772037807;
    double log_marginal_lik = spf.get_log_marginal_likelihood();
    cout << "Estimated log p(y_{1:T}) = " << log_marginal_lik << endl;
    cout << "True log p(y_{1:T}) = " << true_log_marginal_lik << endl;
    cout << "Diff: " << abs(true_log_marginal_lik - log_marginal_lik) << endl;
    
    double probs_T[num_latent_states];
    for (int n = 0; n < options->num_particles; n++)
    {
        probs_T[(*particles)[n]]++;
    }
    
    double truth_T[]{0.41014761778484926, 0.16648987890038008, 0.4233625033147706};
    double total_error = 0.0;
    for (int i = 0; i < num_latent_states; i++)
    {
        probs_T[i] /= options->num_particles;
        double err = abs(probs_T[i] - truth_T[i]);
        total_error += err;
        cout << "Estimated P(x_T = " << i << ") = " << probs_T[i] << endl;
        cout << "True P(x_T = " << i << ") = " << truth_T[i] << endl;
        cout << "Diff: " << err << endl;
    }
    cout << "total error: " << total_error << endl;
    
    // TODO: put this into a unit test framework
}

double *compare_resampling_schemes(long seed, int num_particles, vector<SMCOptions::ResamplingScheme> resamplingSchemes)
{
    int length = 3;
    vector<double> mu({0.1, 0.4, 0.5});
    vector<vector<double> > P({{0.2, 0.4, 0.4}, {0.1, 0.8, 0.1}, {0.9, 0.05, 0.05}});
    vector<vector<double> > Q({{0.4, 0.1, 0.5}, {0.1, 0.1, 0.8}, {0.3, 0.2, 0.5}});
    
    unsigned long num_latent_states = mu.size();
    //unsigned long num_obs_states = Q[0].size();
    
    vector<int> latent;
    vector<int> obs;
    
    gsl_rng* random = generate_random_object(seed);
    
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
    
    DiscreteHMMParams params(mu, P, Q);
    SMCOptions *options = new SMCOptions();
    options->essThreshold = 1;
    double *vars = new double[resamplingSchemes.size()];
    for (unsigned long i = 0; i < resamplingSchemes.size(); i++) {
        options->resampling_scheme = resamplingSchemes[i];
        SMC<int,DiscreteHMMParams> smc(new DiscreteHMM(num_latent_states, obs), options);
        smc.run_smc(params);
        ParticlePopulation<int> *pop = smc.get_curr_population();
        vector<int> *particles = pop->get_particles();
        //vector<double> *log_weights = pop->get_log_weights();
        vector<double> *normalized_weights = pop->get_normalized_weights();
        
        double probs_T[num_latent_states];
        double var_T[num_latent_states];
        for (int n = 0; n < num_particles; n++)
        {
            probs_T[(*particles)[n]] += (*normalized_weights)[n];
        }
        for (int n = 0; n < num_particles; n++)
        {
            double diff = (*normalized_weights)[n] - probs_T[(*particles)[n]];
            var_T[(*particles)[n]] += diff * diff;
        }
        vars[i] = var_T[0] / num_particles;
        cout << options->resampling_scheme << ": " << vars[i] << endl;
    }
    
    return vars;
}

void test_pmcmc()
{
    long seed = 20180402;
    int length = 3;

    vector<double> latent;
    vector<double> obs;

    gsl_rng* random = generate_random_object(seed);

    // generate x_t, y_t where x_t
    // x_1 ~ Normal(0, nu^2)
    // x_t | x_{t-1} ~ Normal(A*x_{t-1} + a, sigma^2)
    // y_t | x_t ~ Normal(B*x_t + b, tau^2)
    
    double mu_A = 0.2; double sd_A = 0.1; // prior over A
    double mu_B = -0.7; double sd_B = 1.1; // prior over B
    double mu_a = 0; double sd_a = 1.4; // prior over a
    double mu_b = 0; double sd_b = 0.7; // prior over b
    double alpha_x0 = 0.51; double beta_x0 = 2.01; // prior over x0 (InvGamma or Gamma)
    double alpha_x = 1.01; double beta_x = 1.01; // prior over x (InvGamma or Gamma)
    double alpha_y = 1.1; double beta_y = 1.1; // prior over y (InvGamma or Gamma)

    double A = gsl_ran_gaussian(random, sd_A) + mu_A;
    double B = gsl_ran_gaussian(random, sd_B) + mu_B;
    double a = gsl_ran_gaussian(random, sd_a) + mu_a;
    double b = gsl_ran_gaussian(random, sd_b) + mu_b;
    double nu = gsl_ran_gamma(random, alpha_x0, beta_x0);
    double sigma = gsl_ran_gamma(random, alpha_x, beta_x);
    double tau = gsl_ran_gamma(random, alpha_y, beta_y);

    double x_t = 0, y_t = 0;
    for (int t = 0; t < length; t++)
    {
        if (t == 0) {
            x_t = gsl_ran_gaussian(random, nu);
        } else {
            x_t = gsl_ran_gaussian(random, sigma) + A * latent[t-1] + a;
        }
        y_t = gsl_ran_gaussian(random, tau) + B * x_t + b;
        latent.push_back(x_t);
        obs.push_back(y_t);
        cout << x_t << "->" << y_t << endl;
    }

    SMCOptions *options = new SMCOptions();
    options->essThreshold = 1;
    options->debug = false;
    options->resampling_scheme = SMCOptions::ResamplingScheme::STRATIFIED;
    options->num_particles = 1000;
    options->resampling_random = generate_random_object(1);
    options->main_random = generate_random_object(1);

    SMC<double, LinearGaussianModelParams> *smc = new SMC<double, LinearGaussianModelParams>(new LinearGaussianModel(&obs), options);
    
    PMCMCOptions *pmcmc_options = new PMCMCOptions();
    pmcmc_options->num_iterations = 10000;
    pmcmc_options->random = generate_random_object(seed * 2);

    ParamProposal<LinearGaussianModelParams> *param_proposal = new LGMParamProposal(mu_A, sd_A,
                                                                                    mu_B, sd_B,
                                                                                    mu_a, sd_a,
                                                                                    mu_b, sd_b,
                                                                                    alpha_x0, beta_x0,
                                                                                    alpha_x, beta_x,
                                                                                    alpha_y, beta_y);
    ParticleMMH<double, LinearGaussianModelParams> pmcmc(pmcmc_options, smc, param_proposal);
    pmcmc.run();
    
    cout << "Truth" << endl;
    cout << "A: " << std::to_string(A) << endl;
    cout << "a: " << std::to_string(a) << endl;
    cout << "B: " << std::to_string(B) << endl;
    cout << "b: " << std::to_string(b) << endl;
    cout << "nu: " << std::to_string(nu) << endl;
    cout << "tau: " << std::to_string(tau) << endl;
    cout << "sigma: " << std::to_string(sigma) << endl;
    
    // output the chain to a file for further analysis
    vector<LinearGaussianModelParams*> *params = pmcmc.get_parameters();
    
    char * dir = getcwd(NULL, 0);
    printf("Current dir: %s", dir);
    
    ofstream out("/Users/seonghwanjun/Google Drive/Research/smc-research/repos/spf/output/test_pmcmc_output.csv", ofstream::out);
    for (size_t i = 0; i < params->size(); i++)
    {
        LinearGaussianModelParams *p = (*params)[i];
        out << p->A << ", " << p->a << ", " << p->B << ", " << p->b << ", " << p->nu << ", " << p->tau << ", " << p->sigma << endl;
    }
    out.close();
}
