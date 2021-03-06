//
//  test.cpp
//  code containing tests for spf
//
//  Created by Seong-Hwan Jun on 2018-04-02.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <vector>
#include <unistd.h>
#include <unordered_map>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <gsl/gsl_randist.h>

#include "csmc.hpp"
#include "discrete_hmm_model.hpp"
#include "discrete_hmm_params.hpp"
#include "normal_normal_model.hpp"
#include "numerical_utils.hpp"
#include "pg.hpp"
#include "pmmh.hpp"
#include "sampling_utils.hpp"
#include "smc.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"
#include "spf.hpp"
#include "sv_model.hpp"
#include "sv_model_gibbs_proposal.hpp"
#include "sv_model_random_walk_proposal.hpp"

/*****
 * global variables for testing
 *****/
double ERR_TOL  = 1e-2;

vector<int> x = {2, 0, 2};
vector<int> y = {1, 0, 0};

vector<double> mu({0.1, 0.4, 0.5});
vector<vector<double> > P({{0.2, 0.4, 0.4}, {0.1, 0.8, 0.1}, {0.9, 0.05, 0.05}});
vector<vector<double> > Q({{0.4, 0.1, 0.5}, {0.1, 0.1, 0.8}, {0.3, 0.2, 0.5}});
DiscreteHMMParams true_params(mu, P, Q);
unsigned long num_latent_states = mu.size();
unsigned long num_emission_states = Q[0].size();

// exact log marginal likelihood computed using sum-product algorithm (outside of this library)
double true_log_marginal_lik = -4.481802772037807;
// exact marginal likelihood of x_3 = {0, 1, 2} computed using sum-product
double truth_T[]{0.41014761778484926, 0.16648987890038008, 0.4233625033147706};

/*****
 * data for SV model
 *****/
SVModelParams true_sv_params(1.0, 0.16, 0.64);

vector<double> sv_x = {-0.089676103448354,
    0.159717226815506,
    0.180403264441257,
    0.25414985739953,
    0.144253401096565,
    0.340106488686879,
    0.404229920781928,
    0.315295339221276,
    0.394951415737954,
    0.507168359988144};

vector<double> sv_y = {-0.0861939850317827,
    0.0338816900378329,
    0.841373193537888,
    -0.668109083552618,
    -0.210869602735491,
    0.207082866424709,
    0.067919515518584,
    1.00321592058082,
    -1.1956525434525,
    -0.321580675697379};

/*****
 * declare helper functions in advance
 *****/
double *compare_resampling_schemes(long seed, unsigned int num_particles, unsigned int length, vector<SMCOptions::ResamplingScheme> resamplingSchemes);

/*****
 * test functions
 *****/
void test_gsl() {
    gsl_rng* random = generate_random_object(123);
    gsl_rng* random_clone = gsl_rng_clone(random);
    double *a = new double[5];
    for (int i = 0; i < 5; i++) {
        double u = gsl_rng_uniform(random);
        a[i] = u;
        cout << u << endl;
    }
    for (int i = 0; i < 5; i++) {
        double u = gsl_rng_uniform(random);
        double v = gsl_rng_uniform(random_clone);
        cout << u << endl;
        cout << a[i] << " - " << v << " = " << (a[i] - v) << endl;
    }
    
    gsl_ran_shuffle(random, a, 5, sizeof(double));
    for (int i = 0; i < 5; i++) {
        cout << a[i] << endl;
    }
    
    int b = 2;
    int idx = b++ % 5;
    cout << b << " " << idx << endl;
}

void test_copy_unordered_map()
{
    unordered_map<string, int> map;
    map = {{"A", 1}, {"B", 2}};
    unordered_map<string, int> copy = map;
    copy["A"] = 0;
    for (auto &elem : copy) {
        cout << elem.first << ", " << elem.second << endl;
    }
    for (auto &elem : map) {
        cout << elem.first << ", " << elem.second << endl;
    }
}

void test_boost()
{
    string s = "a,b, c ,,e,f,";
    vector <string> fields;
    
    cout << "Original = \"" << s << "\"\n\n";
    
    cout << "Split on \',\' only\n";
    boost::split( fields, s, boost::is_any_of( "," ) );
    cout << fields[0] << endl;
}

bool test_numerical_utils()
{
    // test log_add
    double log_w[3];
    double log_sum;
    double log_norm;
    log_w[0] = log(2);
    log_w[1] = log(1);
    log_w[2] = log(3);
    log_sum = log_add(log_w[0], log_w[1]);
    cout << "actual: " << log_sum << " expected: " << log(3) << endl;
    if (abs(log_sum - log(3)) > ERR_TOL) {
        return false;
    }
    
    log_norm = log_add(log_w, 3);
    cout << "actual: " << log_norm << " expected: " << log(6) << endl;
    if (abs(log_norm - log(6)) > ERR_TOL) {
        return false;
    }
    return true;
}

void test_smc(long seed)
{
    cout << "=====Testing SMC=====" << endl;
    SMCOptions options;
    options.main_seed = seed;
    options.ess_threshold = 1;
    options.num_particles = 10000;
    options.resample_last_round = true;

    DiscreteHMM model(y);
    SMC<int, DiscreteHMMParams> smc(model, options);
    smc.run_smc(true_params);
    ParticlePopulation<int> *pop = smc.get_curr_population();
    vector<shared_ptr<int>> *particles = pop->get_particles();
    vector<double> *normalized_weights = pop->get_normalized_weights();

    // retrieve the estimate
    double log_marginal_lik = smc.get_log_marginal_likelihood();
    cout << "Estimate log P(y)= " << log_marginal_lik << "; Expected log P(y)= " << true_log_marginal_lik << endl;
    double diff = abs(exp(true_log_marginal_lik) - exp(log_marginal_lik));
    cout << "Diff: " << diff << endl;
    assert(diff < ERR_TOL);

    vector<double> probs_T(num_latent_states);
    for (int n = 0; n < options.num_particles; n++)
    {
        int state = *particles->at(n).get();
        probs_T[state] += (*normalized_weights)[n];
    }
    
    double total_error = 0.0;
    for (int i = 0; i < num_latent_states; i++)
    {
        double err = abs(probs_T[i] - truth_T[i]);
        total_error += err;
        cout << "Estimated P(x_T = " << i << ") = " << probs_T[i] << endl;
        cout << "True P(x_T = " << i << ") = " << truth_T[i] << endl;
        cout << "Diff: " << err << endl;
    }
    cout << "average error: " << total_error/num_latent_states << endl;
    assert(total_error/num_latent_states < ERR_TOL);
}

bool test_csmc(long seed, size_t num_threads)
{
    cout << "=====Testing conditional SMC=====" << endl;
    SMCOptions options;
    options.main_seed = seed;
    options.resampling_seed = seed+123;
    options.ess_threshold = 0.5;
    options.num_particles = 100000;
    options.num_threads = num_threads;

    DiscreteHMM model(y);
    ConditionalSMC<int, DiscreteHMMParams> csmc(model, options);
    shared_ptr<ParticleGenealogy<int>> genealogy = csmc.initialize(true_params);
    double logZ = csmc.get_log_marginal_likelihood();

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t num_iter = 0; num_iter < 20; num_iter++) {
        genealogy = csmc.run_csmc(true_params, genealogy);
        logZ = csmc.get_log_marginal_likelihood();
        //cout << logZ << endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << elapsed.count() << " seconds with " << num_threads << " threads." << endl;

    unsigned int r = 2;
    vector<int> particles;
//    vector<double> *log_weights = csmc.get_log_weights(r);
    vector<double> normalized_weights;
    double log_norm = csmc.get_log_norm(r);
    double sum = 0.0;
    double norm_w = 0.0;
    for (size_t k = 0; k < options.num_particles; k++) {
        const int &state  = csmc.get_state(r, k);
        double log_w = csmc.get_log_weight(r, k);
        norm_w = exp(log_w - log_norm);
        particles.push_back(state);
        normalized_weights.push_back(norm_w);
        sum += norm_w;
    }

    // retrieve the estimate
    double log_marginal_lik = csmc.get_log_marginal_likelihood();
    cout << "Estimate log P(y)= " << log_marginal_lik << "; Expected log P(y)= " << true_log_marginal_lik << endl;
    double diff = abs(exp(true_log_marginal_lik) - exp(log_marginal_lik));
    cout << "Diff: " << diff << endl;
    if (diff > ERR_TOL) {
        return false;
    }

    vector<double> probs_T(num_latent_states);
    for (int n = 0; n < options.num_particles; n++)
    {
        probs_T[particles[n]] += normalized_weights[n];
    }
    
    double total_error = 0.0;
    for (int i = 0; i < num_latent_states; i++)
    {
        double err = abs(probs_T[i] - truth_T[i]);
        total_error += err;
        cout << "Estimated P(x_T = " << i << ") = " << probs_T[i] << endl;
        cout << "True P(x_T = " << i << ") = " << truth_T[i] << endl;
        cout << "Diff: " << err << endl;
    }
    cout << "average error: " << total_error/num_latent_states << endl;

    if (total_error/num_latent_states > ERR_TOL) {
        return false;
    }

    return true;
}

void test_spf(long seed)
{
    cout << "=====Testing SPF=====" << endl;

    gsl_rng *random = generate_random_object(seed);
    SMCOptions options;
    options.main_seed = gsl_rng_get(random);
    options.resampling_seed = gsl_rng_get(random);
    options.ess_threshold = DOUBLE_INF;
    options.num_particles = 10000;
    options.max_virtual_particles = 20000;
    options.resample_last_round = true;

    DiscreteHMM model(y);
    SPF<int, DiscreteHMMParams> spf(model, options);
    spf.run_spf(true_params);
    ParticlePopulation<int> *pop = spf.get_curr_population();
    vector<shared_ptr<int>> *particles = pop->get_particles();

    double log_marginal_lik = spf.get_log_marginal_likelihood();
    cout << "Estimate log P(y)= " << log_marginal_lik << "; Expected log P(y)= " << true_log_marginal_lik << endl;
    double diff = abs(exp(true_log_marginal_lik) - exp(log_marginal_lik));
    cout << "Diff: " << diff << endl;
    assert(diff < ERR_TOL);

    double probs_T[num_latent_states];
    for (int n = 0; n < options.num_particles; n++)
    {
        int state = *particles->at(n).get();
        probs_T[state]++;
    }
    
    double total_error = 0.0;
    for (int i = 0; i < num_latent_states; i++)
    {
        probs_T[i] /= options.num_particles;
        double err = abs(probs_T[i] - truth_T[i]);
        total_error += err;
        cout << "Estimated P(x_T = " << i << ") = " << probs_T[i] << endl;
        cout << "True P(x_T = " << i << ") = " << truth_T[i] << endl;
        cout << "Diff: " << err << endl;
    }
    cout << "average error: " << total_error/num_latent_states << endl;
    assert(total_error/num_latent_states < ERR_TOL);
}

void test_smc_sampler(long seed)
{
    gsl_rng *random = generate_random_object(seed+3);

    cout << "=====Testing SMC Sampler=====" << endl;
    
    // generate data and compute exact marginal likelihood
    double mu_0 = 2.1;
    double sigma_0 = 1.3;
    double true_mu = gsl_ran_gaussian(random, sigma_0) + mu_0;
    cout << "True mu: " << true_mu << endl;
    double sigma = 1.1;

    // generate data
    size_t num_data = 2;
    vector<double> data;
    double data_log_lik = 0.0;
    double sum_data = 0.0;
    for (size_t n = 0; n < num_data; n++) {
        data.push_back(gsl_ran_gaussian(random, sigma) + true_mu);
        data_log_lik += log(gsl_ran_gaussian_pdf(data[n] - true_mu, sigma));
        sum_data += data[n];
    }
    double log_prior = log(gsl_ran_gaussian_pdf(true_mu - mu_0, sigma_0));

    double var_0 = pow(sigma_0, 2);
    double var = pow(sigma, 2);
    double posterior_mu = 1/(1/var_0 + num_data/var);
    posterior_mu *= ((mu_0/var_0) + (sum_data/var));
    double posterior_var = 1./((1/var_0) + (num_data/var));
    double log_posterior = log(gsl_ran_gaussian_pdf(true_mu - posterior_mu, sqrt(posterior_var)));

    double exact_log_marginal = data_log_lik + log_prior - log_posterior;
    cout << "Exact log marginal: " << exact_log_marginal << endl;

    // run SMC and estimate the marginal likelihood
    SMCOptions options;
    options.main_seed = gsl_rng_get(random);
    options.ess_threshold = 1.0;
    options.num_particles = 10000;
    options.resample_last_round = true;

    NormalNormalHyperParams hp(mu_0, sigma_0);

    // if num_iter = t, it will run for (t+1) iterations, the target for the first iteration will be the prior
    // the target for the second iteration and on will be tempered likelihod * prior
    size_t num_iter = 20;
    NormalNormalModel model(num_iter, 8, data, sigma);
    SMC<NormalNormalState, NormalNormalHyperParams> smc(model, options);
    smc.run_smc(hp);
    double smc_log_marginal = smc.get_log_marginal_likelihood();
    cout << "SMC sampler estimate: " << smc_log_marginal << endl;
    assert(abs(exp(exact_log_marginal) - exp(smc_log_marginal)) < ERR_TOL);
}

void test_resampling_schemes(long seed, unsigned int T)
{
    vector<SMCOptions::ResamplingScheme> resampling_schemes;
    resampling_schemes.push_back(SMCOptions::ResamplingScheme::MULTINOMIAL);
    resampling_schemes.push_back(SMCOptions::ResamplingScheme::SYSTEMATIC);
    resampling_schemes.push_back(SMCOptions::ResamplingScheme::STRATIFIED);
    double vars[3];
    int num_runs = 20;
    gsl_rng *random = generate_random_object(seed);
    for (int i = 0; i < num_runs; i++) {
        double *ret = compare_resampling_schemes(gsl_rng_get(random), 10000, T, resampling_schemes);
        vars[0] += ret[0];
        vars[1] += ret[1];
        vars[2] += ret[2];
    }
    cout << "Multinomial variance: " << vars[0] << endl;
    cout << "Systematic variance: " << vars[1] << endl;
    cout << "Stratified variance: " << vars[2] << endl;
}

double *compare_resampling_schemes(long seed, unsigned int num_particles, unsigned int length, vector<SMCOptions::ResamplingScheme> resamplingSchemes)
{
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
        //cout << x_t << "->" << y_t << endl;
    }
    
    SMCOptions options;
    options.ess_threshold = 1;
    options.num_particles = num_particles;
    options.resample_last_round = false;
    DiscreteHMM model(obs);
    double *vars = new double[resamplingSchemes.size()];
    for (unsigned long i = 0; i < resamplingSchemes.size(); i++) {
        options.resampling_scheme = resamplingSchemes[i];
        SMC<int,DiscreteHMMParams> smc(model, options);
        smc.run_smc(true_params);
        ParticlePopulation<int> *pop = smc.get_curr_population();
        vector<shared_ptr<int>> *particles = pop->get_particles();
        //vector<double> *log_weights = pop->get_log_weights();
        vector<double> *normalized_weights = pop->get_normalized_weights();
        //cout << normalized_weights->size() << ", " << particles->size() << endl;

        double probs_T[num_latent_states];
        double var_T[num_latent_states];
        for (int n = 0; n < options.num_particles; n++)
        {
            int state = *particles->at(n).get();
            probs_T[state] += (*normalized_weights)[n];
        }
        for (int n = 0; n < num_particles; n++)
        {
            int state = *particles->at(n).get();
            double diff = (*normalized_weights)[n] - probs_T[state];
            var_T[state] += diff * diff;
        }
        vars[i] = var_T[0] / num_particles;
        //cout << options->resampling_scheme << ": " << vars[i] << endl;
    }
    
    return vars;
}

void test_pmmh(long seed)
{
    // run PMMH on beta of SVModel
    SMCOptions smc_options;
    smc_options.num_particles = 200;
    SVModel proposal(sv_y);
    ConditionalSMC<double, SVModelParams> smc(proposal, smc_options);
    smc.initialize(true_sv_params);
    double logZ = smc.get_log_marginal_likelihood();
    cout << logZ << endl; // logZ at true params

    PMCMCOptions pmcmc_options(seed, 10000);
    pmcmc_options.burn_in = 1000;
    SVModelRandomWalkProposal param_proposal;
    ParticleMMH<double, SVModelParams> pmmh(pmcmc_options, smc, param_proposal);
    pmmh.run();
    vector<SVModelParams *> *samples = pmmh.get_parameters();
    // compute the posterior mean for beta
    double mean = 0.0;
    size_t count = 0;
    for (size_t i = pmcmc_options.burn_in; i < samples->size(); i+=10) {
        mean += (*samples)[i]->beta;
        count++;
    }
    mean /= count;
    cout << mean << ", " << true_sv_params.beta << endl;
}

void test_pg(long seed)
{
    // run PG on beta of SVModel
    SMCOptions smc_options;
    smc_options.num_particles = 200;
    smc_options.num_threads = 4;
    SVModel proposal(sv_y);
    ConditionalSMC<double, SVModelParams> csmc(proposal, smc_options);

    PMCMCOptions pmcmc_options(seed, 4000);
    pmcmc_options.burn_in = 1000;
    
    double a = 0.01, b = 0.01;
    SVModelGibbsProposal param_proposal(a, b, sv_y);
    
    ParticleGibbs<double, SVModelParams> pg(pmcmc_options, csmc, param_proposal);
    pg.run();
    vector<shared_ptr<SVModelParams> > &samples = pg.get_parameters();
    // compute the posterior mean for beta
    double mean = 0.0;
    size_t count = 0;
    for (size_t i = pmcmc_options.burn_in; i < samples.size(); i+=10) {
        mean += samples[i]->beta;
        count++;
    }
    mean /= count;
    cout << mean << ", " << true_sv_params.beta << endl;
}

void test_open_mp()
{
    omp_set_num_threads(8);
#pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf("%i\n", ID);
    }
}

int main()
{
    test_open_mp();

    long seed = 123;
    test_smc(seed);
    test_spf(seed);
    test_csmc(seed, 4);
    test_csmc(seed, 1);
    test_smc_sampler(seed);

    // implement formal validation framework for PMMH and PG
    //test_pmmh(seed);
    test_pg(seed);
    cout << "All tests passed!" << endl;
    return 0;
}
