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
#include <cmath>
#include <vector>
#include <unistd.h>
#include <unordered_map>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <gsl/gsl_randist.h>

#include "discrete_hmm_model.hpp"
#include "discrete_hmm_params.hpp"
#include "numerical_utils.hpp"
#include "pmmh.hpp"
#include "sampling_utils.hpp"
#include "smc.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"
#include "spf.hpp"
#include "sv_model.hpp"
#include "sv_model_param_proposal.hpp"

/*****
 * global variables for testing
 *****/
vector<int> x = {2, 0, 2};
vector<int> y = {1, 0, 0};

vector<double> mu({0.1, 0.4, 0.5});
vector<vector<double> > P({{0.2, 0.4, 0.4}, {0.1, 0.8, 0.1}, {0.9, 0.05, 0.05}});
vector<vector<double> > Q({{0.4, 0.1, 0.5}, {0.1, 0.1, 0.8}, {0.3, 0.2, 0.5}});
DiscreteHMMParams true_params(mu, P, Q);
unsigned long num_latent_states = mu.size();
unsigned long num_emission_states = Q[0].size();

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

void test_smc(long seed)
{
    cout << "Testing SMC..." << endl;
    SMCOptions *options = new SMCOptions();
    options->main_seed = seed;
    options->essThreshold = 1;
    options->num_particles = 10000;
    options->resample_last_round = true;

    SMC<int, DiscreteHMMParams> smc(new DiscreteHMM(num_latent_states, y), options);
    smc.run_smc(true_params);
    ParticlePopulation<int> *pop = smc.get_curr_population();
    vector<int> *particles = pop->get_particles();
    vector<double> *normalized_weights = pop->get_normalized_weights();

    // true log marginal likelihood computed using sum-product algorithm
    double true_log_marginal_lik = -4.481802772037807;
    // retrieve the estimate
    double log_marginal_lik = smc.get_log_marginal_likelihood();
    cout << "Estimated log p(y_{1:T}) = " << log_marginal_lik << endl;
    cout << "True log p(y_{1:T}) = " << true_log_marginal_lik << endl;
    cout << "Diff: " << abs(true_log_marginal_lik - log_marginal_lik) << endl;

    // true  marginal likelihood of x_3 = {0, 1, 2} computed using sum-product
    double truth_T[]{0.41014761778484926, 0.16648987890038008, 0.4233625033147706};
    vector<double> probs_T(num_latent_states);
    for (int n = 0; n < options->num_particles; n++)
    {
        probs_T[(*particles)[n]] += (*normalized_weights)[n];
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
    cout << "total error: " << total_error << endl;
}

void test_spf(long seed)
{
    cout << "Testing SPF..." << endl;
    
    gsl_rng *random = generate_random_object(seed);
    SMCOptions *options = new SMCOptions();
    options->main_seed = gsl_rng_get(random);
    options->resampling_seed = gsl_rng_get(random);
    options->essThreshold = DOUBLE_INF;
    options->num_particles = 10000;
    options->max_virtual_particles = 20000;
    options->resample_last_round = true;

    SPF<int, DiscreteHMMParams> spf(new DiscreteHMM(num_latent_states, y), options);
    spf.run_spf(true_params);
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
}

void test_resampling_schemes(long seed, unsigned long T)
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
    
    SMCOptions *options = new SMCOptions();
    options->essThreshold = 1;
    options->num_particles = num_particles;
    options->resample_last_round = false;
    double *vars = new double[resamplingSchemes.size()];
    for (unsigned long i = 0; i < resamplingSchemes.size(); i++) {
        options->resampling_scheme = resamplingSchemes[i];
        SMC<int,DiscreteHMMParams> smc(new DiscreteHMM(num_latent_states, obs), options);
        smc.run_smc(true_params);
        ParticlePopulation<int> *pop = smc.get_curr_population();
        vector<int> *particles = pop->get_particles();
        //vector<double> *log_weights = pop->get_log_weights();
        vector<double> *normalized_weights = pop->get_normalized_weights();
        //cout << normalized_weights->size() << ", " << particles->size() << endl;

        double probs_T[num_latent_states];
        double var_T[num_latent_states];
        for (int n = 0; n < options->num_particles; n++)
        {
            probs_T[(*particles)[n]] += (*normalized_weights)[n];
        }
        for (int n = 0; n < num_particles; n++)
        {
            double diff = (*normalized_weights)[n] - probs_T[(*particles)[n]];
            var_T[(*particles)[n]] += diff * diff;
        }
        vars[i] = var_T[0] / num_particles;
        //cout << options->resampling_scheme << ": " << vars[i] << endl;
    }
    
    return vars;
}

void test_pmcmc(long seed)
{
    // run PMMH on the initial distribution mu
    SMCOptions *smc_options = new SMCOptions();
    smc_options->num_particles = 200;
    SVModel *proposal = new SVModel(sv_y);
    SMC<double, SVModelParams> *smc = new SMC<double, SVModelParams>(proposal, smc_options);
    smc->run_smc(true_sv_params);
    double logZ = smc->get_log_marginal_likelihood();
    cout << logZ << endl; // logZ at true params

    PMCMCOptions *pmcmc_options = new PMCMCOptions(seed, 50000);
    pmcmc_options->burn_in = 5000;
    ParamProposal<SVModelParams> *param_proposal = new SVModelParamProposal();
    ParticleMMH<double, SVModelParams> pmmh(pmcmc_options, smc, param_proposal);
    pmmh.run();
    vector<SVModelParams *> *samples = pmmh.get_parameters();
    // compute the posterior mean for beta
    double mean = 0.0;
    size_t count = 0;
    for (size_t i = pmcmc_options->burn_in; i < samples->size(); i+=10) {
        cout << (*samples)[i]->beta << endl;
        mean += (*samples)[i]->beta;
        count++;
    }
    mean /= count;
    cout << mean << ", " << true_sv_params.beta << endl;
}

int main()
{
    long seed = 123;
    test_smc(seed);
    test_spf(seed);

//    seed = 532;
//    unsigned int T = 100;
//    test_resampling_schemes(seed, T);

    test_pmcmc(seed);
    return 0;
}