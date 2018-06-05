/*
 * main.cpp
 *
 *  Created on: 6 Mar 2018
 *      Author: seonghwanjun
 */

#include <unistd.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <math.h>
#include <unordered_map>
#include <vector>
#include <string>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <gsl/gsl_randist.h>

#include "discrete_hmm_model.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"
#include "somatic_mutation.hpp"
#include "smc.hpp"
#include "smc_model.hpp"
#include "test.hpp"
#include "cancer_phylo_param.hpp"
#include "tssb_state.hpp"
#include "tssb_problem_spec.hpp"

using namespace std;

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

void test_find_node(gsl_rng *random, double u, CancerPhyloParameters params)
{
    unordered_map<string, double> node2nu_stick;
    string node_str = "0";
    double nu = 0.0;
    while (true) {
        vector<double> psi_sticks;
        if (!node2nu_stick.count(node_str)) {
            node2nu_stick[node_str] = beta(random, 1, params.alpha_0);
        }
        nu = node2nu_stick[node_str];
        cout << "u: " << u << " nu: " << nu << endl;
        if (u < nu) {
            break;
        }
        u = (u - nu) / (1 - nu);
        cout << "adjusted u: " << u << endl;

        // enumerate over the branching sticks
        size_t j = 0;
        double cum_prod = 1.0;
        double begin = 0.0;
        while (true) {
            if (j >= psi_sticks.size()) {
                // draw new psi-stick
                double psi_j = beta(random, 1, params.gamma);
                psi_sticks.push_back(psi_j);
            }
            double interval_length = cum_prod * psi_sticks[j];
            cout << "psi_" << j << ": " << psi_sticks[j] << endl;
            cout << "interval: " << begin << ", " << (begin + interval_length) << endl;
            if (u > begin && u < (begin + interval_length)) {
                // found the corresponding interval
                u = (u - begin) / interval_length;
                break;
            }
            begin += interval_length;
            cum_prod *= (1 - psi_sticks[j]);
            j++;
        }

        // if j corresponds to a new node, draw a nu stick for it, draw frequency param and update data structures accordingly
        // update local variables node_str, u, nu, psi_sticks
        node_str += to_string(j);
        cout << node_str << endl;
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

vector<SomaticMutation> *read_test_ssm()
{
    //char *dir = getcwd(NULL, 0);
    //printf("Current dir: %s", dir);

    // use the sample data from phylowgs
    ifstream file;
    file.open("/Users/seonghwanjun/Google Drive/Research/smc-research/repos/spf/input/test_ssm_data.txt");
    string line;
    vector<string> fields1;
    vector<string> fields2;
    vector<SomaticMutation> *ret = new vector<SomaticMutation>();
    if (file.is_open()) {
        unsigned int line_no = 0;
        while (getline(file, line)) {
            line_no++;
            if (line_no == 1)
                continue;

            //cout << line << endl;
            boost::split(fields1, line, boost::is_any_of("\t"));
            boost::split(fields2, fields1[2], boost::is_any_of(","));
            vector<unsigned int> a;
            for (size_t i = 0; i < fields2.size(); i++) {
                a.push_back(stoi(fields2[i]));
            }
            boost::split(fields2, fields1[3], boost::is_any_of(","));
            vector<unsigned int> d;
            for (size_t i = 0; i < fields2.size(); i++) {
                d.push_back(stoi(fields2[i]));
            }
            SomaticMutation ssm(fields1[0], a, d);
            //cout << ssm.print() << endl;
            ret->push_back(ssm);
        }
        file.close();
    } else {
        cout << "not here!" << endl;
    }
    return ret;
}

void test_tssb(vector<SomaticMutation> *ssms)
{
    // initialize the problem specification
    long seed = 32;
    gsl_rng* random = generate_random_object(seed);
    int num_particles = 10000;
    SMCOptions *options = new SMCOptions();
    options->essThreshold = 1;
    options->resampling = SMCOptions::ResamplingScheme::STRATIFIED;
    options->resample_last_round = true;
    CancerPhyloParameters params;
    params.gamma = 0.7;
    params.lambda = 1;
    params.alpha_0 = 5;
    params.birth_rate = 0.01;
    params.death_rate = 0.01;
    params.weibull_scale = 1.0;
    params.weibull_shape = 5.0;
    params.sequencing_error_prob = 1e-3;

    TSSBProblemSpecification *problem_spec = new TSSBProblemSpecification(ssms, params);
    SMC<PartialCancerPhylogenyState *> smc(problem_spec, options);
    
    smc.run_smc(random, num_particles);
    
    ParticlePopulation<PartialCancerPhylogenyState *> *pop = smc.get_curr_population();
    vector<PartialCancerPhylogenyState *> *states = pop->get_particles();
    vector<double> *weights = 0;
    if (options->resample_last_round) {
        weights = pop->get_normalized_weights();
    } else {
        weights = pop->get_log_weights();
    }

    unordered_map<double, PartialCancerPhylogenyState*> states_map;
    for (size_t i = 0; i < states->size(); i++) {
        PartialCancerPhylogenyState *state = (*states)[i];
        double loglik = state->get_log_likelihood();
        if (states_map.count(loglik) == 0) {
            states_map[loglik] = state;
            cout << state->get_log_likelihood() << endl;
            cout << state->print() << endl;
        }
    }
}

int main()
{
    //test_gsl();
    //test_boost();
    //test_spf();
    //test_smc();
    //test_resampling_schemes();
    //test_copy_unordered_map();
    /*
    CancerPhyloParameters params;
    params.gamma = 1.0;
    params.alpha0 = 1.7;
    gsl_rng *random = generate_random_object(11);
    test_find_node(random, 0.7, params);
     */
    
    std::cout << std::fixed;
    std::cout << std::setprecision(7);
    
    vector<SomaticMutation> *ssms = read_test_ssm();
    test_tssb(ssms);
    
	return 0;
}


