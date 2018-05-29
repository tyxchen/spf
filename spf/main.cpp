/*
 * main.cpp
 *
 *  Created on: 6 Mar 2018
 *      Author: seonghwanjun
 */

#include <iostream>
#include <iterator>
#include <math.h>
#include <unordered_map>
#include <vector>
#include <gsl/gsl_randist.h>
#include <string>

#include "discrete_hmm_model.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"
#include "smc.hpp"
#include "smc_model.hpp"
#include "test.hpp"
#include "cancer_phylo_param.hpp"

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

void test_tssb()
{
    
}

int main()
{
    //test_gsl();
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
    
    test_tssb();
	return 0;
}


