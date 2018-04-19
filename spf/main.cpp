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
#include "test.hpp"

using namespace std;

int main()
{
    test_spf();
    //test_smc();
    //test_resampling_schemes();
    /*
    gsl_rng* random = get_random(123);
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
    */
	return 0;
}



