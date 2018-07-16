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
#include "smc.hpp"
#include "smc_model.hpp"
#include "../tests/test.hpp"

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

void test_boost()
{
    string s = "a,b, c ,,e,f,";
    vector <string> fields;
    
    cout << "Original = \"" << s << "\"\n\n";
    
    cout << "Split on \',\' only\n";
    boost::split( fields, s, boost::is_any_of( "," ) );
    cout << fields[0] << endl;
}


int main()
{
    //test_gsl();
    //test_boost();
    //test_resampling_schemes();
    //test_copy_unordered_map();
    //test_spf();
    //test_smc();
    test_pmcmc();
	return 0;
}


