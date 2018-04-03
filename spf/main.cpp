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
    test_smc();
	return 0;
}



