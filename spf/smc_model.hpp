/*
 * proposal.hpp
 *
 *  Created on: 6 Mar 2018
 *      Author: seonghwanjun
 */

#ifndef SRC_SMC_MODEL_HPP_
#define SRC_SMC_MODEL_HPP_

#include <gsl/gsl_rng.h>

#include <utility>

template <class T> class ProblemSpecification
{
public:
	virtual unsigned long num_iterations() = 0;
	// propose new sample and return it with log weight
	virtual std::pair<T, double> propose_initial(gsl_rng *random) = 0;
	virtual std::pair<T, double> propose_next(gsl_rng *random, int t, T curr) = 0;
	virtual ~ProblemSpecification() { }
};


#endif /* SRC_SMC_MODEL_HPP_ */
