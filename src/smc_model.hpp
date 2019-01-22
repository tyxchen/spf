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

#include "particle.hpp"

template <class S, class P> class ProblemSpecification
{
public:
	virtual unsigned long num_iterations() = 0;
	// propose new sample and return it with log weight
	virtual std::pair<S, double> *propose_initial(gsl_rng *random, P &params) = 0;
    // design decision to ponder upon:
    // instead of passing in const S &curr, pass S, which means that the original data
    // is not to be modified
    // this approach is nice if default copy constructor can be used (which might be the case for most instances) but otherwise, the user needs to implement a copy constructor
    // for the object
    // alternative is to pass const reference: const S &curr.
    // in this case, the user has to first make a copy of the state and modify the copied state
    // in the function body
	virtual std::pair<S, double> *propose_next(gsl_rng *random, int t, S curr, P &params) = 0;
    virtual ~ProblemSpecification() { }
};

#endif /* SRC_SMC_MODEL_HPP_ */
