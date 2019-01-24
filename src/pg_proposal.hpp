//
//  pg_proposal.hpp
//  spf-lib
//
//  Created by Seong-Hwan Jun on 2018-12-04.
//

#ifndef pg_proposal_h
#define pg_proposal_h

#include <gsl/gsl_rng.h>

#include "particle_genealogy.hpp"

template <class S, class P>
class PGProposal
{
public:
    virtual P *sample_from_prior(gsl_rng *random) = 0;
    virtual P *propose(gsl_rng *random, P *curr, ParticleGenealogy<S> *genealogy) = 0;
    virtual double log_prior(P *curr) = 0; // log p(curr)
};

#endif /* pg_proposal_h */
