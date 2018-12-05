//
//  sv_model_gibbs_proposal.hpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2018-12-04.
//

#ifndef sv_model_gibbs_proposal_hpp
#define sv_model_gibbs_proposal_hpp

#include <vector>

#include "pg_proposal.hpp"
#include "sv_model_params.hpp"

using namespace std;

class SVModelGibbsProposal : public PGProposal<double, SVModelParams>
{
    vector<double> &y;
    double a, b;
public:
    SVModelParams *sample_from_prior(gsl_rng *random);
    SVModelParams *propose(gsl_rng *random, SVModelParams *curr, vector<pair<double, double>> *genealogy);
    double log_prior(SVModelParams *curr); // log p(curr)
    
    SVModelGibbsProposal(double a, double b, vector<double> &y);
};

#endif /* sv_model_gibbs_proposal_hpp */
