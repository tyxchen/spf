//
//  sv_model_gibbs_proposal.cpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2018-12-04.
//

#include <math.h>
#include <gsl/gsl_randist.h>
#include "numerical_utils.hpp"

#include "sv_model_gibbs_proposal.hpp"

SVModelGibbsProposal::SVModelGibbsProposal(double a, double b, vector<double> &y) :
y(y)
{
    this->a = a;
    this->b = b;
}

SVModelParams *SVModelGibbsProposal::sample_from_prior(gsl_rng *random)
{
    // 1. sample x ~ IG(5, 2)
    // 2. beta^2 = 1/x
    // 3. beta = sqrt(1/x)
    double beta = sqrt(1./gsl_ran_gamma(random, a, 1.0/b));
    // initialization of the initial parameter does not have to be from the prior
    //double beta = gsl_ran_flat(random, 0.0, 2.0);
    SVModelParams *param = new SVModelParams(1.0, 0.16, beta);
    return param;
}

SVModelParams *SVModelGibbsProposal::propose(gsl_rng *random, SVModelParams *curr, vector<pair<double, double>> *genealogy)
{
    // sample from posterior distribution:
    size_t T = genealogy->size();
    double sum = 0.0;
    for (size_t t = 0; t < T; t++) {
        double x_t = genealogy->at(t).first;
        double y_t = y[t];
        sum += pow(y_t, 2.0)*exp(-x_t);
    }
    double posterior_a = a + T/2.;
    double posterior_b = b + sum/2.;
    double beta = sqrt(1./gsl_ran_gamma(random, posterior_a, 1.0/posterior_b));
    SVModelParams *param = new SVModelParams(1.0, 0.16, beta);
    return param;
}

double SVModelGibbsProposal::log_prior(SVModelParams *curr)
{
    if (curr->beta < 0 || curr->beta > 2) {
        return DOUBLE_NEG_INF;
    }
    return log(gsl_ran_flat_pdf(curr->beta, 0, 2));
}
