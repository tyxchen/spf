//
//  sv_model.hpp
//  sv
//
//  x_0 ~ N(0, sigma^2/(1-phi^2))
//  x_t | x_{t-1} ~ N(\phi x_t, sigma^2)
//  y_t | x_t ~ N(0, \beta^2 exp(x_t))
//
//  Created by Seong-Hwan Jun on 2018-11-25.
//

#ifndef sv_model_hpp
#define sv_model_hpp

#include <stdio.h>
#include <vector>

#include "particle.hpp"
#include "smc_model.hpp"
#include "sv_model_params.hpp"

using namespace std;

class SVModel : public ProblemSpecification<double, SVModelParams>
{
    vector<double> obs;
public:
    SVModel(vector<double> &obs);
    unsigned long num_iterations();
    double *propose_initial(gsl_rng *random, double &log_w, SVModelParams &params);
    double *propose_next(gsl_rng *random, int t, const double &curr, double &log_w, SVModelParams &params);
    static void generate_data(gsl_rng *random, size_t T, SVModelParams &params, vector<double> &latent, vector<double> &obs);
    ~SVModel();
};

#endif /* sv_model_hpp */
