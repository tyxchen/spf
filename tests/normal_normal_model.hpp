//
//  normal_normal_model.hpp
//  SPF
//
//  Test SMC samplers on a simple Normal-Normal conjugate model
//  - known variance and unknown mean with Normal prior
//  - data likelihood is Normal
//  - normlization constant is known in this case, so we can test SMC samplers
//
//  Created by Seong-Hwan Jun on 2019-01-18.
//

#ifndef normal_normal_model_hpp
#define normal_normal_model_hpp

#include <vector>

#include <gsl/gsl_randist.h>

#include "smc_model.hpp"

#include "normal_normal_params.hpp"
#include "normal_normal_state.hpp"

using namespace std;

class NormalNormalModel : public ProblemSpecification<NormalNormalState, NormalNormalHyperParams>
{
    size_t num_iter;
    size_t num_mh_iter;
    const vector<double> &data;
    double sigma;
    double get_temperature(double t);
public:
    NormalNormalModel(size_t num_iter, size_t num_mh_iter, const vector<double> &data, double sigma);
    unsigned long num_iterations();
    std::pair<NormalNormalState, double> *propose_initial(gsl_rng *random, NormalNormalHyperParams &param);
    std::pair<NormalNormalState, double> *propose_next(gsl_rng *random, int t, NormalNormalState curr, NormalNormalHyperParams &param);
    double compute_log_lik(double mu, double sigma);
};

#endif /* normal_normal_model_hpp */
