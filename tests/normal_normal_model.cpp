//
//  normal_normal_model.cpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-01-18.
//

#include <math.h>

#include "normal_normal_model.hpp"

NormalNormalModel::NormalNormalModel(size_t num_iter, size_t num_mh_iter, const vector<double> &data, double sigma) :
num_iter(num_iter), num_mh_iter(num_mh_iter), data(data), sigma(sigma)
{
}

unsigned long NormalNormalModel::num_iterations()
{
    return num_iter;
}

double NormalNormalModel::get_temperature(double t)
{
    return t/num_iter;
}

std::pair<NormalNormalState, double> *NormalNormalModel::propose_initial(gsl_rng *random, NormalNormalHyperParams &params)
{
    // propose from prior mu ~ N(mu_0, sigma_0)
    double mu = gsl_ran_gaussian(random, params.get_sigma0()) + params.get_mu0();
    NormalNormalState state(mu, sigma);

    double log_lik = compute_log_lik(mu, sigma);
    if (num_iter == 1) {
        // essentially, importance sampling
        return new pair<NormalNormalState, double>(state, log_lik);
    } else {
        return propose_next(random, 0, state, params);
    }
}

std::pair<NormalNormalState, double> *NormalNormalModel::propose_next(gsl_rng *random, int t, NormalNormalState curr, NormalNormalHyperParams &params)
{
    double prev_mu = curr.get_mu(); // save the previous value
    for (size_t i = 0; i < num_mh_iter; i++) {
        double mu_star = gsl_ran_gaussian(random, 0.7) + curr.get_mu();
        double log_lik_new = compute_log_lik(mu_star, sigma);
        double log_lik_old = compute_log_lik(curr.get_mu(), sigma);
        double log_prior_new = log(gsl_ran_gaussian_pdf(mu_star - params.get_mu0(), params.get_sigma0()));
        double log_prior_old = log(gsl_ran_gaussian_pdf(curr.get_mu() - params.get_mu0(), params.get_sigma0()));
        double a_ratio = get_temperature(t+1) * (log_lik_new - log_lik_old) + (log_prior_new - log_prior_old);
        double u = log(gsl_rng_uniform(random));
        if (u < a_ratio) {
            curr.set_mu(mu_star);
        }
    }
    double temperature_diff = get_temperature(t+1) - get_temperature(t);
    double log_alpha = temperature_diff * compute_log_lik(prev_mu, sigma);
    return new pair<NormalNormalState, double>(curr, log_alpha);
}

double NormalNormalModel::compute_log_lik(double mu, double sigma)
{
    double log_lik = 0.0;
    for (size_t i = 0; i < data.size(); i++) {
        log_lik += log(gsl_ran_gaussian_pdf(data[i] - mu, sigma));
    }
    return log_lik;
}

