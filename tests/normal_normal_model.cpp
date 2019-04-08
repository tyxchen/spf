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

shared_ptr<NormalNormalState> NormalNormalModel::propose_initial(gsl_rng *random, double &log_w, NormalNormalHyperParams &params)
{
    // propose from prior mu ~ N(mu_0, sigma_0)
    double mu = gsl_ran_gaussian(random, params.get_sigma0()) + params.get_mu0();

    if (num_iter == 1) {
        // essentially, importance sampling
        log_w = compute_log_lik(data, mu, sigma);
        shared_ptr<NormalNormalState> ret(new NormalNormalState(mu, sigma));
        return ret;
    } else {
        NormalNormalState state(mu, sigma);
        return propose_next(random, 0, state, log_w, params);
    }
}

shared_ptr<NormalNormalState> NormalNormalModel::propose_next(gsl_rng *random, unsigned int t, const NormalNormalState &curr, double &log_w, NormalNormalHyperParams &params)
{
    NormalNormalState *new_state = new NormalNormalState(curr.get_mu(), curr.get_sigma());
    log_w = log_weight(t, curr, params);

    for (size_t i = 0; i < num_mh_iter; i++) {
        double mu_star = gsl_ran_gaussian(random, 0.7) + new_state->get_mu();
        double log_lik_new = compute_log_lik(data, mu_star, sigma);
        double log_lik_old = compute_log_lik(data, new_state->get_mu(), sigma);
        double log_prior_new = log(gsl_ran_gaussian_pdf(mu_star - params.get_mu0(), params.get_sigma0()));
        double log_prior_old = log(gsl_ran_gaussian_pdf(new_state->get_mu() - params.get_mu0(), params.get_sigma0()));
        double a_ratio = get_temperature(t+1, num_iter) * (log_lik_new - log_lik_old) + (log_prior_new - log_prior_old);
        double u = log(gsl_rng_uniform(random));
        if (u < a_ratio) {
            new_state->set_mu(mu_star);
        }
    }
    shared_ptr<NormalNormalState> ret(new_state);
    return ret;
}

double NormalNormalModel::log_weight(unsigned int t, const NormalNormalState &state, const NormalNormalHyperParams &params)
{
    double temperature_diff = get_temperature(t+1, num_iter) - get_temperature(t, num_iter);
    return temperature_diff * compute_log_lik(data, state.get_mu(), sigma);
}

double NormalNormalModel::get_temperature(double t, size_t num_iter)
{
    return t/num_iter;
}

double NormalNormalModel::compute_log_lik(const vector<double> &data, double mu, double sigma)
{
    double log_lik = 0.0;
    for (size_t i = 0; i < data.size(); i++) {
        log_lik += log(gsl_ran_gaussian_pdf(data[i] - mu, sigma));
    }
    return log_lik;
}

