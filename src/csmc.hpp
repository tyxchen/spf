//
//  csmc.hpp
//  spf-lib
//
//  Created by Seong-Hwan Jun on 2018-12-03.
//

#ifndef csmc_h
#define csmc_h

#include <vector>

#include "numerical_utils.hpp"
#include "particle_genealogy.hpp"
#include "resampling.hpp"
#include "sampling_utils.hpp"
#include "smc_model.hpp"
#include "smc_options.hpp"

template <class S, class P>
class ConditionalSMC
{
    SMCOptions *options;
    ProblemSpecification<S, P> *proposal;

    double log_marginal_likelihood = 0;
    double propose(gsl_rng *random, P &params, unsigned int r, ParticleGenealogy<S> *genealogy);
    void resample(gsl_rng *random, unsigned int r, bool fix_last_genealogy, double log_norm);
    ParticleGenealogy<S> *sample_genealogy(const gsl_rng *random);
    double compute_ess(vector<double> &probs);

    // store the particles and ancestor indices in a vector
    vector<vector<S> *> particles;
    vector<vector<double> *> log_weights;
    vector<double> *log_norms;
    vector<vector<int> *> ancestors;

public:
    ConditionalSMC(ProblemSpecification<S, P> *proposal, SMCOptions *options);
    ParticleGenealogy<S> *initialize(P &params);
    ParticleGenealogy<S> *run_csmc(P &params, ParticleGenealogy<S> *genealogy);
    double get_log_marginal_likelihood();
    inline double get_log_norm(unsigned int r) { return (*log_norms)[r]; }
    inline vector<S> *get_particles(unsigned int r) { return particles[r]; }
    inline vector<double> *get_log_weights(unsigned int r){ return log_weights[r]; }
};

template <class S, class P>
ConditionalSMC<S,P>::ConditionalSMC(ProblemSpecification<S, P> *proposal, SMCOptions *options)
{
    this->proposal = proposal;
    this->options = options;
    this->options->init();

    size_t R = proposal->num_iterations();

    // initialize particles, weights, and ancestral indices
    for (size_t r = 0; r < R; r++) {
        particles.push_back(new vector<S>(options->num_particles));
        log_weights.push_back(new vector<double>(options->num_particles));
        log_norms = new vector<double>(R);
        if (r < R - 1) // there are R-1 ancestor indices since no resampling is performed at the last iteration
            ancestors.push_back(new vector<int>(options->num_particles));
    }
}

template <class S, class P>
ParticleGenealogy<S> *ConditionalSMC<S,P>::initialize(P &params)
{
    return run_csmc(params, 0);
}

template <class S, class P>
ParticleGenealogy<S> *ConditionalSMC<S,P>::run_csmc(P &params, ParticleGenealogy<S> *genealogy)
{
    size_t R = proposal->num_iterations();
    log_marginal_likelihood = 0;
    double log_norm = 0.0;
    double log_N = log(options->num_particles);
    for (unsigned int r = 0; r < R; r++)
    {
        log_norm = propose(options->main_random, params, r, genealogy);
        (*log_norms)[r] = log_norm;
        log_marginal_likelihood += (log_norm - log_N);
        
        if (r < R - 1)
            resample(options->resampling_random, r, genealogy != 0, log_norm);
    }

    // sample genealogy
    return sample_genealogy(options->resampling_random);
}

template <class S, class P>
double ConditionalSMC<S,P>::propose(gsl_rng *random, P &params, unsigned int r, ParticleGenealogy<S> *genealogy)
{
    pair<S, double> *ret_val;
    unsigned int N = genealogy == 0 ? options->num_particles : options->num_particles - 1;
    unsigned int parent_idx = 0;
    double log_norm = DOUBLE_NEG_INF;

    vector<S> *samples = new vector<S>();
    vector<double> *logw = new vector<double>();
    for (size_t k = 0; k < N; k++)
    {
        if (r == 0) {
            ret_val = proposal->propose_initial(random, params);
        } else {
            parent_idx = (*ancestors[r-1])[k];
            ret_val = proposal->propose_next(random, r, (*particles[r-1])[parent_idx], params);
        }
        samples->push_back(ret_val->first);
        logw->push_back(ret_val->second);
        log_norm = log_add(log_norm, ret_val->second);
    }
    if (genealogy != 0) {
        // the last particle to be fixed by the given genealogy
        samples->push_back(genealogy->at(r)->first);
        logw->push_back(genealogy->at(r)->second);
        log_norm = log_add(log_norm, genealogy->at(r)->second);
    }
    
    // delete previous particles and log weights
    delete particles[r];
    delete log_weights[r];
    
    particles[r] = samples;
    log_weights[r] = logw;

    return log_norm;
}

// currently supports only multinomial resampling
template <class S, class P>
void ConditionalSMC<S,P>::resample(gsl_rng *random, unsigned int r, bool fix_last_genealogy, double log_norm)
{
    unsigned int N = options->num_particles;

    // normalize the log weights
    vector<double> probs(N);
    normalize(*log_weights[r], probs, log_norm);

    unsigned int num_resamples;
    if (fix_last_genealogy) {
        num_resamples = N-1;
    } else {
        num_resamples = N;
    }
    
    // TODO: instead of using indices, pass in ancestors[r] directly?
    unsigned int *indices = new unsigned int[num_resamples];
    multinomial_resampling(random, &probs, num_resamples, indices);

    for (size_t k = 0; k < N-1; k++)
    {
        (*ancestors[r])[k] = indices[k];
    }
    if (fix_last_genealogy) {
        (*ancestors[r])[N-1] = N-1;
    } else {
        (*ancestors[r])[N-1] = indices[N-1];
    }
    delete [] indices;
}

template <class S, class P>
ParticleGenealogy<S> *ConditionalSMC<S,P>::sample_genealogy(const gsl_rng *random)
{
    size_t R = proposal->num_iterations();
    ParticleGenealogy<S> *genealogy = new ParticleGenealogy<S>(R);

    size_t N = options->num_particles;

    // normalize the log weights
    vector<double> probs(N);
    normalize(*log_weights[R-1], probs, (*log_norms)[R-1]);

    unsigned int idx = multinomial(random, probs);
    for (size_t r = 0; r < R; r++) {
        size_t curr_iter = R - r - 1;
        S particle = (*particles[curr_iter])[idx];
        double log_w = (*log_weights[curr_iter])[idx];
        auto *ret = new pair<S, double>(particle, log_w);
        genealogy->set(curr_iter, ret);
        if (curr_iter > 0) {
            idx = (*ancestors[curr_iter-1])[idx];
        }
    }

    return genealogy;
}

template <class S, class P>
double ConditionalSMC<S,P>::compute_ess(vector<double> &probs)
{
    double ess = 0.0;
    for (size_t k = 0; k < probs.size(); k++) {
        ess += probs[k] * probs[k];
    }
    ess = 1.0/ess;
    ess /= probs.size();
    return ess;
}

template <class S, class P>
double ConditionalSMC<S,P>::get_log_marginal_likelihood()
{
    return log_marginal_likelihood;
}

#endif /* csmc_h */
