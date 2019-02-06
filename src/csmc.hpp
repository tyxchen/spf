//
//  csmc.hpp
//  spf-lib
//
//  Created by Seong-Hwan Jun on 2018-12-03.
//

#ifndef csmc_h
#define csmc_h

#include <memory>
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
    vector<vector<shared_ptr<S> >> particles;
    vector<vector<double>> log_weights;
    vector<vector<int>> ancestors;
    vector<double> log_norms;
    vector<double> probs; // helper variable to avoid having to re-allocate space in memory for resampling

public:
    ConditionalSMC(ProblemSpecification<S, P> *proposal, SMCOptions *options);
    ParticleGenealogy<S> *initialize(P &params);
    ParticleGenealogy<S> *run_csmc(P &params, ParticleGenealogy<S> *genealogy);
    double get_log_marginal_likelihood();
    inline double get_log_norm(unsigned int r) { return log_norms[r]; }
    const S& get_state(unsigned int r, unsigned int k);
    double get_log_weight(unsigned int r, unsigned int k);
    ~ConditionalSMC();
};

template <class S, class P>
ConditionalSMC<S,P>::ConditionalSMC(ProblemSpecification<S, P> *proposal, SMCOptions *options) :
log_norms(proposal->num_iterations()), probs(options->num_particles)
{
    this->proposal = proposal;
    this->options = options;
    this->options->init();

    // initialize and allocate memory to store log weights, ancestral indices, log norms, and helper variable (probs)
    size_t R = proposal->num_iterations();
//    log_norms.reserve(R);
//    probs.reserve(options->num_particles);
    //log_norms = new vector<double>(R);
    //probs = new vector<double>(options->num_particles);

    for (size_t r = 0; r < R; r++) {
        particles.push_back(vector<shared_ptr<S>>(options->num_particles));
        log_weights.push_back(vector<double>(options->num_particles));
        if (r < R - 1) // there are R-1 ancestor indices since no resampling is performed at the last iteration
            ancestors.push_back(vector<int>(options->num_particles));
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
        log_norms[r] = log_norm;
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
    unsigned int N = genealogy == 0 ? options->num_particles : options->num_particles - 1;
    unsigned int parent_idx = 0;
    double log_norm = DOUBLE_NEG_INF;

    // invariant condition: particle_weight_pairs.size() == r or == R
    if (particles.size() != r && particles.size() != proposal->num_iterations()) {
        cerr << "Bug in the implementation of cSMC: " << particles.size() << endl;
        exit(-1);
    }
//    bool particles_initialized = true; // use push_back to insert to states
    // allocate memory to store the particles
//    if (particles.size() == r) {
//        particles.push_back(new vector<shared_ptr<S> >());
//        particles_initialized = false;
//    }

    vector<shared_ptr<S> > &particles_at_r = particles[r];
    vector<double> &log_w = log_weights[r];

    S *s;
    for (size_t k = 0; k < N; k++)
    {
        if (r == 0) {
            s = proposal->propose_initial(random, log_w[k], params);
        } else {
            parent_idx = ancestors[r-1][k];
            S *parent_particle = particles[r-1].at(parent_idx).get();
            s = proposal->propose_next(random, r, *parent_particle, log_w[k], params);
        }

        particles_at_r[k].reset(s);
//        if (!particles_initialized) {
//            particles_at_r.push_back(shared_ptr<S>(s));
//        } else {
//            particles_at_r[k].reset(s);
//        }
        log_norm = log_add(log_norm, log_w[k]);
    }
    if (genealogy != 0) {
        // the last particle is to be fixed by the given genealogy
        particles_at_r[N] = genealogy->get_state_ptr_at(r);
//        if (!particles_initialized) {
//            particles_at_r.push_back(genealogy->get_state_ptr_at(r));
//        } else {
//            particles_at_r[N] = genealogy->get_state_ptr_at(r);
//        }
        log_w[N] = genealogy->get_log_weight_at(r);
        log_norm = log_add(log_norm, log_w[N]);
    }

    return log_norm;
}

// currently supports only multinomial resampling
template <class S, class P>
void ConditionalSMC<S,P>::resample(gsl_rng *random, unsigned int r, bool fix_last_genealogy, double log_norm)
{
    unsigned int N = options->num_particles;

    // normalize the log weights
    normalize(log_weights[r], probs, log_norm);

    unsigned int num_resamples = fix_last_genealogy ? N - 1: N;
    
    // TODO: instead of using indices, pass in ancestors[r] directly?
    unsigned int *indices = new unsigned int[num_resamples];
    multinomial_resampling(random, &probs, num_resamples, indices);

    for (size_t k = 0; k < N-1; k++)
    {
        ancestors[r][k] = indices[k];
    }
    if (fix_last_genealogy) {
        ancestors[r][N-1] = N-1;
    } else {
        ancestors[r][N-1] = indices[N-1];
    }
    delete [] indices;
}

template <class S, class P>
ParticleGenealogy<S> *ConditionalSMC<S,P>::sample_genealogy(const gsl_rng *random)
{
    size_t R = proposal->num_iterations();

    // normalize the log weights
    normalize(log_weights[R-1], probs, log_norms[R-1]);
    unsigned int idx = multinomial(random, probs);
    
    ParticleGenealogy<S> *genealogy = new ParticleGenealogy<S>();
    for (size_t r = 0; r < R; r++) {
        size_t curr_iter = R - r - 1;
        genealogy->set(particles[curr_iter].at(idx), log_weights[curr_iter].at(idx));
        if (curr_iter > 0) {
            idx = ancestors[curr_iter-1][idx];
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

template <class S, class P>
const S& ConditionalSMC<S,P>::get_state(unsigned int r, unsigned int k)
{
    return *particles[r].at(k).get();
}

template <class S, class P>
double ConditionalSMC<S,P>::get_log_weight(unsigned int r, unsigned int k)
{
    return log_weights[r].at(k);
}

template <class S, class P>
ConditionalSMC<S,P>::~ConditionalSMC()
{
}

#endif /* csmc_h */
