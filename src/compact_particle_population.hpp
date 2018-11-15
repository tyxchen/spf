//
//  compact_particle_population.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-04-17.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef SRC_COMPACT_PARTICLE_POPULATION_H
#define SRC_COMPACT_PARTICLE_POPULATION_H

#include <math.h>

#include "numerical_utils.hpp"

using namespace std;

template <class P>
class CompactParticlePopulation
{
    double log_sum_of_weights;
    double log_sum_of_square_weights;
    unsigned int num_particles;
public:
    CompactParticlePopulation();
    ~CompactParticlePopulation();
    void add_weight(double logw);
    inline unsigned int get_num_particles() { return num_particles; }
    double get_log_sum_weights();
    double get_log_sum_of_square_weights();
    double ess();
    double logZ();
};

template <class P>
CompactParticlePopulation<P>::CompactParticlePopulation()
{
    log_sum_of_weights = DOUBLE_NEG_INF;
    log_sum_of_square_weights = DOUBLE_NEG_INF;
    num_particles = 0;
}

template <class P>
CompactParticlePopulation<P>::~CompactParticlePopulation()
{
    
}

template <class P>
void CompactParticlePopulation<P>::add_weight(double logw)
{
    log_sum_of_weights = log_add(log_sum_of_weights, logw);
    log_sum_of_square_weights = log_add(log_sum_of_square_weights, 2 * logw);
    num_particles++;
}

template <class P>
double CompactParticlePopulation<P>::get_log_sum_weights()
{
    return log_sum_of_weights;
}

template <class P>
double CompactParticlePopulation<P>::get_log_sum_of_square_weights()
{
    return log_sum_of_square_weights;
}

template <class P>
double CompactParticlePopulation<P>::ess()
{
    return exp(2 * log_sum_of_weights - log_sum_of_square_weights);
}

template <class P>
double CompactParticlePopulation<P>::logZ()
{
    return log_sum_of_weights - log(num_particles);
}
#endif /* SRC_COMPACT_PARTICLE_POPULATION_H */
