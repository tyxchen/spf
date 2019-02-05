//
//  particle_genealogy.hpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-01-24.
//

#ifndef particle_genealogy_hpp
#define particle_genealogy_hpp

#include <vector>

using namespace std;

template <class S>
class ParticleGenealogy
{
    vector<shared_ptr<S>> state_genealogy;
    vector<double> log_weights_genealogy;
public:
    ParticleGenealogy();
    size_t size();
    void set(shared_ptr<S> state, double log_w);
    const S &get_state_at(size_t r);
    double get_log_weight_at(size_t r);
    shared_ptr<S> get_state_ptr_at(size_t r);
};

template <class S>
ParticleGenealogy<S>::ParticleGenealogy()
{
    
}

template <class S>
size_t ParticleGenealogy<S>::size()
{
    return state_genealogy.size();
}

template <class S>
void ParticleGenealogy<S>::set(shared_ptr<S> state, double log_w)
{
    state_genealogy.push_back(state);
    log_weights_genealogy.push_back(log_w);
}

template <class S>
const S &ParticleGenealogy<S>::get_state_at(size_t r)
{
    return *state_genealogy[r].get();
}

template <class S>
double ParticleGenealogy<S>::get_log_weight_at(size_t r)
{
    return log_weights_genealogy[r];
}

template <class S>
shared_ptr<S> ParticleGenealogy<S>::get_state_ptr_at(size_t r)
{
    return state_genealogy[r];
}

#endif /* particle_genealogy_hpp */
