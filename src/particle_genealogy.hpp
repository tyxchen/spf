//
//  particle_genealogy.hpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-01-24.
//

#ifndef particle_genealogy_hpp
#define particle_genealogy_hpp

#include <iostream>
#include <memory>
#include <vector>

using namespace std;

template <class S>
class ParticleGenealogy
{
    vector<shared_ptr<S> > state_genealogy;
    vector<double> log_weights_genealogy;
public:
    ParticleGenealogy(size_t length);
    size_t size() const;
    void set(size_t r, shared_ptr<S> state, double log_w);
    const S &get_state_at(size_t r) const;
    double get_log_weight_at(size_t r) const;
    shared_ptr<S> get_state_ptr_at(size_t r);
};

template <class S>
ParticleGenealogy<S>::ParticleGenealogy(size_t length) :
state_genealogy(length), log_weights_genealogy(length)
{
    
}

template <class S>
size_t ParticleGenealogy<S>::size() const
{
    return state_genealogy.size();
}

template <class S>
void ParticleGenealogy<S>::set(size_t r, shared_ptr<S> state, double log_w)
{
    if (r > state_genealogy.size()) {
        cerr << "Error: particle genealogy construction out of index." << endl;
        exit(-1);
    }
    state_genealogy[r] = state;
    log_weights_genealogy[r] = log_w;
}

template <class S>
const S &ParticleGenealogy<S>::get_state_at(size_t r) const
{
    return *state_genealogy[r].get();
}

template <class S>
double ParticleGenealogy<S>::get_log_weight_at(size_t r) const
{
    return log_weights_genealogy[r];
}

template <class S>
shared_ptr<S> ParticleGenealogy<S>::get_state_ptr_at(size_t r)
{
    return state_genealogy[r];
}

#endif /* particle_genealogy_hpp */
