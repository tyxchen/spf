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
    size_t len;
    vector<pair<S, double> *> *genealogy = 0;
public:
    ParticleGenealogy(size_t len);
    pair<S, double> *at(size_t r);
    void set(size_t r, pair<S, double> *ret);
    inline size_t size() { return len; }
};

template <class S>
ParticleGenealogy<S>::ParticleGenealogy(size_t len) :
len(len)
{
    genealogy = new vector<pair<S, double> *>(len);
}

template <class S>
pair<S, double> *ParticleGenealogy<S>::at(size_t r)
{
    return genealogy->at(r);
}

template <class S>
void ParticleGenealogy<S>::set(size_t r, pair<S, double> *ret)
{
    (*genealogy)[r] = ret;
}

#endif /* particle_genealogy_hpp */
