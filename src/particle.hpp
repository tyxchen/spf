//
//  particle.hpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-01-14.
//

#ifndef particle_hpp
#define particle_hpp

template <class S>
class Particle
{
    S *s; // pointer to the state object
    double log_alpha; // log of weight update function
public:
    S *get_state();
    double get_log_alpha();
    
    Particle(S *s, double log_alpha);
    static Particle<S> *make_particle(S *s, double log_alpha);
};

template <class S>
Particle<S>::Particle(S *s, double log_alpha)
{
    this->s = s;
    this->log_alpha = log_alpha;
}

template <class S>
S *Particle<S>::get_state()
{
    return s;
}

template <class S>
double Particle<S>::get_log_alpha()
{
    return log_alpha;
}

template <class S>
static Particle<S> *make_particle(S *s, double log_alpha)
{
    Particle<S> *particle = new Particle<S>(s, log_alpha);
    return particle;
}

#endif /* particle_hpp */
