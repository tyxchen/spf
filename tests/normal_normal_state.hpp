//
//  normal_normal_state.hpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-01-18.
//

#ifndef normal_normal_state_hpp
#define normal_normal_state_hpp

class NormalNormalState
{
    double mu;
    double sigma;
public:
    NormalNormalState(double mu, double sigma);
    inline double get_mu() { return mu; }
    inline double get_sigma() { return sigma; }
    void set_mu(double mu);
};

#endif /* normal_normal_state_hpp */
