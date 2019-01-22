//
//  normal_normal_state.cpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-01-18.
//

#include "normal_normal_state.hpp"

NormalNormalState::NormalNormalState(double mu, double sigma) :
mu(mu), sigma(sigma)
{
    
}

void NormalNormalState::set_mu(double mu)
{
    this->mu = mu;
}
