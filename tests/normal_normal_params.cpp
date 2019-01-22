//
//  normal_normal_params.cpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-01-18.
//

#include "normal_normal_params.hpp"

NormalNormalHyperParams::NormalNormalHyperParams(double mu_0, double sigma_0) :
mu_0(mu_0), sigma_0(sigma_0)
{
    
}

double NormalNormalHyperParams::get_mu0()
{
    return mu_0;
}

double NormalNormalHyperParams::get_sigma0()
{
    return sigma_0;
}
