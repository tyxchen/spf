//
//  normal_normal_params.hpp
//  SPF
//
//  Created by Seong-Hwan Jun on 2019-01-18.
//

#ifndef normal_normal_params_hpp
#define normal_normal_params_hpp

class NormalNormalHyperParams
{
    double mu_0;
    double sigma_0;
public:
    NormalNormalHyperParams(double mu_0, double sigma_0); // accept hyper parameters
    double get_mu0();
    double get_sigma0();
};

#endif /* normal_normal_params_hpp */
