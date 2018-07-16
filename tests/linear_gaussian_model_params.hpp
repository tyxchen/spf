//
//  linear_gaussian_model_params.hpp
//  spf
//
//  A, B are vectors of dimension d
//  a, b are scalars
//
//  Created by Seong-Hwan Jun on 2018-07-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef linear_gaussian_model_params_hpp
#define linear_gaussian_model_params_hpp

#include <stdio.h>
#include <vector>

#include "param.hpp"

using namespace std;

class LinearGaussianModelParams : public Parameters
{
public:
    double A;
    double B;
    double a;
    double b;
    double nu;
    double tau;
    double sigma;
    string to_string();
};

#endif /* linear_gaussian_model_params_hpp */
