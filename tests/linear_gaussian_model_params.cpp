//
//  linear_gaussian_model_params.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-07-15.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include "linear_gaussian_model_params.hpp"

string LinearGaussianModelParams::to_string()
{
    string str = "A: " + std::to_string(A) + "\n";
    str += "a: " + std::to_string(a) + "\n";
    str += "B: " + std::to_string(B) + "\n";
    str += "b: " + std::to_string(b) + "\n";
    str += "nu: " + std::to_string(nu) + "\n";
    str += "tau: " + std::to_string(tau) + "\n";
    str += "sigma: " + std::to_string(sigma) + "\n";
    return str;
}
