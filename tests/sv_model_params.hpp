//
//  sv_model_params.h
//  sv
//
//  Created by Seong-Hwan Jun on 2018-11-25.
//

#ifndef sv_model_params_hpp
#define sv_model_params_hpp

class SVModelParams
{
public:
    double phi;
    double beta;
    double sigma;
    inline SVModelParams(double phi, double sigma, double beta) {
        this->phi = phi;
        this->sigma = sigma;
        this->beta = beta;
    }
};

#endif /* sv_model_params_hpp */
