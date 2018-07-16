//
//  param.h
//  spf
//
//  abstract class (interface) for parameters (all parameters are to inherit this class)
//
//  Created by Seong-Hwan Jun on 2018-07-16.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef param_h
#define param_h

#include <string>

using namespace std;

class Parameters
{
public:
    virtual string to_string() = 0;
};

#endif /* param_h */
