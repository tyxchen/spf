//
//  somatic_mutation.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-25.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include "somatic_mutation.hpp"

SomaticMutation::SomaticMutation(string ssm_id, unsigned int a, unsigned int d)
{
    this->a = a;
    this->d = d;
    this->id = ssm_id;
}

SomaticMutation::SomaticMutation(const SomaticMutation &src)
{
    this->a = src.a;
    this->d = src.d;
    this->id = src.id;
}

