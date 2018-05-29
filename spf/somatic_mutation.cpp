//
//  somatic_mutation.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-25.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include "somatic_mutation.hpp"

SomaticMutation::SomaticMutation(string ssm_id, vector<unsigned int> &a, vector<unsigned int> &d)
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

unsigned int SomaticMutation::get_a(size_t sample_idx)
{
    return a[sample_idx];
}

unsigned int SomaticMutation::get_d(size_t sample_idx)
{
    return d[sample_idx];
}

size_t SomaticMutation::num_samples()
{
    return a.size();
}
