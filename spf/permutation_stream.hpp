//
//  permutation_stream.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-04-17.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef permutation_stream_h
#define permutation_stream_h

#include <vector>
#include <gsl/gsl_rng.h>

using namespace std;

class PermutationStream
{
    gsl_rng *random = 0;
    long seed;
    unsigned int *indices;
    unsigned int size;
    unsigned int num_calls;
public:
    PermutationStream(unsigned int size);
    unsigned int pop();
    void reset();
    void set_seed(long seed);
    inline gsl_rng *get_random() { return random; }
    ~PermutationStream();
};

PermutationStream::PermutationStream(unsigned int size)
{
    this->size = size;
    this->indices = new unsigned int[size];
    this->num_calls = 0;
    this->random = generate_random_object(seed);
    for (int i = 0; i < size; i++) {
        indices[i] = i;
    }
}

void PermutationStream::reset()
{
    gsl_rng_set(random, seed);
    this->num_calls = 0;
    for (int i = 0; i < size; i++) {
        indices[i] = i;
    }
}

void PermutationStream::set_seed(long seed)
{
    this->seed = seed;
    reset();
}

unsigned int PermutationStream::pop()
{
    int idx = num_calls++ % size;
    if (idx == 0) {
        // shuffle
        gsl_ran_shuffle(random, indices, size, sizeof(unsigned int));
    }
    return indices[idx];
    
}

PermutationStream::~PermutationStream()
{
    if (!random)
        gsl_rng_free(random);
    delete indices;
}

#endif /* permutation_stream_h */
