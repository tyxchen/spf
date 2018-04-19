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
    gsl_rng *random;
    long seed;
    unsigned int *indices;
    unsigned int size;
    unsigned int num_calls;
public:
    PermutationStream(unsigned int size, long seed);
    unsigned int pop();
    void reset();
    inline gsl_rng *get_random() { return random; }
    ~PermutationStream();
};

PermutationStream::PermutationStream(unsigned int size, long seed)
{
    this->seed = seed;
    this->size = size;
    this->indices = new unsigned int[size];
    reset();
}

void PermutationStream::reset()
{
    this->random = generate_random_object(seed);
    this->num_calls = 0;
    for (int i = 0; i < size; i++) {
        indices[i] = i;
    }

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
    delete indices;
}

#endif /* permutation_stream_h */
