//
//  somatic_mutation.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-23.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef somatic_mutation_h
#define somatic_mutation_h

#include <string>

using namespace std;

class SomaticMutation
{
    unsigned int a;
    unsigned int d;
public:
    string id;
    SomaticMutation(string id, unsigned int a, unsigned int d);
    SomaticMutation(const SomaticMutation &src);
    unsigned int get_a() { return a; }
    unsigned int get_d() { return d; }
    bool operator==(const SomaticMutation &other) const
    {
        return (id == other.id);
    }
};
namespace std
{
    template <>
    class hash<SomaticMutation> {
    public:
        size_t operator()(const SomaticMutation& obj) const {
            return hash<string>()(obj.id);
        }
    };
}

#endif /* somatic_mutation_h */
