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
#include <vector>

using namespace std;

class SomaticMutation
{
    vector<unsigned int> a;
    vector<unsigned int> d;
public:
    string id;
    SomaticMutation(string id, vector<unsigned int> &a, vector<unsigned int> &d);
    SomaticMutation(const SomaticMutation &src);
    unsigned int get_a(size_t sample_idx);
    unsigned int get_d(size_t sample_idx);
    size_t num_samples();
    string print() {
        string ret = id + " ";
        /*
        for (int i = 0; i < a.size(); i++) {
            ret += to_string(a[i]) + "/" + to_string(d[i]) + " ";
        }
         */
        return ret;
    }
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
