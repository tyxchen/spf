//
//  Node.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-28.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef Node_hpp
#define Node_hpp

#include <string>
#include <unordered_map>
#include "somatic_mutation.hpp"

using namespace std;

class Node
{
    unordered_map<SomaticMutation, unsigned int> *cn_ref = 0;
    unordered_map<SomaticMutation, unsigned int> *cn_var = 0;
    bool operator==(const Node &other) const
    {
        return (name == other.name);
    }
public:
    string name;
    Node(string name);
    Node(const Node &src);
    void set_cn_profile(SomaticMutation &datum, unsigned int cn_ref, unsigned int cn_var);
    ~Node();
};
namespace std
{
    template <>
    class hash<Node> {
    public:
        size_t operator()(const Node& obj) const {
            return hash<string>()(obj.name);
        }
    };
}

#endif /* Node_hpp */
