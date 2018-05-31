//
//  Node.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-28.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include "node.hpp"

Node::Node(string name)
{
    this->name = name;
    this->cn_ref = new unordered_map<SomaticMutation, unsigned int>();
    this->cn_var = new unordered_map<SomaticMutation, unsigned int>();
}

Node::Node(const Node &src)
{
    this->name = src.name;
    this->cn_ref = src.cn_ref; // note: assignment operator for unordered_map is overriden to copy the elements over
    this->cn_var = src.cn_var;
}

Node::~Node()
{
    delete cn_ref;
    delete cn_var;
}

void Node::set_cn_profile(SomaticMutation &datum, unsigned int cnr, unsigned int cnv)
{
    (*cn_ref)[datum] = cnr;
    (*cn_var)[datum] = cnv;
}
