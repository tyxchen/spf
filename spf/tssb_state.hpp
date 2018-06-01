//
//  TSSBState.hpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-23.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#ifndef TSSBState_h
#define TSSBState_h

#include <unordered_map>
#include <vector>
#include <gsl/gsl_rng.h>
#include "cancer_phylo_param.hpp"
#include "node.hpp"
#include "somatic_mutation.hpp"

using namespace std;

class PartialCancerPhylogenyState
{
    vector<SomaticMutation> *unassigned_data_points = 0;
    unordered_map<string, Node *> *str2node = 0;
    unordered_map<SomaticMutation, Node *, hash<SomaticMutation>> *datum2node = 0;

    string form_node_string(string curr_node_str, int branch);
    string get_parent_string(string curr_node_str);
    double assign_data_point_helper(gsl_rng *random, double u, int idx, CancerPhyloParameters &params);
    void sample_frequency(gsl_rng *random, size_t num_samples, Node *curr_node, Node *parent_node);
    double compute_log_likelihood(gsl_rng *random, SomaticMutation &datum, string node_str, CancerPhyloParameters &params);
public:
    PartialCancerPhylogenyState() = default;
    PartialCancerPhylogenyState(vector<SomaticMutation> *data_points);
    PartialCancerPhylogenyState(const PartialCancerPhylogenyState &src); // make a deep copy
    ~PartialCancerPhylogenyState();
    double assign_data_point(gsl_rng *random, CancerPhyloParameters &params);
    string print();
};

#endif /* TSSBState_h */
