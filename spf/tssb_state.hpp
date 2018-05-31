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
    unordered_map<string, Node *, hash<Node>> *instantiated_nodes = 0;
    unordered_map<SomaticMutation, string, hash<SomaticMutation>> *assigned_data_points = 0;
    unordered_map<string, vector<SomaticMutation>> *node2data = 0;
    unordered_map<string, vector<double>> *node2freq = 0;
    unordered_map<string, double> *node2nu_stick = 0;
    unordered_map<string, vector<double>> *node2psi_sticks = 0;
    double assign_data_point_helper(gsl_rng *random, double u, int idx, CancerPhyloParameters &params);
    void sample_frequncy(gsl_rng *random, size_t num_samples, string curr_node_str, string parent_node_str);
    double compute_log_likelihood(gsl_rng *random, SomaticMutation &datum, string node_str, CancerPhyloParameters &params);
public:
    PartialCancerPhylogenyState() = default;
    PartialCancerPhylogenyState(vector<SomaticMutation> *data_points);
    PartialCancerPhylogenyState(PartialCancerPhylogenyState &src); // make a deep copy
    double assign_data_point(gsl_rng *random, CancerPhyloParameters &params);
    string print();
};

#endif /* TSSBState_h */
