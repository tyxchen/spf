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
#include <queue>
#include <vector>
#include <gsl/gsl_rng.h>
#include "cancer_phylo_param.hpp"
#include "node.hpp"
#include "somatic_mutation.hpp"
#include "numerical_utils.hpp"

using namespace std;

class PartialCancerPhylogenyState
{
    vector<SomaticMutation> *unassigned_data_points = 0;
    unordered_map<string, Node *> *str2node = 0;
    unordered_map<SomaticMutation, string, hash<SomaticMutation>> *datum2node = 0;
    vector<double> *eta = 0;
    double loglik = 0;
    size_t num_samples = 0;

    string form_node_string(string curr_node_str, size_t branch);
    string get_parent_string(string curr_node_str);
    void assign_data_point(gsl_rng *random, double u, SomaticMutation &datum, CancerPhyloParameters &params);
    void sample_frequency_helper(gsl_rng *random, SomaticMutation &datum, Node *curr_node);
    double compute_log_likelihood(gsl_rng *random, CancerPhyloParameters &params);
    double compute_log_likelihood_helper(gsl_rng *random, SomaticMutation &datum, CancerPhyloParameters &params);
    vector<Node *> *get_children_nodes(Node *curr_node);
public:
    PartialCancerPhylogenyState() = default;
    PartialCancerPhylogenyState(vector<SomaticMutation> *data_points);
    PartialCancerPhylogenyState(const PartialCancerPhylogenyState &src); // make a deep copy
    ~PartialCancerPhylogenyState();
    double sample_next_state(gsl_rng *random, CancerPhyloParameters &params);
    inline double get_log_likelihood() { return loglik; }
    string print();
};

#endif /* TSSBState_h */
