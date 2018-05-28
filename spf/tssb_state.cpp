//
//  tssb_state.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-23.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <cmath>
#include <queue>
#include <gsl/gsl_randist.h>
#include "tssb_state.hpp"
#include "sampling_utils.hpp"

PartialCancerPhylogenyState::PartialCancerPhylogenyState(vector<SomaticMutation> data_points)
{
    this->unassigned_data_points = &data_points; // for vector, assignment results in deep copy of the elements
}

PartialCancerPhylogenyState::PartialCancerPhylogenyState(PartialCancerPhylogenyState &src)
{
    // make a deep copy of the state variables
    // vectors and maps will copy over the elements via overloading of = operator
    this->unassigned_data_points = src.unassigned_data_points;
    this->assigned_data_points = src.assigned_data_points;
    this->node2nu_stick = src.node2nu_stick;
    this->node2psi_sticks = src.node2psi_sticks;
    this->node2data = src.node2data;
    this->node2freq = src.node2freq;
    this->instantiated_nodes = src.instantiated_nodes;
}

double PartialCancerPhylogenyState::assign_data_point_helper(gsl_rng *random, double u, int idx, CancerPhyloParameters &params)
{
    // find the node corresponding to u (similar to algorithm of Adams et. al. (2010))
    // generate copy number variation via forward simulation using birth death process
    // compute and return the likelihood
    
    // 1. assigning the data to a node:
    // start with the root, i.e., node_str = "0"
    // enumerate the branching sticks i in node2psi_sticks[node_str] and compute the corresponding interval, check if u falls in the interval
    // extend node_str as node_str += to_string(i);
    // new psi-sticks and nu-sticks are drawn as needed (i.e., new nodes being created)
    string node_str = "0";
    double nu = 0.0;
    while (true) {
        if (!(*node2nu_stick).count(node_str)) {
            // draw new nu stick and sample frequency param for this node
            (*node2nu_stick)[node_str] = beta(random, 1, params.alpha(node_str)); // nu-stick
            string parent_node_str = node_str.substr(0, node_str.length() - 1);
            (*node2freq)[node_str] = uniform(random, 0.0, (*node2freq)[parent_node_str]); // sample frequency from the prior
            (*node2freq)[parent_node_str] -= (*node2freq)[node_str]; // update the parent's frequency
            // TODO: sample the frequencies outside of the outer-while using say nested SMC?
        }
        nu = (*node2nu_stick)[node_str];
        if (u < nu) {
            // found the node!!
            break;
        }
        u = (u - nu) / (1 - nu);

        // enumerate over the branching sticks
        vector<double> psi_sticks = (*node2psi_sticks)[node_str];
        size_t j = 0;
        double cum_prod = 1.0;
        double begin = 0.0;
        while (true) {
            if (j >= psi_sticks.size()) {
                // draw new psi-stick
                // TODO: perform update of data structures at this point?
                double psi_j = beta(random, 1, params.gamma);
                psi_sticks.push_back(psi_j);
            }
            double interval_length = cum_prod * psi_sticks[j];
            if (u > begin && u < (begin + interval_length)) {
                // found the corresponding interval
                u = (u - begin) / interval_length;
                break;
            }
            begin += interval_length;
            cum_prod *= (1 - psi_sticks[j]);
            j++;
        }

        // if j corresponds to a new node, draw a nu stick for it, draw frequency param and update data structures accordingly
        // update local variables node_str, u, nu, psi_sticks
        node_str += to_string(j);
    }

    // assign the data point to the node
    SomaticMutation datum = (*unassigned_data_points)[idx];
    (*assigned_data_points)[datum] = node_str;
    (*node2data)[node_str].push_back(datum);
    (*unassigned_data_points).erase((*unassigned_data_points).begin() + idx);

    // compute the likelihood by simulating the copy number profile and all descendents with nu-stick instantiated
    return compute_log_likelihood(random, datum, node_str, params);
}

double PartialCancerPhylogenyState::compute_log_likelihood(gsl_rng *random, SomaticMutation &datum, string node_str, CancerPhyloParameters &params)
{
    double prevalence = 0.0;
    double prob_obs_ref = 0.0; // probability of observing reference allele
    queue<string> q;
    size_t len;
    q.push(node_str);
    // note:
    // we simulate the number of reference alleles from root to the parent node of node_str
    // initially there are 2 reference alleles
    // we subtract by 1 since one of them turns to variant
    // for variant allele, it starts off with 1
    unsigned int cn_ref = sample_birth_death_process(random, 2, node_str.length()-1, params.birth_rate, params.death_rate) - 1;
    unsigned int cn_var = 1;
    while (q.size() > 0) {
        string curr_node_str = q.front();
        q.pop();

        if (!instantiated_nodes->count(curr_node_str)) { // check to ensure that node is already instantiated, if not, instantiate it
            (*instantiated_nodes)[node_str] = new Node(curr_node_str);
        }

        Node curr_node(*(*instantiated_nodes)[node_str]); // make a copy of the Node object
        cn_ref = sample_birth_death_process(random, cn_ref, 1, params.birth_rate, params.death_rate);
        cn_var = sample_birth_death_process(random, cn_var, 1, params.birth_rate, params.death_rate);
        curr_node.set_cn_profile(datum, cn_ref, cn_var);
        (*instantiated_nodes)[node_str] = &curr_node; // update the instantiated_nodes to point to the newly modified copy

        prevalence += (*node2freq)[curr_node_str];
        prob_obs_ref += (*node2freq)[curr_node_str] * cn_ref / (cn_ref + cn_var);

         // enumerate over the psi-sticks for the current node
        len = this->node2psi_sticks->size();
        for (size_t i = 0; i < len; i++) {
            // if nu-stick is sampled, then there is at least one data point assigned to this node or one of its descendants
            if (node2nu_stick->count(curr_node_str + to_string(i))) {
                q.push(curr_node_str + to_string(i)); // add to the queue
            }
        }
    }
    prob_obs_ref += (1 - prevalence);
    return log2(gsl_ran_binomial_pdf(datum.get_a(), prob_obs_ref, datum.get_d()));

    // TODO: use nested SMC idea to generate the copy number variation?
}

std::pair<PartialCancerPhylogenyState *, double> PartialCancerPhylogenyState::assign_data_point(gsl_rng *random, CancerPhyloParameters &params)
{
    // make a copy of the current state
    PartialCancerPhylogenyState *new_state = new PartialCancerPhylogenyState(*this);

    // 1. sample a data point to be assigned
    // 2. sample U ~ Uniform(0, 1);
    // 3. find the node corresponding to U = u, generate cnv profile along the way and also, generate new nodes and sticks as needed
    int idx = discrete_uniform(random, unassigned_data_points->size());
    double u = uniform(random);
    double logw = assign_data_point_helper(random, u, idx, params);
    return make_pair(new_state, logw);
}
