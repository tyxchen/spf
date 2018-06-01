//
//  tssb_state.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-23.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <cmath>
#include <queue>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <gsl/gsl_randist.h>

#include "tssb_state.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"

PartialCancerPhylogenyState::PartialCancerPhylogenyState(vector<SomaticMutation> *data_points)
{
    this->unassigned_data_points = data_points; // for vector, assignment results in deep copy of the elements
    this->str2node = new unordered_map<string, Node*>();
    this->datum2node = new unordered_map<SomaticMutation, Node *, hash<SomaticMutation>>();
}

PartialCancerPhylogenyState::PartialCancerPhylogenyState(const PartialCancerPhylogenyState &src)
{
    // make a deep copy of the state variables
    this->unassigned_data_points = new vector<SomaticMutation>(*src.unassigned_data_points);
    this->str2node = new unordered_map<string, Node*>(*src.str2node);
    this->datum2node = new unordered_map<SomaticMutation, Node *, hash<SomaticMutation>>(*src.datum2node);
}

PartialCancerPhylogenyState::~PartialCancerPhylogenyState()
{
    delete unassigned_data_points;
    delete str2node;
    delete datum2node;
}

double PartialCancerPhylogenyState::assign_data_point(gsl_rng *random, CancerPhyloParameters &params)
{
    // 1. sample a data point to be assigned
    // 2. sample U ~ Uniform(0, 1);
    // 3. find the node corresponding to U = u, generate cnv profile along the way and also, generatne new nodes and sticks as needed
    int idx = discrete_uniform(random, unassigned_data_points->size());
    double u = uniform(random);
    double logw = assign_data_point_helper(random, u, idx, params);
    (*unassigned_data_points).erase((*unassigned_data_points).begin() + idx); // delete the idx from the unassigned_data_points
    return logw;
}

double PartialCancerPhylogenyState::assign_data_point_helper(gsl_rng *random, double u, int idx, CancerPhyloParameters &params)
{
    // find the node corresponding to u (similar to algorithm of Adams et. al. (2010))
    // generate copy number variation via forward simulation using birth death process
    // compute and return the likelihood

    // assigning the data to a node:
    // start with the root, i.e., node_str = "0"
    // enumerate the branching sticks i in node2psi_sticks[node_str] and compute the corresponding interval, check if u falls in the interval
    // extend node_str as node_str += to_string(i);
    // new psi-sticks and nu-sticks are drawn as needed (i.e., new nodes being created)
    string node_str = "0";
    Node *node = 0;
    double nu = 0.0;
    while (true) {
        if (!str2node->count(node_str)) { // check if Node object corr to node_str exists
            // Node object does not exist yet, instantiate one
            // Draw and set nu-stick
            // Sample cluster frequency
            // Update str2node
            node = new Node(node_str);

            node->set_nu_stick(beta(random, 1, params.alpha(node_str)));

            string parent_node_str = get_parent_string(node_str);
            Node *parent_node = 0;
            if (parent_node_str != "")
                parent_node = (*str2node)[parent_node_str];
            sample_frequency(random, (*unassigned_data_points)[idx].num_samples(), node, parent_node);
        } else {
            node = new Node(*(*str2node)[node_str]); // make a copy of the node since its frequency will be modified in place
        }
        (*str2node)[node_str] = node;

        nu = node->get_nu_stick();
        if (u < nu) {
            // found the node!!
            break;
        }
        
        // shrink u relative to the current segment
        u = (u - nu) / (1 - nu);

        // find sub branch by enumerate over the branching sticks
        unsigned int j = node->find_branch(random, u, params);

        // update local variable node_str
        node_str = form_node_string(node_str, j);
    }

    // assign the data point to the node by updating datum2node
    SomaticMutation datum = (*unassigned_data_points)[idx];
    (*datum2node)[datum] = node;

    // compute the likelihood by simulating the copy number profile
    return compute_log_likelihood(random, datum, node_str, params);
}

void PartialCancerPhylogenyState::sample_frequency(gsl_rng *random, size_t num_samples, Node *curr_node, Node *parent_node)
{
    if (curr_node->name == "0") { // sample the initial set of frequencies for the root node, this call takes place for the first data point
        for (size_t i = 0; i < num_samples; i++) {
            curr_node->add_frequency(uniform(random, 0.0, 1.0));
        }
    } else {
        // sample from prior: Uniform(0, parent_freq)
        // no need to consider sibling frequencies because they have already been subtracted out
        for (size_t i = 0; i < num_samples; i++) {
            double parent_freq = parent_node->get_frequency(i);
            double freq = uniform(random, 0.0, parent_freq); // sample frequency from the prior
            curr_node->add_frequency(freq);
            parent_node->set_frequency(i, parent_freq - freq);
        }
    }
}

double PartialCancerPhylogenyState::compute_log_likelihood(gsl_rng *random, SomaticMutation &datum, string node_str, CancerPhyloParameters &params)
{
    size_t num_samples = datum.num_samples();
    double prevalence[num_samples];
    double prob_obs_ref[num_samples]; // probability of observing reference allele
    queue<string> q;
    size_t len;
    q.push(node_str);
    // note:
    // we simulate the number of reference alleles from root to the parent node of node_str
    // initially there are 2 reference alleles
    // we subtract by 1 since one of them turns to variant
    // for variant allele, it starts off with 1
    int num_copies = 2 + sample_birth_death_process(random, 2, node_str.length()-1, params.birth_rate, params.death_rate);
    unsigned int cn_ref = num_copies > 0 ? num_copies - 1 : 0;
    unsigned int cn_var = num_copies > 0 ? 1 : 0;
    while (q.size() > 0) {
        string curr_node_str = q.front();
        q.pop();

        if (!str2node->count(curr_node_str)) {
            cerr << "Error: there is a bug in the program!" << endl;
            exit(-1);
        }
        Node *curr_node = (*str2node)[curr_node_str];

        cn_ref += sample_birth_death_process(random, cn_ref, 1, params.birth_rate, params.death_rate);
        cn_var += sample_birth_death_process(random, cn_var, 1, params.birth_rate, params.death_rate);
        curr_node->set_cn_profile(datum, cn_ref, cn_var);

        for (size_t i = 0; i < num_samples; i++) {
            prevalence[i] += curr_node->get_frequency(i);
            prob_obs_ref[i] += curr_node->get_frequency(i) * cn_ref / (cn_ref + cn_var);
        }

         // enumerate over the psi-sticks for the current node
        vector<double> *psi_sticks = curr_node->get_psi_sticks();
        len = psi_sticks->size();
        string next_node_str;
        for (size_t i = 0; i < len; i++) {
            // if node is instantiated already, then there is at least one data point assigned to the next_node_str or one of its descendant
            // therefore, we need to generate copy number for it
            next_node_str = curr_node_str + to_string(i);
            if (str2node->count(next_node_str)) {
                q.push(next_node_str); // add to the queue
            }
        }
    }
    double logw = 0.0;
    for (size_t i = 0; i < num_samples; i++) {
        prob_obs_ref[i] += (1 - prevalence[i]);
        double logp = log_binomial_pdf(datum.get_a(i), prob_obs_ref[i], datum.get_d(i));
        logw += logp;
    }
    return logw;

    // TODO: use nested SMC idea to generate the copy number variation?
}

string PartialCancerPhylogenyState::form_node_string(string curr_node_str, int branch)
{
    return curr_node_str + "_" + to_string(branch);
}

string PartialCancerPhylogenyState::get_parent_string(string curr_node_str)
{
    // split the string using "_"
    vector<string> fields;
    boost::split(fields, curr_node_str, boost::is_any_of("_"));
    string ret = "";
    size_t num_iter = fields.size() - 1;
    for (size_t i = 0; i < num_iter; i++)
    {
        ret += fields[i];
        if (i < num_iter - 1) {
            ret += "_";
        }
    }
    return ret;
}

string PartialCancerPhylogenyState::print()
{
    string ret = "state: \n";
    for (unordered_map<string,Node*>::iterator it = str2node->begin(); it != str2node->end(); ++it)
    {
        Node *node = it->second;
        ret += node->print() + "\n";
    }
//    for (unordered_map<SomaticMutation,Node*>::iterator it = datum2node->begin(); it != datum2node->end(); ++it)
//    {
//        SomaticMutation ssm = it->first;
//        Node *node = it->second;
//        ret += ssm.id + ": " + node->name + "\n";
//    }
    return ret;
}
