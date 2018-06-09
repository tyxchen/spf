//
//  tssb_state.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-23.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <cmath>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <gsl/gsl_randist.h>

#include "tssb_state.hpp"
#include "numerical_utils.hpp"
#include "sampling_utils.hpp"

PartialCancerPhylogenyState::PartialCancerPhylogenyState(vector<SomaticMutation> *data_points)
{
    num_samples = data_points->at(0).num_samples();
    this->unassigned_data_points = data_points; // for vector, assignment results in deep copy of the elements
    this->str2node = new unordered_map<string, Node*>();
    this->datum2node = new unordered_map<SomaticMutation, string, hash<SomaticMutation>>();
    this->datum2cnr = new unordered_map<SomaticMutation, double, hash<SomaticMutation>>();
    this->datum2cnv = new unordered_map<SomaticMutation, double, hash<SomaticMutation>>();
    this->eta = new vector<double>(num_samples);
    for (size_t i = 0; i < num_samples; i++) {
        (*eta)[i] = 1.0;
    }
}

PartialCancerPhylogenyState::PartialCancerPhylogenyState(const PartialCancerPhylogenyState &src)
{
    // make a deep copy of the state variables
    num_samples = src.num_samples;
    this->unassigned_data_points = new vector<SomaticMutation>(*src.unassigned_data_points);
    this->str2node = new unordered_map<string, Node*>(*src.str2node);
    this->datum2node = new unordered_map<SomaticMutation, string, hash<SomaticMutation>>(*src.datum2node);
    this->datum2cnr = new unordered_map<SomaticMutation, double, hash<SomaticMutation>>(*src.datum2cnr);
    this->datum2cnv = new unordered_map<SomaticMutation, double, hash<SomaticMutation>>(*src.datum2cnv);
    this->eta = new vector<double>(*src.eta);
}

PartialCancerPhylogenyState::~PartialCancerPhylogenyState()
{
    delete unassigned_data_points;
    delete str2node;
    delete datum2node;
    delete eta;
}

double PartialCancerPhylogenyState::sample_next_state(gsl_rng *random, CancerPhyloParameters &params)
{
    // 1. sample a data point to be assigned
    // 2. sample U ~ Uniform(0, 1);
    // 3. find the node corresponding to U = u, generate cnv profile along the way and also, generatne new nodes and sticks as needed
    int idx = discrete_uniform(random, unassigned_data_points->size());
    SomaticMutation &datum = (*unassigned_data_points)[idx];
    double u = uniform(random);
    assign_data_point(random, u, datum, params); // calls sample_frequency as necessary
    this->loglik = compute_log_likelihood(random, datum, params);
    (*unassigned_data_points).erase((*unassigned_data_points).begin() + idx); // delete the idx from the unassigned_data_points
    return loglik;
}

void PartialCancerPhylogenyState::assign_data_point(gsl_rng *random, double u, SomaticMutation &datum, CancerPhyloParameters &params)
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
    vector<Node*> *path = new vector<Node*>();
    while (true) {
        if (!str2node->count(node_str)) { // check if Node object corresponding to node_str exists in the map
            // Node object does not exist yet, instantiate one
            node = new Node(node_str, num_samples);
            // Draw and set nu-stick
            node->set_nu_stick(beta(random, 1, params.alpha(node_str)));
        } else {
            // make a copy of the node since its frequency will be modified in place
            node = new Node(*(*str2node)[node_str]);
        }
        // update str2node[node_str] (even if it already exists) since node object will change
        (*str2node)[node_str] = node;
        path->push_back(node);

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
    (*datum2node)[datum] = node_str;
    // update SSMs assigned to the node
    node->add_ssm(datum);
    if (!node->is_prevalence_sampled()) {
        sample_frequency_helper(random, path, node, datum);
    }
}

void PartialCancerPhylogenyState::sample_frequency_helper(gsl_rng *random, vector<Node*> *path, Node *curr_node, SomaticMutation &datum)
{
    vector<double> phi_parent(num_samples, 1.0);
    double phi = 0.0;
    for (vector<Node*>::iterator it = path->begin(); it != path->end(); ++it)
    {
        Node *node = *it;
        if (!node->is_prevalence_sampled()) {
            for (size_t i = 0; i < num_samples; i++) {
                phi = uniform(random, 0.0, phi_parent[i]);
                node->add_prevalence(phi);
            }
        }
        for (size_t i = 0; i < num_samples; i++) {
            phi_parent[i] = node->get_prevalence(i);
        }
    }
}

double PartialCancerPhylogenyState::compute_log_likelihood(gsl_rng *random, SomaticMutation &datum, CancerPhyloParameters &params)
{
    return compute_log_likelihood_helper(random, datum, params);

    // enumerate over SomaticMutation in datum2node
    /*
    double logw = 0.0;
    for (unordered_map<SomaticMutation, string>::iterator it = datum2node->begin(); it != datum2node->end(); ++it)
    {
        SomaticMutation datum = it->first;
        logw += compute_log_likelihood_helper(random, datum, params);
    }
    return logw;
     */
}

double PartialCancerPhylogenyState::compute_log_likelihood_helper(gsl_rng *random, SomaticMutation &datum, CancerPhyloParameters &params)
{
    double logw = 0.0;

    Node *node = (*str2node)[(*datum2node)[datum]];

    // propose the average ploidy from Weibull distribution
    double cnr = weibull(random, params.weibull_scale, params.weibull_shape);
    double cnv = weibull(random, params.weibull_scale, params.weibull_shape);
    (*datum2cnr)[datum] = cnr;
    (*datum2cnv)[datum] = cnv;

    // compute the likelihood
    double prob_obs_ref = 0.0; // probability of observing the reference allele
    double prob_obs_var = 0.0; // probability of observing the variant allele
    double prevalence = 0.0;
    for (size_t idx = 0; idx < num_samples; idx++) {
        prevalence = node->get_prevalence(idx);
        if (cnr > 0) {
            prob_obs_ref = prevalence * cnr / (cnr + cnv);
        }
        if (cnv > 0) {
            prob_obs_var = prevalence * cnv / (cnr + cnv);
        }
        prob_obs_ref = (1 - prevalence + prob_obs_ref) * (1-params.sequencing_error_prob) + prob_obs_var * params.sequencing_error_prob;
        prob_obs_var = (1 - prevalence + prob_obs_ref) * params.sequencing_error_prob + prob_obs_var * (1 - params.sequencing_error_prob);
        double logp = log_binomial_pdf(datum.get_a(idx), prob_obs_ref, datum.get_d(idx));
        logp += log_binomial_pdf(datum.get_d(idx) - datum.get_a(idx), prob_obs_var, datum.get_d(idx));
        logw += logp;
    }

    return logw;
}

string PartialCancerPhylogenyState::form_node_string(string curr_node_str, size_t branch)
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

vector<Node *> *PartialCancerPhylogenyState::get_children_nodes(Node *curr_node)
{
    vector<Node *> *children = new vector<Node *>();
    vector<double> *psi_sticks = curr_node->get_psi_sticks();
    size_t len = psi_sticks->size();
    string curr_node_str = curr_node->name;
    string next_node_str;
    for (size_t i = 0; i < len; i++) {
        // if node is instantiated already, then there is at least one data point assigned to the next_node_str or one of its descendant
        // therefore, we need to generate copy number for it
        next_node_str = form_node_string(curr_node_str, i);
        if (str2node->count(next_node_str)) {
            children->push_back((*str2node)[next_node_str]);
        }
    }

    return children;
}

string PartialCancerPhylogenyState::print()
{
    string ret = "state: " + to_string(get_log_likelihood()) + "\n";
    // print node specific infor
    for (unordered_map<string,Node*>::iterator it = str2node->begin(); it != str2node->end(); ++it)
    {
        Node *node = it->second;
        if (node != 0 && node->get_ssms()->size() > 0) {
            ret += node->print() + "\n";
        }
    }
    // print SSM specific info
    for (unordered_map<SomaticMutation,string>::iterator it = datum2node->begin(); it != datum2node->end(); ++it)
    {
        SomaticMutation ssm = it->first;
        Node *node = (*str2node)[it->second];
        ret += "[" + ssm.id + "\n";
        ret += "(" + node->name + ", cnr=" + to_string((*datum2cnr)[ssm]) + ", cnv=" + to_string((*datum2cnv)[ssm]) + ")\n";
        ret += "]\n";
    }
    return ret;
}

/*
 double PartialCancerPhylogenyState::compute_log_likelihood_helper(gsl_rng *random, SomaticMutation &datum, CancerPhyloParameters &params)
 {
 // compute the probability of observing the reference allele
 // case 1: if datum is assigned in previous iteration, new descendant nodes may have been introduced
 // case 2: if datum is assigned in this iteration, need to simulate CNV for the current node as well as descendant nodes
 
 // generate cn_profile for all descendant nodes and compute the prevalence and hence the probability of observing reference allele
 
 //cout << datum.id << "=" << (*datum2node)[datum]->name << endl;
 
 Node *curr_node = 0;
 string curr_node_str = (*datum2node)[datum];
 size_t num_samples = datum.num_samples();
 double prevalence[num_samples];
 double prob_obs_ref[num_samples]; // probability of observing the reference allele
 double prob_obs_var[num_samples]; // probability of observing the variant allele
 queue<string> q;
 size_t len;
 q.push(curr_node_str);
 while (q.size() > 0) {
 curr_node_str = q.front();
 
 if (!str2node->count(curr_node_str)) {
 cerr << "Error: there is a bug in the program!" << endl;
 exit(-1);
 }
 curr_node = new Node(*(*str2node)[curr_node_str]);
 unsigned int cn_ref = 0;
 unsigned int cn_var = 0;
 // check if CN profile exists for datum in node
 if (!curr_node->cnprofile_exists(datum)) {
 // generate cn profile
 // note:
 // first check if there is cn profile generated for datum at the parent node
 string parent_node_str = get_parent_string(curr_node_str);
 Node *parent_node = (*str2node)[parent_node_str];
 if (parent_node_str != "" && parent_node->cnprofile_exists(datum)) {
 cn_ref = parent_node->get_cnr(datum);
 cn_var = parent_node->get_cnv(datum);
 } else {
 // we simulate the number of reference alleles from root to the parent node first
 // initially there are 2 reference alleles
 // we subtract by 1 since one of them turns to variant somewhere along the edge between the parent and curr node
 // for variant allele, it starts off with 1
 int num_copies = 2 + sample_birth_death_process(random, 2, curr_node_str.length() - 1, params.birth_rate, params.death_rate);
 cn_ref = num_copies > 0 ? num_copies - 1 : 0;
 cn_var = num_copies > 0 ? 1 : 0;
 }
 
 // simulate the cn from parent_node to curr_node
 cn_ref += sample_birth_death_process(random, cn_ref, 1, params.birth_rate, params.death_rate);
 cn_var += sample_birth_death_process(random, cn_var, 1, params.birth_rate, params.death_rate);
 curr_node->set_cn_profile(datum, cn_ref, cn_var);
 } else {
 cn_ref = curr_node->get_cnr(datum);
 cn_var = curr_node->get_cnv(datum);
 }
 
 // update the data structure
 (*str2node)[curr_node_str] = curr_node;
 
 if (curr_node->is_freq_sampled()) {
 for (size_t i = 0; i < num_samples; i++) {
 double freq = curr_node->get_frequency(i);
 prevalence[i] += freq;
 if (cn_ref + cn_var == 0) {
 prob_obs_ref[i] += 0.0;
 prob_obs_var[i] += 0.0;
 } else {
 prob_obs_ref[i] += freq * cn_ref / (cn_ref + cn_var);
 prob_obs_var[i] += freq * cn_var / (cn_ref + cn_var);
 }
 }
 }
 
 // enumerate over the psi-sticks for the current node
 vector<double> *psi_sticks = curr_node->get_psi_sticks();
 len = psi_sticks->size();
 string next_node_str;
 for (size_t i = 0; i < len; i++) {
 // if node is instantiated already, then there is at least one data point assigned to the next_node_str or one of its descendant
 // therefore, we need to generate copy number for it
 next_node_str = form_node_string(curr_node_str, i);
 if (str2node->count(next_node_str)) {
 q.push(next_node_str); // add to the queue
 }
 }
 
 q.pop();
 }
 double logw = 0.0;
 //Node *node = (*datum2node)[datum];
 for (size_t i = 0; i < num_samples; i++) {
 prob_obs_ref[i] += (1 - prevalence[i]);
 double logp = log_binomial_pdf(datum.get_a(i), prob_obs_ref[i], datum.get_d(i));
 logp += log_binomial_pdf(datum.get_d(i) - datum.get_a(i), prob_obs_var[i], datum.get_d(i));
 //cout << "Sample " << i << ": p=" << prob_obs_ref[i] << ", freq=" << prevalence[i] << ", data=" << datum.get_a(i) << "/" << datum.get_d(i) << ", cnr=" << node->get_cnr(datum) << ", cnv=" << node->get_cnv(datum) << ", logp=" << logp << endl;
 logw += logp;
 }
 return logw;
 
 // TODO: use nested SMC idea to generate the copy number variation?
 }
 */
