//
//  Node.cpp
//  spf
//
//  Created by Seong-Hwan Jun on 2018-05-28.
//  Copyright © 2018 Seong-Hwan Jun. All rights reserved.
//

#include <iostream>

#include "node.hpp"
#include "sampling_utils.hpp"

Node::Node(string name, size_t num_samples)
{
    this->name = name;
    this->ssms = new vector<SomaticMutation>();
    //this->cluster_freq = new vector<double>();
    this->psi_sticks = new vector<double>();
    this->prevalence = new vector<double>();
}

Node::Node(const Node &src)
{
    this->nu = src.nu;
    this->name = src.name;
    this->ssms = new vector<SomaticMutation>(*src.ssms);
    //this->cluster_freq = new vector<double>(*src.cluster_freq);
    this->prevalence = new vector<double>(*src.prevalence);
    this->psi_sticks = new vector<double>(*src.psi_sticks);
}

Node::~Node()
{
    delete ssms;
    //delete cluster_freq;
    delete prevalence;
    delete psi_sticks;
}

void Node::set_nu_stick(double nu)
{
    this->nu = nu;
}

double Node::get_nu_stick()
{
    return this->nu;
}

/*
void Node::add_frequency(double freq)
{
    cluster_freq->push_back(freq);
}

void Node::set_frequency(size_t idx, double freq)
{
    (*cluster_freq)[idx] = freq;
    (*prevalence)[idx] += freq;
}

double Node::get_frequency(size_t idx)
{
    return cluster_freq->at(idx);
}

bool Node::is_freq_sampled()
{
    return cluster_freq->size() > 0;
}
*/

unsigned int Node::find_branch(gsl_rng *random, double u, CancerPhyloParameters &params)
{
    unsigned int j = 0;
    double cum_prod = 1.0;
    double begin = 0.0;
    while (true) {
        if (j >= psi_sticks->size()) {
            // draw new psi-stick
            double psi_j = beta(random, 1, params.gamma);
            psi_sticks->push_back(psi_j);
        }
        double interval_length = cum_prod * (*psi_sticks)[j];
        if (u > begin && u < (begin + interval_length)) {
            // found the corresponding interval
            u = (u - begin) / interval_length;
            break;
        }
        begin += interval_length;
        cum_prod *= (1 - (*psi_sticks)[j]);
        j++;
    }
    return j;
}

vector<double> *Node::get_psi_sticks()
{
    return psi_sticks;
}

void Node::add_ssm(SomaticMutation &datum)
{
    ssms->push_back(datum);
}

vector<SomaticMutation> *Node::get_ssms()
{
    return ssms;
}

double Node::get_prevalence(size_t idx)
{
    return (*prevalence)[idx];
}

void Node::add_prevalence(double u)
{
    prevalence->push_back(u);
}

void Node::update_prevalence(size_t idx, double u)
{
    (*prevalence)[idx] += u;
}

bool Node::is_prevalence_sampled()
{
    return prevalence->size() > 0;
}

string Node::print()
{
    // print the node name, SSMs, and copy number variation
    string ret = this->name;
    ret += "( ";
    for (auto ssm : *ssms)
    {
        ret += ssm.id + " ";
    }
    SomaticMutation ssm = this->ssms->at(0);
    size_t len = ssm.num_samples();
    for (size_t i = 0; i < len; i++)
        ret += ", " + to_string(this->get_prevalence(i));
    ret += ")";
    return ret;
}
