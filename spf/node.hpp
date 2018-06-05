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
#include <vector>

#include <gsl/gsl_randist.h>

#include "cancer_phylo_param.hpp"
#include "somatic_mutation.hpp"

using namespace std;

class Node
{
    vector<SomaticMutation> *ssms = 0;
    vector<double> *cluster_freq = 0;
    vector<double> *prevalence = 0;
    vector<double> *psi_sticks = 0;
    double nu = 0.0;
    bool operator==(const Node &other) const
    {
        return (name == other.name);
    }
public:
    string name;
    Node(string name, size_t num_samples);
    Node(const Node &src);
    ~Node();

    void set_nu_stick(double nu);
    double get_nu_stick();
    
    void add_frequency(double freq);
    void set_frequency(size_t idx, double freq);
    double get_frequency(size_t idx);
    bool is_freq_sampled();
    
    double get_prevalence(size_t idx);
    void add_prevalence(double u);
    void update_prevalence(size_t idx, double u);
    
    vector<double> *get_psi_sticks();
    
    unsigned int find_branch(gsl_rng *random, double u, CancerPhyloParameters &params);
    
    /*
    void set_cn_profile(SomaticMutation &datum, unsigned int cn_ref, unsigned int cn_var);
    unsigned int get_cnr(SomaticMutation &datum);
    unsigned int get_cnv(SomaticMutation &datum);
    bool cnprofile_exists(SomaticMutation &datum);
     */
    
    void add_ssm(SomaticMutation &datum);
    vector<SomaticMutation> *get_ssms();
    
    string print();
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
