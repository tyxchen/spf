/*
 * numerical_utils.hpp
 *
 *  Created on: 6 Mar 2018
 *      Author: seonghwanjun
 */

#ifndef SRC_NUMERICAL_UTILS_HPP_
#define SRC_NUMERICAL_UTILS_HPP_

#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

const double DOUBLE_INF = std::numeric_limits<double>::infinity();
const double DOUBLE_NEG_INF = -std::numeric_limits<double>::infinity();

void normalize_destructively(vector<double> &log_weights);
void normalize_destructively(double *log_weights, int size);
double log_add(double x, double y);
double log_add(double *x, int size);
double log_add(vector<double> x);

template <typename T>
void print_vector(const vector<T> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
}

#endif /* SRC_NUMERICAL_UTILS_HPP_ */
