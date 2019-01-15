/*
 * numerical_utils.cpp
 *
 *  Created on: 7 Mar 2018
 *      Author: seonghwanjun
 */

#include <math.h>
#include <gsl/gsl_sf.h>

#include "numerical_utils.hpp"

double log_binomial_pdf(const unsigned int k, const double p, const unsigned int n)
{
    if (k > n)
    {
        return 0;
    }
    else
    {
        double logP;
        
        if (p == 0)
        {
            logP = (k == 0) ? 0 : DOUBLE_NEG_INF;
        }
        else if (p == 1)
        {
            logP = (k == n) ? 0 : DOUBLE_NEG_INF;
        }
        else
        {
            double ln_Cnk = gsl_sf_lnchoose (n, k);
            logP = ln_Cnk + k * log (p) + (n - k) * log1p(-p);
        }
        
        return logP;
    }
}

double normalize(vector<double> &log_weights, vector<double> &weights)
{
    // compute the lognorm
    double log_norm = log_add(log_weights);
    normalize(log_weights, weights, log_norm);
    return log_norm;
}

void normalize(vector<double> &log_weights, vector<double> &weights, double log_norm)
{
    for (unsigned int i = 0; i < log_weights.size(); i++)
    {
        weights[i] = exp(log_weights[i] - log_norm);
    }
}

double normalize_destructively(vector<double> &log_weights)
{
    // compute the lognorm
    double log_norm = log_add(log_weights);
    for (unsigned int i = 0; i < log_weights.size(); i++)
    {
        log_weights[i] = exp(log_weights[i] - log_norm);
    }
    return log_norm;
}

double normalize_destructively(double *log_weights, int size)
{
	// compute the lognorm
	double log_norm = log_add(log_weights, size);
	for (int i = 0; i < size; i++)
	{
		log_weights[i] = exp(log_weights[i] - log_norm);
	}
    return log_norm;
}

double log_add(double x, double y)
{
	// make x the max
	if (y > x) {
		double temp = x;
		x = y;
		y = temp;
	}
	// now x is bigger
	if (x == DOUBLE_NEG_INF) {
		return x;
	}
	double negDiff = y - x;
	if (negDiff < -20) {
		return x;
	}
	return x + log(1.0 + exp(negDiff));
}

double log_add(vector<double> x)
{
    double max = DOUBLE_NEG_INF;
    double maxIndex = 0;
    for (unsigned int i = 0; i < x.size(); i++)
    {
        if (x[i] > max) {
            max = x[i];
            maxIndex = i;
        }
    }
    if (max == DOUBLE_NEG_INF) return DOUBLE_NEG_INF;
    // compute the negative difference
    double threshold = max - 20;
    double sumNegativeDifferences = 0.0;
    for (unsigned int i = 0; i < x.size(); i++) {
        if (i != maxIndex && x[i] > threshold) {
            sumNegativeDifferences += exp(x[i] - max);
        }
    }
    if (sumNegativeDifferences > 0.0) {
        return max + log(1.0 + sumNegativeDifferences);
    } else {
        return max;
    }
    
}

double log_add(double *x, int size)
{
	double max = DOUBLE_NEG_INF;
	double maxIndex = 0;
	for (int i = 0; i < size; i++)
	{
		if (x[i] > max) {
			max = x[i];
			maxIndex = i;
		}
	}
	if (max == DOUBLE_NEG_INF) return DOUBLE_NEG_INF;
	  // compute the negative difference
	  double threshold = max - 20;
	  double sumNegativeDifferences = 0.0;
	  for (int i = 0; i < size; i++) {
	    if (i != maxIndex && x[i] > threshold) {
	      sumNegativeDifferences += exp(x[i] - max);
	    }
	  }
	  if (sumNegativeDifferences > 0.0) {
	    return max + log(1.0 + sumNegativeDifferences);
	  } else {
	    return max;
	  }

}

void add(double *x, double c, size_t size)
{
    for (size_t s = 0; s < size; s++)
    {
        x[s] += c;
    }
}

void multiply(double *x, double c, double *ret, size_t size)
{
    for (size_t s = 0; s < size; s++)
    {
        ret[s] = x[s] * c;
    }
}


