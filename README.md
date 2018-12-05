# Library for Sequential Monte Carlo Algorithms

[![Build Status](https://travis-ci.org/junseonghwan/spf.svg?branch=master)](https://travis-ci.org/junseonghwan/spf)

# Currently supports:
    - SMC Sampler with
        + Adaptive resampling
        + Stratified resampling
        + Multinomial resampling
        + Systematic resampling
    - Conditional SMC with multinomial resampling
    - Particle MCMC:
        - PMMH
        - PG

# Installation (on OSX)
+ Dependencies: `cmake, gsl, boost`
    - `brew install cmake`
    - `brew install gsl`
    - `brew install boost`
+ Then,
```bash
git clone https://github.com/junseonghwan/spf.git
cd spf
mkdir build
cd build
cmake ..
make install
```
The last line will install `libSPF.a` to /usr/local/bin` and header files to `/usr/local/include/spf`.

# Usage
Example usage can be found under `tests/` directory with `DiscreteHMM` and `SVModel`. Additionally, `tests/test.cpp` contains example code.

## Defining state space and parameters
State space refers to the type of the values taken by the latent states X(t). For DiscreteHMM, the state space is `int` and for SVModel, it is `double`. But it could be any user defined type. If there is a single continuous parameter, it can simply be a `double`. In most settings, there are multiple parameters; therefore, it is recommended for the user to create a class representing the parameters. 

## Proposal for SMC
First extend `ProblemSpecification`. Implement the functions and constructor as necessary.
+ `propose_initial`: returns initial state and its log weight.
+ `propose_next`: returns next state and its log weight, given the current state.
+ `num_iterations`: returns the number of iterations for SMC.
Note that the class should accept reference to observations through the constructor as the observations are neccessary for computing the log weight of the samples. 

## SMCOptions
Instantiate this class and set various options.
```cpp
#include "smc_options.hpp"

SMCOptions smc_options;
smc_options.num_particles = 10000;
smc_options.ess_threshold = 1.0;
smc_options.main_seed = 532;
smc_options.resampling_seed = 8234532;
smc_options.track_population = true;
smc_options.init();
```
Note: `init` function generates random object for the seend provided. This class may change over time.

## Running SMC
Assuming that the parameters are known, we can infer the latent states using `SMC`. It offers features to get log marginal likelihood estimate in addition to retrieving the particle population at each iteration of the SMC.
```cpp
DiscreteHMMParams
ProblemSpecification<int, DiscreteHMMParams> *smc_model = new DiscreteHMM(num_latent_states, y);
SMC<int, DiscreteHMMParams> smc(smc_model, smc_options);
smc.run_smc(params);
double log_marginal = smc.get_log_marginal_likelihood();
ParticlePopluation<int> *pop0 = smc.get_population(0); // returns the particle population for X(0)
```
ParticlePopulation offers functions for getting the state values and log weights.

## PMCMCOptions
To be completed.

## Running PMCMC
To be completed.
