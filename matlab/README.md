# SysID (MATLAB prototype)

Identify a system using the Markov chain Monte Carlo (MCMC) simulation.

Author: Bowei (Bobbie) Wu, 2019

## Code structure

* `main.m` - main function containing a use case example
  
* `buildTestCase.m` - function for building a testing `model`. Key components include:
  * `model.y`: data (See `buildModel.m` for how to define data and model function)
  * `model.modelFun`: a model function G mapping `th`, the parameters, to G(`th`), the quantity of interest (QoI).
  * `model.mu_th`, `model.sig_th`: mean and standard deviation of the prior p(`th`)
  * `model.mu_eps`, `model.sig_eps`: mean and standard deviation of the likelihood p(`eps`) = p(y -G(`th`))
  * `model.propose`: proposal algorithm for MCMC
  * `model.posterior`: compute posterior distribution

* `mcmc/` - folder containing functions for MCMC simulation
  * `metropolis_hastings.m` implements the [Metropolis-Hastings](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm) algorithm
  * `posterior.m` computes posterior probability (supports various likelihood functions)
  * `propose.m` proposal algorithm (supports annealing)
  
* `solvers/` - folder containing PDE solvers for the example systems to be learned
  * `reactDiffuse1d.m` - 1-species reaction-diffusion equation in 1-D
  * `reactDiffuse1d2sp.m` - 2-species reaction-diffusion equation in 1-D
  * `convectReactDiffuse1d.m` - convection-reaction-diffusion equation in 1-D