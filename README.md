# SysID

Identify a system using the Markov chain Monte Carlo (MCMC) simulation.

Author: Bowei (Bobbie) Wu

## Code structure

* `main.py` - main function containing a use case example
* `mcmc.py` - class for MCMC simulation
  * Implements the [Metropolis-Hastings](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm) algorithm
  * Attributes:
    * `mcmc.model`: a statistical inference model
    * `mcmc.data`: data to learn from
    * `mcmc.niter`: number of MCMC iterations to implement
  * See `main.py` for the routine of implementing a MCMC model
* `model.py` - class for building statistical inference models for a problem
  * Attributes:
    * `Model.modelFun`: a model function G mapping the quantity of interest (QoI) `th` to the data it generates ​G(`th`)
    * `Model.mu_th`, `Model.cov_th`: mean and covariance of the prior p(`th`)
    * `Model.mu_eps`, `Model.cov_eps`: mean and covariance of the deviation of a given guess G(`th`) from the actual data
  * See `mcmc_test_cases.py` for how to define a model
* `mcmc_test_cases.py` - contain test cases for MCMC inference