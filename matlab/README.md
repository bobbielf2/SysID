# SysID (MATLAB prototype)

Identify a system using the Markov chain Monte Carlo (MCMC) simulation.

Author: Bowei (Bobbie) Wu, 2019

## Code structure

* `main.m` - main function containing a use case example
* `buildTestCase.m` - function for building a testing `model`. Key components include:
  * `model.y`: data (See `buildModel.m` for how to define data)
  * `model.modelFun`: a model function G mapping `th`, the parameters, to G(`th`), the quantity of interest (QoI). (See `buildModel.m` for how to define a model function)
  * `model.mu_th`, `model.sig_th`: mean and standard deviation of the prior p(`th`)
  * `model.mu_eps`, `model.sig_eps`: mean and standard deviation of the likelihood p(`eps`) = p(y -G(`th`))
  * `model.propose`: proposal algorithm for MCMC
  * `model.posterior`: compute posterior distribution
* `model.qoi`: information about the Quantities of Interest (QoI), see `qoi/qoiInit.m`
  * `model.eqtype`: type of system being solved, see `solvers/Solver.m`
* `mcmc/` - folder containing functions for MCMC simulation
  * `metropolis_hastings.m` implements the [Metropolis-Hastings](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm) algorithm
  * `posterior.m` computes posterior probability (supports various likelihood functions)
  * `propose.m` proposal algorithm (supports annealing)
* `solvers/` - folder containing PDE solvers for the example systems to be learned
  * `Solvers.m` - wrapper function, directing MCMC to the desired solver based on `model.eqtype`
  * `cahnhilliard1d.m` - 1-species Cahn-Hilliard equation in 1-D. Used an implicit-explicit scheme time-stepping, evolved via a Newton solver.
  * `reactDiffuse1d.m` - 1-species reaction-diffusion equation in 1-D
  * `reactDiffuse1d2sp.m` - 2-species reaction-diffusion equation in 1-D
  * `convectReactDiffuse1d.m` - convection-reaction-diffusion equation in 1-D
* `qoi/` - folder containing different types of Quantities of Interest (QoI)
  * `QoI.m` - wrapper function, directing MCMC to the desired QoI based on `model.qoiType`
  * `sizeDistribution.m` - compute PDF of size distribution of a solution field
  * `matrixComposition.m` - compute matrix composition (i.e. average of local mins) of a solution field. Can also compute average of local max if wanted.
  * `thresholdArea.m` - compute the "total area of the field that is above a threshold". This is a functional QoI (area vs threshold).