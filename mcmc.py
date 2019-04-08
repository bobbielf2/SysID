import numpy as np
from scipy.stats import multivariate_normal as mvn
from scipy.stats import laplace


class mcmc:

    def __init__(self, niter, data, model, isSparse):
        self.niter  = niter
        self.data   = data
        self.model  = model
        self.th     = np.zeros([niter, model.dim])  # quantity of interest (QoI)
        self.p      = np.zeros(niter)   # relative poterior probability
        self.th_mc  = np.zeros([niter, model.dim])  # store the Markov chain {theta_t,t=1,2,...}
        self.isSparse = isSparse

    '''
    compute posterior probability
    '''
    def posterior(self, th):
        mu_th = self.model.mu_th
        cov_th = self.model.cov_th
        mu_eps = self.model.mu_eps
        cov_eps = self.model.cov_eps

        # compute model prediction corresp. to guess th
        predict = self.model.modelFun(th)
        if any(predict == np.nan): # if the prediction not make sense, posterior = 0
            return 0
        # data = predict + eps, where eps is assumed normally distr'd
        epsilon = np.linalg.norm( self.data -  predict)

        if self.isSparse:  # use sparse inducing prior?
            p_th = laplace.pdf(th, loc=mu_th, scale=cov_th).prod()
        else:
            # generic prior (assume std normal distr)
            p_th = mvn.pdf(th, mean=mu_th, cov=cov_th)
        # likelihood (assume std normal distr for now)
        p_eps = mvn.pdf(epsilon, mean=mu_eps, cov=cov_eps)
        # posterior ~ likelihood * prior
        p = p_eps * p_th
        return p

    '''
    proposal algorithm
    '''
    def propose(self, th_t):
        # proposal algorithm, given the current theta_t, propose the next theta
        cov = self.model.cov_th / np.power(1.5, 2)  # std deviation for the proposal distr
        th = mvn.rvs(th_t, cov)
        return th

    '''
    Metropolis-Hastings iterations
    '''
    def metropolis_hastings(self):
        # initial guess & posterior prob
        self.th[0]  = mvn.rvs(self.model.mu_th, self.model.cov_th)
        self.p[0]   = self.posterior(self.th[0])

        th_t    = self.th[0]    # current theta, theta_t
        p_t     = self.p[0]     # current p
        for i in range(1, self.niter):

            if np.mod(i, np.floor(self.niter / 10)) == 1 or (self.niter - i < 10):
                print('Iter', i, 'out of', self.niter) # print iteration info

            self.th[i]  = self.propose(th_t)            # propose new sample based on previous
            self.p[i]   = self.posterior(self.th[i])    # calculate posterior probability
            alpha       = min([1, self.p[i] / p_t]) if p_t != 0 else 1 if self.p[i] != 0 else 0     # acceptance probability

            if np.random.rand() <= alpha:  # accept or not
                th_t = self.th[i]
                p_t = self.p[i]

            self.th_mc[i] = th_t


if __name__ == "__main__":
    from mcmc_test_cases import testCase
    case = 2
    noise_level = 0.02
    test = testCase(case, noise_level)
    th_t = test.th[0]
    result_t = test.model.modelFun(test.th[0])
    print('Test case:', case)
    print('Data:', test.data)
    print('Initial guess:', th_t)
    print('Result from the guess:', result_t)
    print('Error: |data - result| =', np.linalg.norm(test.data - result_t))