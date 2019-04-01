import numpy as np
from react_diffuse import reactDiffuse1d
from model import Model
from mcmc import mcmc

def testFun(th, case):
    if case == 1:
        y = np.sin(th)
    elif case == 2:
        U = reactDiffuse1d(np.r_[th[0], th[0], th[1]])
        y = U.reshape(-1)
    elif case == 3:
        U = reactDiffuse1d(th)
        y = U.reshape(-1)
    return y


def genData(case, noise_level):
    if case == 1:
        y = np.sin(2)
        y = y + np.random.randn() * noise_level  # additive noise
        #y = y * (1 + np.random.randn() * noise_level) # multiplicative noise
    elif case in [2, 3]:
        U = reactDiffuse1d(np.array([1, 1, 0]))
        U = U + np.random.standard_normal(
            U.shape) * noise_level  # additive noise
        #U = U .* (1 + np.random.randn(U.shape)*noise_level); # multiplicative noise
        y = U.reshape(-1)
    return y


'''
define test cases
'''
def testCase(case, noise_level):
    if case == 1: # solve equation: y = sin(theta)
        niter   = 9000  # num of MCMC interations
        dim     = 1     # model dimensions
        mu_th   = 0     # mean & std for prior
        cov_th  = np.power(np.pi, 2)
        mu_eps  = 0     # mean & std for likelihood
        cov_eps = np.power(0.02, 2)
    elif case == 2:  # system id: theta = prefactors in the PDE to be identified
        niter   = 3000  # num of MCMC interations
        dim     = 2  # model dimensions
        mu_th   = np.array([0,0])  # mean & std for prior
        cov_th  = np.power(np.eye(dim) * 0.3, 2)
        mu_eps  = 0  # mean & std for likelihood
        cov_eps = np.power(0.5, 2)
    elif case == 3:  # system id: theta = prefactors in the PDE to be identified
        niter   = 5000  # num of MCMC interations
        dim     = 3  # model dimensions
        mu_th   = np.array([0, 0, 0])  # mean & std for prior
        cov_th  = np.power(np.eye(dim) * 0.3, 2)
        mu_eps  = 0  # mean & std for likelihood
        cov_eps = np.power(0.2, 2)

    data        = genData(case, noise_level)    # data
    modelFun    = lambda th: testFun(th, case)  # model

    # Model initialization
    model = Model(dim)
    model.buildModel(modelFun, mu_th, cov_th, mu_eps, cov_eps)
    test = mcmc(niter, data, model)

    return test

if __name__ == "__main__":
    case = 1
    noise_level = 0.02
    test = testCase(case, noise_level)
    print('Test case: ', case)
    print('Model attributes: ', test.model.attrList())
    print('Exact data sin(2) =', test.model.modelFun(2))
    print('Noisy data sin(2) ≈', test.data)