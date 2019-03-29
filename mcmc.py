import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import multivariate_normal as mvn
#from scipy.linalg import toeplitz
from scipy.linalg import solve_banded
from scipy import sparse
from mpl_toolkits.mplot3d import Axes3D

'''
Metropolis-Hastings iterations
'''
def metropolis_hastings(model):
    niter   = model['niter']
    th      = model['th']
    p       = model['p']

    th_t    = th[0]     # current theta, theta_t
    p_t     = p[0]      # current p
    th_T    = th.copy() # store the Markov chain {theta_t,t=1,2,...}

    for i in range(1,niter):

        if np.mod(i,np.floor(niter/10)) == 1 or (niter - i < 10):
            print('Iter',i)

        th[i]   = model['propose'](th_t)    # propose new sample based on previous
        p[i]    = model['posterior'](th[i]) # calculate posterior probability
        alpha   = min([1, p[i]/p_t])        # acceptance probability

        if np.random.rand() <= alpha:           # accept or not
            th_t = th[i]
            p_t  = p[i]

        th_T[i]  = th_t

    model['th']     = th
    model['p']      = p
    model['th_T']   = th_T
    return model

'''
test cases
'''

def testCase(test_case,noise_level):
    # Initialization
    model = {}
    model['y']          = genData(test_case,noise_level)    # data
    model['modelFun']   = lambda th: myModel(th,test_case)  # model

    if test_case == 1: # solve equation: y = sin(theta)
        niter = 9000
        model.update({
            'mu_th':    0,  # mean & std for prior
            'sig_th':   np.pi,
            'mu_eps':   0,  # mean & std for likelihood
            'sig_eps':  0.02,
            'niter':    niter,                             # num of MCMC interations
            'th':       np.zeros(niter),                              # initial theta sample
            'p':        np.zeros(niter)                   # posterior probs
        })
    elif test_case == 2:  # system id: theta = prefactors in the PDE to be identified
        niter = 3000
        model.update({
            'mu_th':    np.array([0,0]),  # mean & std for prior
            'sig_th':   np.array([1,1])*0.3,
            'mu_eps':   0,    # mean & std for likelihood
            'sig_eps':  0.5,
            'niter':    niter,                          # num of MCMC interations
            'th':       [np.zeros(2)]*niter,            # initial theta sample
            'p':        np.zeros(niter)                 # posterior probs
        })
    elif test_case == 3:  # system id: theta = prefactors in the PDE to be identified
        model.update({
            'mu_th':    np.array([0,0,0]),  # mean & std for prior
            'sig_th':   np.array([1,1,1])*0.3,
            'mu_eps':   0,  # mean & std for likelihood
            'sig_eps':  0.2,
            'niter':    5000,                             # num of MCMC interations
            'th':       [[0,0,0]],                               # initial theta sample
            'p':        [postProb(model['th'][0],model)]  # initial posterior probs
        })
    model['p'][0]       = postProb(model['th'][0], model)  # initial posterior probs
    model['propose']    = lambda th_t: genSamp(th_t,model) # proposal algorithm
    model['posterior']  = lambda th: postProb(th,model)    # posterior
    return model

'''
generate data
'''
def genData(test_case,noise_level):
    if test_case == 1:
        y = np.sin(2)
        y = y  + np.random.randn()*noise_level # additive noise
        #y = y * (1 + np.random.randn() * noise_level) # multiplicative noise
    elif test_case in [2, 3]:
        U = reactDiffuse1d(np.array([1,1,0]))
        U = U + np.random.standard_normal(U.shape)*noise_level # additive noise
        #U = U .* (1 + np.random.randn(U.shape)*noise_level); # multiplicative noise
        y = U.reshape(-1)
    return y

'''
my model functions
'''
def myModel(th,test_case):
    if test_case == 1:
        y = np.sin(th)
    elif test_case == 2:
        U = reactDiffuse1d(np.r_[th[0], th[0], th[1]])
        y = U.reshape(-1)
    elif test_case == 3:
        U = reactDiffuse1d(th)
        y = U.reshape(-1)
    return y

'''
compute posterior probability
'''
def postProb(th,model):
    y_th    = model['modelFun'](th)                 # compute y(theta) from my model
    epsilon = np.linalg.norm(model['y'] - y_th)     # y = y_th + eps, where eps is assumed normally distr'd
    # prior (assume std normal distr for now)
    p_th = mvn.pdf(th, mean = model['mu_th'], cov = np.power(model['sig_th'],2) )
    # likelihood (assume std normal distr for now)
    p_eps = mvn.pdf(epsilon, mean = model['mu_eps'], cov = np.power(model['sig_eps'], 2) )
    # posterior ~ likelihood * prior
    p = p_eps * p_th
    return p

'''
proposal algorithm
'''
def genSamp(th_t,model):
    # proposal algorithm, given the current theta_t, propose the next theta
    sig = model['sig_th']/1.5 # std deviation for the proposal distr
    th = mvn.rvs(th_t,np.power(sig,2))
    return th

'''
Reaction-diffusion equation (1D) solver
'''
def reactDiffuse1d(th):

    D = 1 # diffisivity
    L, m = 5.0, 100 # domain = [-L,L], using m subdivisions
    T, n = 0.5, 10 # time = [0,T], using n time steps
    dx = L*2/m
    dt = T/n

    U = np.zeros([m+1,n+1]) # store solutions

    # form iteration matrix
    a = D*dt/(dx*dx)
    #r = np.r_[1 + 2 * a, -a, np.zeros(m - 1)]
    # #A = toeplitz(r)
    #A[0,1] = A[m,m-1] = -2*a
    A_bands = [
        np.r_[0, -2 * a, [-a] * (m - 1)],
        np.ones(m + 1) * (1 + 2 * a),
        np.r_[[-a] * (m - 1), -2 * a, 0]
    ] # banded matrix

    x = np.linspace(-L,L,m+1)
    U[:, 0] = 0.05 * np.exp(-5 * x * x)  # initial condition

    for i in range(n):
        R = react(U[:, i], th)
        U[:, i + 1] = solve_banded((1,1), A_bands, U[:, i] + R)
        #U[:,i+1] = np.linalg.solve(A, U[:,i] + R)

    return U

def react(u,th):
    # reaction term
    R = th[0] * u - th[1] * np.power(u, 2) + th[2] * np.power(u, 3) + 0 * np.power(u, 4)
    return R




if __name__ == "__main__":
    test_case   = 1                                 # test case num
    noise_level = 0.0                               # noise (percentage)
    model       = testCase(test_case,noise_level)   # generate test case
    model       = metropolis_hastings(model)        # Metropolis-Hastings Iteration

    th = np.array(model['th'])  # all the theta tried
    p = model['p']              # posterior prob corresp to th
    i_burn  = np.s_[int(p.size/4):p.size] # burn-in
    th_T = np.array(model['th_T'])
    if test_case == 1:
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax1.scatter(th, p, s=5, c='r')
        ax2 = fig.add_subplot(122)
        ax2.plot(np.r_[i_burn],th_T[i_burn])
        plt.show()
    elif test_case == 2:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(th[:,0], th[:,1], p, s=5, c='r')
        plt.show()
