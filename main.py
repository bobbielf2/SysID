import numpy as np
from mcmc_test_cases import testCase
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

case        = 2                                 # test case num
noise_level = 0.02                              # noise (percentage)
test        = testCase(case, noise_level)       # generate test case

# MCMC: Metropolis-Hastings Iteration
test.metropolis_hastings()

th = test.th  # all the theta tried
p = test.p    # posterior prob corresp to th
i_burn  = np.s_[int(p.size/4):p.size] # burn-in
th_mc = test.th_mc
if case == 1:
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.scatter(th, p, s=5, c='r')
    ax2 = fig.add_subplot(122)
    ax2.plot(np.r_[i_burn],th_mc[i_burn])
    plt.show()
elif case == 2:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(th[:,0], th[:,1], p, s=5, c='r')
    plt.show()
