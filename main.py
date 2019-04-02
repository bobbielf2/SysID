import numpy as np
from mcmc_test_cases import testCase
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

case        = 1                                 # test case num
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
    # find modes
    from scipy.signal import find_peaks
    th = th.transpose()[0] # convert to 1d vector
    ind = np.argsort(th)  # sort p in ascending order of th
    th_mode = th[ind]
    p_mode = p[ind]
    loc, _  = find_peaks(p_mode)       # find modes of th based on local max of p
    ind2 = np.argsort(p_mode[loc])[::-1]  # sort modes in descending order
    th_mode = th_mode[loc][ind2]
    ax1.set_title("predicted modes (in desc order): " +
                  ("{:.2f}, " * len(th_mode)).format(*th_mode))
    ax2 = fig.add_subplot(122)
    ax2.plot(np.r_[i_burn],th_mc[i_burn])
    plt.show()
elif case == 2:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(th[:,0], th[:,1], p, s=5, c='r')
    mean_th = th_mc[i_burn, :].mean(axis=0)
    ax.set_title("predicted mean: [" +
                  ("{:.2f} " * len(mean_th)).format(*mean_th) + "]")
    plt.show()
elif case == 3:
    import matplotlib.ticker as mtick
    fig = plt.figure()
    ax1 = fig.add_subplot(221, projection='3d')
    ax1.scatter(th[i_burn, 0], th[i_burn, 1], p[i_burn], s=5, c='r')
    mean_th = th_mc[i_burn, :].mean(axis=0)
    ax1.set_title("predicted mean: [" +
                  ("{:.2f} " * len(mean_th)).format(*mean_th) + "]")
    ax1.set_xlabel(r'$\theta_1$')
    ax1.set_ylabel(r'$\theta_2$')
    ax1.zaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax2 = fig.add_subplot(222, projection='3d')
    ax2.scatter(th[i_burn, 1], th[i_burn, 2], p[i_burn], s=5, c='r')
    ax2.set_xlabel(r'$\theta_2$')
    ax2.set_ylabel(r'$\theta_3$')
    ax2.zaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax3 = fig.add_subplot(223, projection='3d')
    ax3.scatter(th_mc[i_burn,0],th_mc[i_burn,1],th_mc[i_burn,0],s=5,c='r')
    # ax3.set_aspect('equal')
    ax3.set_xlabel(r'$\theta_1$')
    ax3.set_ylabel(r'$\theta_2$')
    ax3.set_zlabel(r'$\theta_3$')
    ax4 = fig.add_subplot(224, projection=None)
    ax4.plot(th_mc) # plot Markov chain
    plt.show()
