import numpy as np
from scipy.linalg import solve_banded

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
        try:
            U[:, i + 1] = solve_banded((1, 1), A_bands, U[:, i] + R)
        except:
            U.fill(np.nan)
            return U
        #U[:,i+1] = np.linalg.solve(A, U[:,i] + R)
    return U

def react(u,th):
    # reaction term
    R = th[0] * u - \
        th[1] * np.power(u, 2) + \
        th[2] * np.power(u, 3) + \
        0 * np.power(u, 4)
    return R