import numpy as np

'''
Model to be learned by MCMC
'''
class Model:
    def __init__(self,
                   modelFun,
                   mu_th=0.0,
                   cov_th=1.0,
                   mu_eps=0.0,
                   cov_eps=1.0):
        # problem dimension
        self.dim = 1 if isinstance(mu_th, (int, float)) else mu_th.size
        self.mu_th = mu_th  # prior mean
        self.cov_th = cov_th  # prior std
        self.mu_eps = mu_eps  # likelihood mean
        self.cov_eps = cov_eps  # likelihood std
        self.modelFun = modelFun  # model

        # save the following for the MCMC module?
        # self.th = np.zeros([self.niter,
        #                     self.dim])  # quantity of interest (QoI)
        # self.p = np.zeros(self.niter)  # relative poterior probability
        # self.data = data  # data

    def attrList(self):  # list of public attributes
        return [attr for attr in dir(self) if not attr.startswith('__')]


if __name__ == "__main__":
    # test fun
    myfun = lambda th: np.sin(th)
    noisy_data = lambda th: myfun(th) + np.random.randn() * 0.02

    # build model
    model = Model(myfun)
    print('Attributes of Model: ', model.attrList())
    print('Exact data sin(2) =', model.modelFun(2))
    print('Noisy data sin(2) ≈', noisy_data(2))
