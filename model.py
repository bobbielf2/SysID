import numpy as np

'''
Model to be learned by MCMC
'''
class Model:
    def __init__(self, dim=1):
        self.dim = dim  # problem dimension

    def attrList(self):  # list of public attributes
        return [attr for attr in dir(self) if not attr.startswith('__')]

    def buildModel(self,
                   modelFun,
                   mu_th=None,
                   cov_th=None,
                   mu_eps=0.0,
                   cov_eps=1.0):
        if mu_th is None:
            mu_th = 0.0 if self.dim == 1 else np.zeros(self.dim)
        if cov_th is None:
            cov_th = 1.0 if self.dim == 1 else np.eye(self.dim)

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


if __name__ == "__main__":
    # test fun
    myfun = lambda th: np.sin(th)
    noisy_data = lambda th: myfun(th) + np.random.randn() * 0.02

    # build model
    model = Model()
    print('Attributes before builModel: ', model.attrList())
    model.buildModel(myfun)
    print('Attributes after builModel: ', model.attrList())
    print('Exact data sin(2) =', model.modelFun(2))
    print('Noisy data sin(2) ≈', noisy_data(2))
