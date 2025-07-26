import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.datasets import load_digits


class PCA_:
    '''
        A homemade class that does PCA, reproduce sklearn.decomposition.PCA
        do PCA on 2D array
        nrow variables, ncol features
        ncomponent <= nfeature

        Parameters
        ----------
            same as attributes

        Attributes
        ----------
            data: 2D np.ndarray, each row represents a variable, each col represents a feature
            nvar: nrow of data
            nfeat: ncol of data
            ncomponent: how many features to keep, <= nfeat
        
        Methods
        -------
            center: center input data, PCA applies on centered data
            cov: covariance matrix of centered input data
            cov2: covariance matrix of centered input data
            pca: eig decompose covariance matrix to perform PCA
            pca2: svd to perform PCA    
        
    '''
    __slots__ = ["data","nvar","nfeat","ncomponent"]

    def __init__(self, data:np.ndarray, nvar:int, nfeat:int, ncomponent:int):
        if nfeat < ncomponent:
            raise ValueError("ncomponent should not be larger than nfeat")
        self.data:np.ndarray = data
        self.nvar:int = nvar
        self.nfeat:int = nfeat
        self.ncomponent:int = ncomponent

    def center(self)->np.ndarray:
        '''centered input'''
        mean = np.mean(self.data,axis=0)
        res = self.data - mean[None,:]
        return res    
    
    def cov(self)->np.ndarray:
        '''covariance mat, <=> np.cov(data, rowvar=False)'''
        data_centered = self.center()
        res = data_centered.T @ data_centered / (self.nvar - 1.)
        return res
    
    def cov2(self)->np.ndarray:
        '''<=> cov'''
        res = np.zeros((self.nfeat, self.nfeat))
        centered_arr = self.center()
        for i in range(self.nfeat):
            for j in range(self.nfeat):
                res[i,j] = np.sum(centered_arr[:,i] * centered_arr[:,j]) / (self.nvar - 1)
        return res

    def void(self)->None:
        cov_mat = self.cov()
        vals, vecs = np.linalg.eigh(cov_mat)
        mask = np.argsort(vals)
        vals = vals[mask]
        tmp = vals[-1:-self.ncomponent-1:-1]
        print(tmp)

    def pca(self)->np.ndarray:
        '''svd on centered input'''
        means = np.mean(self.data, axis=0)
        U, sigma, Vh = np.linalg.svd(self.center(), full_matrices=False)
        if self.ncomponent == 1:
            axis = Vh[0]
            res = (self.center() @ axis[:,None]) * axis
            for i in range(len(means)):
                res[:,i] += means[i]
        else:
            axis = Vh[:self.ncomponent]
            res = (self.center() @ axis.T) @ axis
            for i in range(len(means)):
                res[:,i] += means[i]
        return res
    
    def pca2(self)->np.ndarray:
        '''or eig decompose covariance matrix of centered input'''
        means = np.mean(self.data, axis=0)
        vals, vecs = np.linalg.eigh(self.cov())
        mask = np.argsort(vals)
        vals = vals[mask]
        vecs = vecs[:,mask]
        if self.ncomponent == 1:
            axis = vecs[:,-1]
            res = (self.center() @ axis[:,None]) * axis
            for i in range(len(means)):
                res[:,i] += means[i]
        else:
            axis = vecs[:,-self.ncomponent:]
            res = (self.center() @ axis) @ axis.T
            for i in range(len(means)):
                res[:,i] += means[i]
        return res


def linearLeastSquareMethod(data:np.ndarray)->np.ndarray:
    '''
        different from PCA
        input data should be 2D array
    '''
    nrow, ncol = data.shape
    x, y = data[:,0], data[:,1]
    X = np.zeros((nrow, ncol))
    X[:,0] = 1.
    X[:,1] = data[:,0]
    beta = np.linalg.inv(X.T @ X) @ X.T @ y[:,None]
    y_pred = X @ beta
    new_data = np.hstack((x[:,None], y_pred))
    
    return new_data




def main():
    # generate random data
    nvar = 500
    nfeat = 2
    x0 = np.linspace(0, 10, nvar, dtype=float)
    y0 = x0 + np.random.rand(nvar) * 5.2 - np.random.rand(nvar) * 5.2
    # z0 = x0 + np.random.rand(nvar) * 2.2 - np.random.rand(nvar) * 2.2
    data = np.hstack((x0[:,None], y0[:,None])) # , z0[:,None]))

    # sklearn PCA result
    pca = PCA(n_components=1, svd_solver="full")
    pca.fit(data)
    reconstructed_sklearn = pca.inverse_transform(pca.transform(data))

    # my class result
    test = PCA_(data, nvar, nfeat, 1)
    res = test.pca()
    res2 = test.pca2()
    
    # are they same
    print(np.allclose(res, res2))
    print(np.allclose(res, reconstructed_sklearn))

    # plt.figure(figsize=(6,6))
    # ax = plt.figure().add_subplot(projection='3d')
    # ax.axis('equal')
    # ax.scatter(*data.T, s=2)
    # ax.scatter(*reconstructed_sklearn.T, s=2)
    # ax.scatter(*res2.T, s=50, alpha=0.5)
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')
    # plt.tight_layout()
    # plt.show()

    plt.scatter(data[:,0], data[:,1], s=2)
    plt.scatter(res[:,0], res[:,1], s=2)
    plt.show()


def main2():
    ncomponent = 5
    raster = load_digits().data[39].reshape((8,8))
    test = PCA_(raster, 8, 8, ncomponent)
    test_res = test.pca()

    pp = PCA(n_components=ncomponent)
    pp.fit(raster)
    pp_res = pp.inverse_transform(pp.transform(raster))
    
    print(np.allclose(pp_res, test_res))

def main3():
    ncomponent = 5
    raster = load_digits().data[39].reshape((8,8))
    pp = PCA(n_components=0.9)
    pp.fit(raster)
    # pp_res = pp.inverse_transform(pp.transform(raster))
    print(pp.explained_variance_)
    # print(pp.explained_variance_ / np.sum(pp.explained_variance_) * 100)

    test = PCA_(raster, 8, 8, ncomponent)
    test.void()


def main4()->None:
    file = "raster2.dat"
    nrow = 400
    ncol = 600

    raster_data = np.loadtxt(file, dtype=float)
    rgb_data = raster_data.reshape((nrow, ncol, 3))

    noise = np.random.normal(raster_data, 1.)
    noise_rgb = noise.reshape((nrow, ncol, 3))

    ins = PCA_(noise, nrow, ncol, 3)
    raster_data_pca = ins.pca()
    rgb_data_pca = raster_data_pca.reshape((nrow, ncol, 3))
    
    fig, ax = plt.subplots(1,3)
    ax[0].imshow(noise_rgb, cmap="binary")
    ax[0].set_title("noise")

    ax[1].imshow(rgb_data_pca, cmap="binary")
    ax[1].set_title("pca processed")

    ax[2].imshow(rgb_data, cmap="binary")
    ax[2].set_title("raw")
    plt.show()


def lsm_pca_comp():
    nvar = 500
    nfeat = 2
    x0 = np.linspace(0, 10, nvar, dtype=float)
    y0 = x0 + np.random.rand(nvar) * 5.2 - np.random.rand(nvar) * 5.2
    data = np.hstack((x0[:,None], y0[:,None]))
    lsm = linearLeastSquareMethod(data)

    test = PCA_(data, nvar, nfeat, 1)
    res = test.pca()

    pp = PCA(n_components=1)
    pp.fit(data)
    ppp = pp.inverse_transform(pp.transform(data))

    fig, ax = plt.subplots(1,1)
    ax.scatter(data[:,0], data[:,1], s=2)
    ax.scatter(lsm[:,0], lsm[:,1], s=2)
    ax.scatter(res[:,0], res[:,1], s=2)
    ax.scatter(ppp[:,0], ppp[:,1], s=2)
    ax.legend(["raw","least square method","me pca", "PCA"])
    plt.show()



    
if __name__ == "__main__":
    # lsm_pca_comp()
    main()