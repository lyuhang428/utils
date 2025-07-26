# k-means is more appropriate for clusters that are isotropic and normally distributed (i.e. spherical gaussians). 
# For anisotropical data it can lead to wrong clustering. 
# sklearn KMeans scales like O(1), but this self-implentation scales linearly

import numpy as np
import matplotlib.pyplot as plt
import sklearn.cluster
from sklearn.datasets import make_blobs, make_circles, make_moons
from typing import Any


def timeit(func:callable) -> callable:
    def wrapper(*args, **kwargs) -> Any:
        import time
        begin = time.time()
        res = func(*args, **kwargs)
        print(f"{func.__name__} took {time.time()-begin:.6f} sec")
        return res
    return wrapper


class KMeans:
    def __init__(self, X:np.ndarray, k:int=3, max_iter:int=200, init:str|None='kmeans++', ninit:int=10):
        self.X = X
        self.k = k
        self.max_iter = max_iter
        self.nsample, self.nfeat = X.shape
        self.init = init
        self.ninit = ninit
        self.centroids = self.init_centroids_kmeanpp() if self.init == 'kmeans' else self.init_centroids()
        self.labels = np.zeros((self.nsample, ), dtype=int)

    def init_centroids(self) -> np.ndarray:
        rng = np.random.default_rng()
        centroids = rng.choice(self.X, self.k)
        return centroids

    def init_centroids_kmeanspp(self) -> np.ndarray:
        centroids = []
        rng = np.random.default_rng()
        init_idx = rng.integers(self.nsample)
        centroids.append(self.X[init_idx])

        for _ in range(1, self.k):
            # Compute the minimum squared distance of each point to any centroid
            dists = np.min(np.sum((self.X[:, None, :] - np.array(centroids)[None, :, :]) ** 2, axis=2), axis=1)
            probs = dists / dists.sum()
            next_idx = rng.choice(self.nsample, p=probs)
            centroids.append(self.X[next_idx])

        return np.array(centroids)
    
    def assign_cluster(self):
        '''assign labels'''
        dists = np.sum((self.X[:, None, :] - self.centroids[None, :, :]) ** 2, axis=2)
        self.labels = np.argmin(dists, axis=1)

    def update_centroids(self) -> np.ndarray:
        new_centroids = np.zeros((self.k, self.nfeat))
        for i in range(self.k):
            mask = (self.labels == i)
            new_centroids[i] = self.X[mask].mean(axis=0)
        return new_centroids
    
    def fit(self) -> "KMeans":
        best_inertia = float('inf')
        best_model:"KMeans"|None = None

        for __ in range(self.ninit):
            model = KMeans(self.X, self.k, self.max_iter, 'kmeans++', ninit=self.ninit)
            
            for _ in range(self.max_iter):
                model.assign_cluster()
                new_centroids = model.update_centroids()
                if np.allclose(new_centroids, model.centroids):
                    break
                model.centroids = new_centroids

            dists = np.sum((model.X[:,None,:] - model.centroids[None,:,:])**2, axis=2)
            model.labels = np.argmin(dists, axis=1)
            inertia = np.sum((model.X - model.centroids[model.labels])**2)
            
            if inertia <= best_inertia:
                best_inertia = inertia
                best_model = model

        return best_model

    def predict(self, x:np.ndarray) -> int:
        '''return clustering label'''
        dists = []
        for centroid in self.centroids:
            dists.append(np.sum((x - centroid)**2))
        return np.argmin(dists)




def kmeans_benchmark():
    try:
        X = np.loadtxt('gauss_blob_testdata.dat', dtype=float, comments='#')
        yy = np.loadtxt('gauss_blob_datalabel.dat', dtype=int, comments='#')
    except FileNotFoundError:
        X, yy = make_blobs(n_samples=4000, n_features=2, centers=4)  
    finally:
        np.savetxt('gauss_blob_testdata.dat', X, fmt="%.6f")
        np.savetxt('gauss_blob_datalabel.dat', yy, fmt="%d")  

    ncluster = 5
    X, yy = make_blobs(n_samples=10000, n_features=2, centers=ncluster, cluster_std=[1., 2.5, 0.5, .2, 1.1], random_state=42)
    # X, yy = make_circles(n_samples=10000)
    # X, yy = make_moons(n_samples=5000)

    k = sklearn.cluster.KMeans(n_clusters=ncluster).fit(X)
    label_ = k.labels_
    centroids_ = k.cluster_centers_

    kmeanspp = KMeans(X, k=ncluster, init=None,ninit=10).fit()
    labelspp = kmeanspp.labels

    fig, ax = plt.subplots(1,3, figsize=(16,4))
    ax[0].scatter(X[:,0], X[:,1], c=yy, label='sklearn make blob', s=1)
    ax[1].scatter(X[:,0], X[:,1], c=label_, label='sklearn KMeans', s=1)
    ax[2].scatter(X[:,0], X[:,1], c=labelspp, label='me++', cmap='tab20', s=1)
    
    ax[1].scatter(centroids_[:,0], centroids_[:,1], marker='x', c='black', s=40)
    ax[2].scatter(kmeanspp.centroids[:,0], kmeanspp.centroids[:,1], marker='x', c='black', s=40)
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.show()




if __name__ == '__main__':
    kmeans_benchmark()


