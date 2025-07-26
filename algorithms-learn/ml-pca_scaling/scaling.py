# -*- encoding: utf-8 -*-
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_wine
import numpy as np

def scaling(data:np.ndarray)->np.ndarray:
    '''
        mean 0, stdvar 1
        if one feature is constant, its stdvar is 0 and scaling gives NAN or INF. Const feature should be ruled out first
        always scale input data before doing PCA; sklearn.decomposition.PCA centers data but not scales them
    '''
    res = data - np.mean(data, axis=0)[None,:]
    res /= np.std(res, axis=0)[None,:]
    if not np.allclose(np.std(res, axis=0), 1.):
        raise ValueError("scaling failed!")
    return res

def main():
    data = load_wine()
    xx = data.data

    # center to mean 0, and scale to stdvar 1
    scale = StandardScaler()
    scale.fit(xx)
    scaled = scale.transform(xx) # return is np.ndarray

    cc = scaling(xx)

    print(np.allclose(cc, scaled))

if __name__ == "__main__":
    main()