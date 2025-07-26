# -*- coding: utf-8 -*-
import numpy as np

class LeastSquareMethod:
    def __init__(self, order:int, x:list|np.ndarray, y:list|np.ndarray):
        self.order = order
        self.x = x
        self.y = y
        self.nmat = self.get_normal_mat()
        self.betas = self.get_coef()
        self.y_fitted = self.get_predication()
        self.r2 = self.get_r2()

    def get_normal_mat(self) -> np.ndarray:
        mat = np.zeros((len(self.x), self.order+1))
        for i in range(self.order+1):
            if i == 0:
                mat[:,0] = 1.
            else:
                mat[:,i] = self.x**i
        return mat

    def get_coef(self) -> np.ndarray:
        betas = np.linalg.inv(self.nmat.T @ self.nmat) @ self.nmat.T @ self.y
        return betas

    def get_predication(self) -> np.ndarray:
        y_fitted = np.zeros(self.y.shape)
        for beta, ord in zip(self.betas, range(self.order+1)):
            if ord == 0:
                y_fitted += beta
                continue
            y_fitted += beta * self.x**ord
        return y_fitted

    def get_r2(self) -> float:
        y_mean = self.y.mean()
        ss_res = np.sum((self.y - self.y_fitted)**2)
        ss_tot = np.sum((self.y - y_mean)**2)
        return 1. - ss_res / ss_tot


def lsm_benchmark():
    from numpy.polynomial.polynomial import Polynomial
    import matplotlib.pyplot as plt

    order = 7
    n = 200
    xmin = -5.
    xmax = 5.
    x = np.linspace(xmin, xmax, n, endpoint=True)
    # y = np.exp(-5. * x**2)
    y = 1. / (1 + 2 * x**2)
    
    lsm = LeastSquareMethod(order=order, x=x, y=y)
    y_fitted = lsm.y_fitted
    print(lsm.r2)

    z1 = np.polyfit(x, y, order)
    p1 = np.poly1d(z1)
    np_fitted = np.array([p1(ii) for ii in x])

    pp = Polynomial.fit(x, y, deg=order)
    pp_fitted = pp(x)

    print(np.allclose(y_fitted, np_fitted))
    print(np.allclose(y_fitted, pp_fitted))

    fig, ax = plt.subplots()
    ax.plot(x,y,label='data')
    ax.plot(x,y_fitted, label='fitted model')
    ax.plot(x,np_fitted, '--', label='np')
    ax.plot(x, np_fitted - y_fitted, label='diff')
    ax.legend(frameon=False)
    plt.show()




if __name__ == '__main__':
    lsm_benchmark()
