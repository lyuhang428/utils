# -*- encoding: utf-8 -*-
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt

class ExampleData:
    '''Lorentian example data'''
    data:np.ndarray = np.array([[-1, 0.038], 
                     [-0.8, 0.058], 
                     [-0.60, 0.10], 
                     [-0.4,0.20], 
                     [-0.2, 0.50], 
                     [0, 1], 
                     [0.2, 0.5], 
                     [0.4, 0.2], 
                     [0.6, 0.1], 
                     [0.8, 0.058], 
                     [1, 0.038]])
    
    def __init__(self):
        pass

    @classmethod
    def show_data(cls):
        fig, ax = plt.subplots(1, 1, figsize=(5,5))
        ax.scatter(cls.data[:,0], cls.data[:,1], facecolor='none', edgecolors='#800000', marker='d', s=100)
        ax.plot(cls.data[:,0], cls.data[:,1], linewidth=2.2, c='#000080')
        plt.tight_layout()
        plt.show()


class CubicSpline:
    def __init__(self, data: np.ndarray):
        self.data = data
        self.coef = self.get_coef()

    def get_mat(self, bc: str = "not-a-knot") -> np.ndarray:
        n = len(self.data) - 1  # number of intervals
        mat = np.zeros((4 * n, 4 * n))
        rhs = np.zeros(4 * n)

        for i in range(n):
            x = self.data[i, 0]
            mat[2 * i, 4 * i : 4 * i + 4] = [1, x, x**2, x**3]
            rhs[2 * i] = self.data[i, 1]

            x = self.data[i + 1, 0]
            mat[2 * i + 1, 4 * i : 4 * i + 4] = [1, x, x**2, x**3]
            rhs[2 * i + 1] = self.data[i + 1, 1]

        for i in range(n - 1):
            x = self.data[i + 1, 0]

            row = 2 * n + 2 * i
            mat[row, 4 * i : 4 * i + 4] = [0, 1, 2 * x, 3 * x**2]
            mat[row, 4 * (i + 1) : 4 * (i + 1) + 4] = [0, -1, -2 * x, -3 * x**2]

            row = 2 * n + 2 * i + 1
            mat[row, 4 * i : 4 * i + 4] = [0, 0, 2, 6 * x]
            mat[row, 4 * (i + 1) : 4 * (i + 1) + 4] = [0, 0, -2, -6 * x]

        if bc == "not-a-knot":
            x0 = self.data[1, 0]
            mat[-2, :4] = [0, 0, 0, 6]
            mat[-2, 4:8] = [0, 0, 0, -6]

            xn = self.data[-2, 0]
            mat[-1, -8:-4] = [0, 0, 0, 6]
            mat[-1, -4:] = [0, 0, 0, -6]
        elif bc == "natural":
            mat[-2, :4] = [0, 0, 2, 6 * self.data[0, 0]]
            mat[-1, -4:] = [0, 0, 2, 6 * self.data[-1, 0]]
        elif bc == "clamped":
            raise NotImplementedError("Clamped boundary condition not implemented")

        return mat, rhs

    def get_coef(self) -> np.ndarray:
        mat, rhs = self.get_mat()
        coef = np.linalg.solve(mat, rhs)  # numerically stabler than inv
        return coef

    def interpolate(self, x: float) -> float:
        if x < self.data[0, 0] or x > self.data[-1, 0]:
            raise ValueError("x out of bounds")
        idx = bisect_left(self.data[:, 0], x) - 1
        idx = max(0, idx)  # Ensure idx >= 0
        a, b, c, d = self.coef[4 * idx : 4 * idx + 4]
        return a + b * x + c * x**2 + d * x**3
    

if __name__ == "__main__":
    gaussian = lambda x: np.exp(-2. * x**2)
    x = np.arange(-2, 2+0.2, 0.2)
    y = np.array([gaussian(ii) for ii in x])
    data = np.hstack((x[:,None], y[:,None]))

    cc = CubicSpline(data)
    vals = np.arange(-2, 2, 0.01)
    res = []
    for val in vals:
        res.append(cc.interpolate(val))

    fig, ax = plt.subplots(1, 1, figsize=(5,5))
    ax.scatter(data[:,0], data[:,1], s=100, marker='o', facecolor='none', edgecolors='black', label="data")
    ax.scatter(vals, res, s=10, marker='d', facecolor='none', edgecolors='red', label="interpolation")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.show()



