import numpy as np
import matplotlib.pyplot as plt

def main():
    '''using normal equation as common least square method, but specifically designed for L-J potential, only x^6 and x^12 are considered'''
    x = np.linspace(1., 6., 500)
    x_inv = 1. / x
    x_inv6 = x_inv * x_inv * x_inv * x_inv * x_inv * x_inv
    x_inv12 = x_inv6**2
    y = (x_inv12 - 0.8 * x_inv6) * 5.

    noise = np.random.normal(0., 0.02, x.shape)
    y_noise = y + noise

    # normal equation fitting
    X = np.hstack((x_inv6[:,None], x_inv12[:,None]))
    X_cov = X.T @ X
    Y = np.dot(X.T, y_noise)
    A, B = np.dot(np.linalg.inv(X_cov), Y)

    y_fitted = A * x_inv6 + B * x_inv12

    fig, ax = plt.subplots()
    # ax.scatter(x, y, label='LJ')
    ax.plot(x, y_noise, label='LJ + noise', c='#000080', alpha=0.5)
    ax.scatter(x, y_fitted, label='LJ fitted', s=5)

    plt.legend()
    plt.tight_layout()
    plt.show()



if __name__ == '__main__':
    main()
    
