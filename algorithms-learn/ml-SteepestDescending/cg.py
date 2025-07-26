import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def load_data() -> tuple[np.ndarray, np.ndarray]:
    '''data available at https://archive.ics.uci.edu/dataset/9/auto+mpg.

    returns test data for linear fitting for least square method fitting, gradient descending fitting, and conjugated gradient fitting

    Parameters
    ----------
    no parameters

    Returns
    -------
    X: 'mpg' data

    y: 'acceleration' data
    '''

    head = ['mpg', 'cylinders', 'displacement', 'horsepower', 'weight', 'acceleration', 'model year', 'origin', 'car name']
    df = pd.read_csv("auto-mpg.data", names=head, sep=r"\s+")
    X = df['mpg'].to_numpy()
    y = df['acceleration'].to_numpy()

    return X, y



def lsm(X:np.ndarray, y:np.ndarray) -> tuple[float, float, float, float]:
    '''least square method linear fitting benchmark, y = kx + b
        
    Parameters
    ----------
    X: 1D x data

    y: 1D y data, usually y and x have linear relation

    
    Returns
    -------
    beta0: intercept

    beta1: slope

    loss: loss function 

    r2: R^2
    '''

    X = np.hstack((np.ones(len(X))[:,None], X[:,None])) # add intercept offset 
    cov = X.T @ X                                       # covariance matrix
    beta = np.dot(np.linalg.inv(cov), np.dot(X.T, y))   # intercept and slope
    
    y_pred = (X @ beta[:,None]).flatten()

    loss = np.sum((y_pred - y)**2) / len(X)             # loss function
    r2 = 1. - np.sum((y - y_pred)**2) / np.sum((y - y.mean())**2) # R^2

    return *beta, loss, r2


def gd(X:np.ndarray, y:np.ndarray) -> tuple[float, float, list[float]]:
    '''gradient descending, steepest descending

    Parameters
    ----------
    X: 1D x data

    y: 1D y data

    Returns
    -------
    beta0: intercept

    beta1: slope

    losses: loss function changes during optimization
    '''

    X_mean, X_std = X.mean(), X.std()
    y_mean, y_std = y.mean(), y.std()

    X = (X - X_mean) / X_std                            # centered and standardized
    X = np.hstack((np.ones(len(X))[:,None], X[:,None])) # add offset column
    y = (y - y_mean) / y_std                            # centered and standardized

    initial = np.array([1., 1.])
    alpha = 0.2
    maxiter = 2000
    tol = 1e-10
    counter = 0
    loss = float('inf')
    losses = []
    nsample = len(X)

    while True:
        grad = 2. * (X.T @ (X @ initial[:,None] - y[:,None])).flatten() / nsample
        initial -= alpha * grad

        y_pred = (X @ initial).flatten()
        loss_new = np.sum((y_pred - y)**2) / nsample
        losses.append(loss_new)
        diff = loss_new - loss

        if abs(diff) <= tol:
            print(f"\nConverged after {counter} steps, loss {loss_new:.6f}, gradient {grad}")
            break
        loss = loss_new
        # print(f"Iteration {counter:<4}, loss {loss_new:<10}, difference {diff:.8e}")
        counter += 1
        if counter >= maxiter:
            break

    beta1 = (y_std / X_std) * initial[1]
    beta0 = initial[0] * y_std + y_mean - beta1 * X_mean
    print(f"beta0: {beta0:.6f}", f"beta1: {beta1:.6f}")

    return beta0, beta1, losses
    


def cg(X:np.ndarray, y:np.ndarray) -> tuple[float, float, list[float]]:
    '''conjugate gradient. 

    CG solves linear equation Ax = b, A should be positive definite matrix. 

    For linear fitting A = X^T X, A is symmetry, and X^T y = b

    Attention! This function is not generic!

    Parameters
    ----------
    X: 1D x data

    y: 1D y data

    Returns
    -------
    beta0: intecept

    beta1: slope

    losses: loss function during optimization

    Examples
    --------
    >>> beta = cg()
    >>> print(beta)
    1.0 2.0
    '''

    X_mean, X_std = X.mean(), X.std()
    y_mean, y_std = y.mean(), y.std()
    X = (X - X_mean) / X_std
    X = np.hstack((np.ones(len(X))[:, None], X[:, None]))
    y = (y - y_mean) / y_std

    n_samples = len(X)
    maxiter = 2000
    tol = 1e-10
    count = 0
    losses = []
    A = X.T @ X
    b = X.T @ y
    beta = np.array([0., 0.]) # initial value
    rold = b - A @ beta       # initial residue
    d = rold.copy()           # initial search direction

    while True:
        rsold = rold.T @ rold
        Ad = A @ d
        losses.append(rsold / n_samples)

        alpha = rsold / (d.T @ Ad)
        beta += alpha * d
        rnew = rold - alpha * Ad
        d = rnew + (rnew.T @ rnew) / rsold * d
        rold = rnew

        if np.sqrt(rsold) <= tol:
                break

        if count >= maxiter:
            break
        count += 1

    beta1 = (y_std / X_std) * beta[1]
    beta0 = beta[0] * y_std + y_mean - beta1 * X_mean

    return beta0, beta1, losses


def benchmarking():
    X, y = load_data()

    beta0, beta1, _, _ = lsm(X, y)
    beta0_sd, beta1_sd, _ = gd(X, y)
    beta0_cg, beta1_cg, losses = cg(X, y)

    print(beta0, beta1)
    print(beta0_sd, beta1_sd)
    print(beta0_cg, beta1_cg)



def solve2():
    '''Solve Ax = b. Here `A` matrix is symmetric by design. Though `A^-1 y = x` also does the job'''
    size = 5
    tmp = np.random.random((size, size))
    A = tmp + tmp.T # symmetry
    _x = np.random.normal(0., 1., (size,)) # standard solution
    y = A @ _x

    maxiter = 2000
    tol = 1e-10
    count = 0
    b = y.copy()
    x = np.zeros((size,)) # initial guess
    rold = b - A @ x   # initial residue
    d = rold.copy()       # initial search direction

    while True:
        rsold = rold.T @ rold
        Ad = A @ d

        alpha = rsold / (d.T @ Ad)
        x += alpha * d
        rnew = rold - alpha * Ad
        d = rnew + (rnew.T @ rnew) / rsold * d
        rold = rnew

        if np.sqrt(rsold) <= tol:
            print(f"Converged after {count} iterations")
            break

        if count >= maxiter:
            break
        count += 1

    print(np.allclose(x, _x))


def solve():
    '''Solve Ax = b. Here A is not symmetric and may not be positive definite'''
    size = 500
    _A = np.random.random((size, size))
    _x = np.random.normal(0., 1., (size,))
    _y = _A @ _x

    A = _A.T @ _A
    y = _A.T @ _y

    maxiter = 2000
    tol = 1e-10
    count = 0
    b = y.copy()
    x = np.zeros((size,)) # initial guess
    rold = b - A @ x   # initial residue
    d = rold.copy()       # initial search direction

    while True:
        rsold = rold.T @ rold
        Ad = A @ d

        alpha = rsold / (d.T @ Ad)
        x += alpha * d
        rnew = rold - alpha * Ad
        d = rnew + (rnew.T @ rnew) / rsold * d
        rold = rnew

        if np.sqrt(rsold) <= tol:
            print(f"Converged after {count} iterations")
            break

        if count >= maxiter:
            break
        count += 1

    print(np.allclose(x, _x))




if __name__ == '__main__':
    # benchmarking()

    solve()

