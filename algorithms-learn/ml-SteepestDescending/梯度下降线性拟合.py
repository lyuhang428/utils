import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def lsm() -> tuple[float, float, float, float]:
    # least square method fitting benchmark
    head = ['mpg', 'cylinders', 'displacement', 'horsepower', 'weight', 'acceleration', 'model year', 'origin', 'car name']
    df = pd.read_csv("auto-mpg.data", names=head, sep=r"\s+")
    X = df['mpg'].to_numpy()
    X = np.hstack((np.ones(len(X))[:,None], X[:,None]))
    
    y = df['acceleration'].to_numpy()
    
    cov = X.T @ X
    beta = np.dot(np.linalg.inv(cov), np.dot(X.T, y))
    
    y_pred = (X @ beta[:,None]).flatten()

    loss = np.sum((y_pred - y)**2) / len(X) # loss function
    r2 = 1. - np.sum((y - y_pred)**2) / np.sum((y - y.mean())**2)

    return *beta, loss, r2


def gdm() -> tuple[float, float, list[float]]:
    # gradient descending, steepest descending
    head = ['mpg', 'cylinders', 'displacement', 'horsepower', 'weight', 'acceleration', 'model year', 'origin', 'car name']
    df = pd.read_csv("auto-mpg.data", names=head, sep=r"\s+")

    X = df['mpg'].to_numpy()
    y = df['acceleration'].to_numpy()

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
    


if __name__ == '__main__':
    head = ['mpg', 'cylinders', 'displacement', 'horsepower', 'weight', 'acceleration', 'model year', 'origin', 'car name']
    df = pd.read_csv("auto-mpg.data", names=head, sep=r"\s+")

    X = df['mpg'].to_numpy()
    y = df['acceleration'].to_numpy()

    beta = [0., 0.]
    beta[0], beta[1], loss, r2 = lsm()
    
    beta0_, beta1_, losses = gdm()

    fig, ax = plt.subplots(1,2, figsize=(16,6))
    ax[0].scatter(X, y, s=2, label='raw data')
    ax[0].plot(X, X * beta[1] + beta[0], label=f'least square method {beta[0]:.6f}, {beta[1]:.6f}', c='#000080', linewidth=1)
    ax[0].plot(X, X * beta1_ + beta0_, label=f'gradient descending {beta0_:.6f}, {beta1_:.6f}', c='green', linewidth=1)
    ax[1].plot(losses, label='loss function during training', c='black')
    ax[1].legend()
    ax[0].legend()
    plt.tight_layout()
    plt.show()

