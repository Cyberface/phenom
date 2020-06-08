from sklearn.preprocessing import StandardScaler
import numpy as np

def make_scalers(X):
    """
    returns a list of sklearn.preprocessing.StandardScaler
    of input matrix.

    N, M = X.shape
    N = number of entries
    M = dimensionality
    """
    N, M = X.shape

    scalers = []

    # make a scaler for each dimension
    for i in range(M):
        scalers.append(StandardScaler())
        scalers[i].fit(X[:, i].reshape(-1, 1))

    return scalers


def apply_scaler(X, scalers):
    """
    applies the scalers to X
    N, M = X.shape
    N = number of entries
    M = dimensionality
    """
    N, M = X.shape

    X_scaled = np.zeros(shape=(N, M))

    for i in range(M):
        X_scaled[:, i] = scalers[i].transform(
            X[:, i].reshape(-1, 1)).reshape(1, -1)

    return X_scaled


def apply_inverse_scaler(X_scaled, scalers):
    """
    applies the inverse scalers to X_scaled
    N, M = X.shape
    N = number of entries
    M = dimensionality
    """
    N, M = X_scaled.shape

    X = np.zeros(shape=(N, M))

    for i in range(M):
        X[:, i] = scalers[i].inverse_transform(
            X_scaled[:, i].reshape(-1, 1)).reshape(1, -1)

    return X

def save_scalers(scalers, filename):
    """
    scalers {instance of StandardScaler}
    filename {str}
    """
    np.save(filename, scalers)

def load_scalers(filename, allow_pickle=True):
    """
    returns instance of StandardScaler
    """
    scalers = np.load(filename, allow_pickle=allow_pickle)
    return scalers
    


