import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 16}) 
import numpy as np

import pickle
import os

from . import scale


def save_model(model, filename):
    """
    model {instance of sklearn.pipeline.make_pipeline}
    filename {str}
    """
    with open(filename, 'wb') as f:
        pickle.dump(model, f, protocol=pickle.HIGHEST_PROTOCOL)
        
def load_model(filename):
    """
    returns model {instance of sklearn.pipeline.make_pipeline}
    """
    with open(filename, 'rb') as f:
        obj = pickle.load(f)
    return obj

def fit(X, y, dataname, parname, model, scaleX=False, scaleY=False, showplot=True, verbose=False):
    """
    dataname {str}: used to make output directory
    
    model {instance of sklearn.pipeline.make_pipeline} or really anything with `fit(X,y)` and `predict(X)` methods
    
    """
    
    outdir = dataname
    print(f"making output dir: {outdir}")
    os.makedirs(f"{outdir}", exist_ok=True)
    
    subdir = os.path.join(outdir, parname)
    print(f"making output sub-dir: {subdir}")
    os.makedirs(f"{subdir}", exist_ok=True)
    
    if scaleX:
        X = X.copy()
        X_scalers = scale.make_scalers(X)
        X = scale.apply_scaler(X, X_scalers)
        scale.save_scalers(X_scalers, os.path.join(subdir, "X_scalers"))
    if scaleY:
        y = y.copy()
        Y_scalers = scale.make_scalers(y)
        y = scale.apply_scaler(y, Y_scalers)
        scale.save_scalers(Y_scalers, os.path.join(subdir, "Y_scalers"))

    if verbose:
        print(f"X.shape: {X.shape}")
        print(f"y.shape: {y.shape}")

        print("X data")
        print(X)
        print("y data")
        print(y)
        
    model.fit(X, y)
    
    print("saving model")
    filename = os.path.join(subdir, "model.pickle")
    print(f"saving model: {filename}")
    save_model(model, filename)
    
    yhat = model.predict(X).reshape(-1, 1)

    if scaleY:
        y = scale.apply_inverse_scaler(y, Y_scalers)
        yhat = scale.apply_inverse_scaler(yhat, Y_scalers)

    if scaleX:
        X = scale.apply_inverse_scaler(X, X_scalers)
    
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    labels = ['eta', 'cos(theta)', 'chi']
    fig.suptitle(parname)
    for i in range(len(labels)):

        axes[i].scatter(X[:,i], y, label='training')
        axes[i].scatter(X[:,i], yhat, label='final fit', marker='x', s=100)
        axes[i].set_xlabel(labels[i])
    axes[0].legend()
#     axes[3].scatter(range(len(y)), 100*(y-yhat)/y)
    axes[3].scatter(range(len(y)), (y-yhat))
#     axes[3].set_title('% difference')
    axes[3].set_title('difference')
    
    
    axes[4].scatter(range(len(y)), y)
    axes[4].scatter(range(len(y)), yhat, marker='x', s=100)
    axes[4].set_title("data vs model")
    axes[4].set_xlabel("case num")
    
#     axes[3].scatter(range(len(y)), (y-yhat))
#     axes[3].set_title('difference')

#     axes[3].set_ylim(-30,30)

    fig.savefig(os.path.join(subdir, "fit-scatter.png"), bbox_inches='tight')
    if showplot:
        plt.show()
    plt.close()
    
    return model
