from numpy import loadtxt
from scipy.interpolate import interp1d

def read_psd_from_txt(filename):
    """
    read psd from text file
    input:
        filename : string
            path to two column file (f, PSD)
    """
    data = loadtxt(filename)
    return data

def interpolate_psd(psd):
    """
    read two column array and return interpolant
    """
    f = psd[:,0]
    y = psd[:,1]
    return interp1d(f, y)

def read_and_interp_psd_from_txt(filename):
    """convenient function
    to load and interpolate a psd in one command"""
    psd = read_psd_from_txt(filename)
    return interpolate_psd(psd)
