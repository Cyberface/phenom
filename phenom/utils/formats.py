from numpy import linspace
from scipy.interpolate import interp1d

class MakeWaveformSeries(object):
    """docstring for MakeWaveformSeries"""
    def __init__(self, freqs, hptilde, hctilde, f_min=None, f_max=None, df=None):
        """
        assuming uniformly spaced!
        """

        if f_min is None:
            self.f_min = freqs[0]
        else:
            self.f_min = f_min
        if f_max is None:
            self.f_max = freqs[-1]
        else:
            self.f_max = f_max



        if df is None:
            #TODO: Add option to interpolate to unitformly spaced?
            self.flist_Hz = freqs
            self.hptilde = hptilde
            self.hctilde = hctilde
            self.df = self.flist_Hz[1] - self.f_min
            self.npts = len(self.flist_Hz)
        else:
            self.df = df
            interp_hptilde = interp1d(freqs, hptilde)
            interp_hctilde = interp1d(freqs, hctilde)
            self.npts = (self.f_max-self.f_min)/self.df
            self.flist_Hz = linspace(self.f_min, self.f_max, self.npts)
            self.hptilde = interp_hptilde(self.flist_Hz)
            self.hctilde = interp_hctilde(self.flist_Hz)
