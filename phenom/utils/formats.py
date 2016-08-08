from numpy import arange
from scipy.interpolate import interp1d

class MakeWaveformSeries(object):
    """docstring for MakeWaveformSeries"""
    def __init__(self, freqs, hptilde, hctilde, df=None):
        """
        assuming uniformly spaced!
        """
        self.f_min = freqs[0]
        self.f_max = freqs[-1]

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
            self.flist_Hz = arange(self.f_min, self.f_max, self.df)
            self.hptilde = interp_hptilde(self.flist_Hz)
            self.hctilde = interp_hctilde(self.flist_Hz)
            self.npts = len(self.flist_Hz)

            
