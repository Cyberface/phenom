from numpy import linspace, angle, unwrap, absolute, exp
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
            #setup frequencies
            self.df = df
            self.npts = (self.f_max-self.f_min)/self.df
            self.flist_Hz = linspace(self.f_min, self.f_max, self.npts)
            #get the amplitude and phase of each polarisation and then
            #interpolate the amplitude and phase and recombine to make strain.
            #this is more robust than interpolating the complex strain.
            #get amp
            hp_amp = self.getamp(hptilde)
            hc_amp = self.getamp(hctilde)
            #get phase
            hp_phase = self.getphase(hptilde)
            hc_phase = self.getphase(hctilde)
            #inteprolate amp
            interp_hp_amp = interp1d(freqs, hp_amp)
            interp_hc_amp = interp1d(freqs, hc_amp)
            #interpolate phase
            interp_hp_phase = interp1d(freqs, hp_phase)
            interp_hc_phase = interp1d(freqs, hc_phase)
            #sample amplitde
            new_hp_amp = interp_hp_amp(self.flist_Hz)
            new_hc_amp = interp_hc_amp(self.flist_Hz)
            #sample phase
            new_hp_phase = interp_hp_phase(self.flist_Hz)
            new_hc_phase = interp_hc_phase(self.flist_Hz)
            #compute strain from amp and phase
            self.hptilde = self.strain(new_hp_amp, new_hp_phase)
            self.hctilde = self.strain(new_hc_amp, new_hc_phase)


    def getamp(self, z):
        """z : array of complex numbers"""
        return absolute(z)

    def getphase(self, z):
        """z : array of complex numbers"""
        return unwrap(angle(z))

    def strain(self, amp, phase):
        return amp * exp(-1.j * phase)
