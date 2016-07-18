
class MakeWaveformSeries(object):
    """docstring for MakeWaveformSeries"""
    def __init__(self, freqs, hptilde, hctilde):
        """
        assuming uniformly spaced!
        """
        #TODO: Add option to interpolate to unitformly spaced?
        self.flist_Hz = freqs
        self.hptilde = hptilde
        self.hctilde = hctilde
        self.f_min = self.flist_Hz[0]
        self.f_max = self.flist_Hz[-1]
        self.df = self.flist_Hz[1] - self.f_min
        self.npts = len(self.flist_Hz)


#example
#php = phenom.PhenomP()
#php_ws = MakeWaveformSeries(php.flist_Hz, php.hptilde, php.hctilde)

#external data example

# nr_wf = MakeWaveformSeries(nr_frequs, nr_amp, nr_phase)
# nr_wf = MakeWaveformSeries(nr_frequs, nr_hptilde, nr_hctilde)

# nr_wf.hptilde
