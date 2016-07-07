



class Match(object):
    """docstring for Match"""
    def __init__(self):
        print "init Match class"
        pass

    def match(self, ph1, ph2):
        """computes the overlap between normalised
        templates, optimised over time and phase.
        Templates must be in frequency domain.
        ph1, ph2 = instances of PhenomD class"""
        import numpy as np
        ph1.IMRPhenomDGenerateFD()
        ph2.IMRPhenomDGenerateFD()
        flist = ph1.flist_Hz
        h1 = ph1.htilde
        h2 = ph2.htilde
        assert(len(h1) == len(h2))
        n = len(h1)
        h1abs = np.abs(h1)
        h2abs = np.abs(h2)
        norm1 = np.dot(h1abs, h1abs)
        norm2 = np.dot(h2abs, h2abs)
        intergrad = h1 * h2.conj()
        zpf = 5
        integrand_zp = np.concatenate([np.zeros(n*zpf), intergrad, np.zeros(n*zpf)])
        csnr = np.asarray(np.fft.fft(integrand_zp)) # numpy.fft = Mma iFFT with our conventions
        return np.max(np.abs(csnr)) / np.sqrt(norm1*norm2)
