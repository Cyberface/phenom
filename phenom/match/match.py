from phenom.utils.utils import setmask, findindex



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

    ################################################
    # TODO: the input to this function, ph1 and ph2
    # should be members of a `FrequencySeries` class
    # with attributes flist, deltaF etc.
    def my_match(self, ph1, ph2, fmin, fmax):
        """computes the overlap between normalised
        templates, optimised over time and phase.
        Templates must be in frequency domain.
        ph1, ph2 = instances of PhenomD class.
        Input data must be in Hz"""
        import numpy as np
        #generate strains
        ph1.IMRPhenomDGenerateFD()
        ph2.IMRPhenomDGenerateFD()

        #first thing to do is to find the
        #common frequency range.

        if fmin == 0:
            fmin = max(ph1.flist_Hz[0], ph2.flist_Hz[0])
        if fmax == 0:
            fmax = min(ph1.flist_Hz[-1], ph2.flist_Hz[-1])

        fminIndex_1 = findindex(ph1.flist_Hz, fmin)
        fminIndex_2 = findindex(ph2.flist_Hz, fmin)

        fmaxIndex_1 = findindex(ph1.flist_Hz, fmax)
        fmaxIndex_2 = findindex(ph2.flist_Hz, fmax)

        mask1, _, _ = setmask(ph1.flist_Hz, fminIndex_1, fmaxIndex_1 )
        mask2, _, _ = setmask(ph1.flist_Hz, fminIndex_2, fmaxIndex_2 )

        flist_Hz_1 = ph1.flist_Hz[mask1]
        flist_Hz_2 = ph2.flist_Hz[mask2]

        assert(len(flist_Hz_1) == len(flist_Hz_2))
        assert((flist_Hz_1==flist_Hz_2).all())

        #now we have constructed a unique frequency list
        #between both input waveforms we can safely assign
        #one frequency list 'flist' for the rest of the calculation.
        flist = flist_Hz_1

        # get h1 and h2 only over fmin, fmax
        h1 = ph1.htilde[mask1]
        h2 = ph2.htilde[mask2]

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
