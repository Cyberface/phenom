from phenom.utils.utils import setmask, findindex
import numpy as np


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
    def my_match(self, ph1, ph2, fmin=0, fmax=0, psd_fun=None):
        """computes the overlap between normalised
        templates, optimised over time and phase.
        Templates must be in frequency domain.
        ph1, ph2 = instances of PhenomD class.
        Input data must be in Hz"""
        #generate strains
        ph1.IMRPhenomDGenerateFD()
        ph2.IMRPhenomDGenerateFD()

        #first thing to do is to find the
        #common frequency range.

        # print ""
        # print "initial f"
        # print ph1.flist_Hz[0], ph1.flist_Hz[-1]
        # print ph2.flist_Hz[0], ph2.flist_Hz[-1]

        # the maximum and minimum frequency for the
        # match integral must enclosed in at least
        # one of the frequency lists of the input waveforms.
        if fmax > max(ph1.flist_Hz[-1], ph2.flist_Hz[-1]):
            raise ValueError('fmax is greater than maximum end frequency')

        if fmin < min(ph1.flist_Hz[0], ph2.flist_Hz[0]):
            raise ValueError('fmin is less than minimum starting frequency')

        # setting default values for match frequency integral
        if fmin == 0:
            fmin = max(ph1.flist_Hz[0], ph2.flist_Hz[0])
        if fmax == 0:
            fmax = min(ph1.flist_Hz[-1], ph2.flist_Hz[-1])

        #if fmin is less than the highest starting frequency
        #then set fmin equal to the highest starting frequency
        if fmin <= max(ph1.flist_Hz[0], ph2.flist_Hz[0]):
            fmin = max(ph1.flist_Hz[0], ph2.flist_Hz[0])
        #if fmax is greater than the smallest ending frequency
        #then set fmax equal to the smallest ending frequency
        if fmax >= min(ph1.flist_Hz[-1], ph2.flist_Hz[-1]):
            fmax = min(ph1.flist_Hz[-1], ph2.flist_Hz[-1])

        # print ""
        # print "fmin = ", fmin
        # print "fmax = ", fmax

        #find indeces of frequency arrays in the frequency range [fmin, fmax]
        #min
        fminIndex_1 = findindex(ph1.flist_Hz, fmin)
        fminIndex_2 = findindex(ph2.flist_Hz, fmin)
        #max
        fmaxIndex_1 = findindex(ph1.flist_Hz, fmax)
        fmaxIndex_2 = findindex(ph2.flist_Hz, fmax)
        #helper function to get masks
        mask1, _, _ = setmask( ph1.flist_Hz, fminIndex_1, fmaxIndex_1 )
        mask2, _, _ = setmask( ph2.flist_Hz, fminIndex_2, fmaxIndex_2 )
        #masked frequency arrays
        flist_Hz_1 = ph1.flist_Hz[mask1]
        flist_Hz_2 = ph2.flist_Hz[mask2]

        try:
            assert(len(flist_Hz_1) == len(flist_Hz_2))
        except:
            print ""
            print "frequency ranges"
            print "flist_Hz_1[0] = ", flist_Hz_1[0]
            print "flist_Hz_1[-1] = ", flist_Hz_1[-1]
            print "flist_Hz_2[0] = ", flist_Hz_2[0]
            print "flist_Hz_2[0] = ", flist_Hz_2[-1]
            raise ValueError('lengths not equal. Check frequency ranges above, they should be the same')
        try:
            assert((flist_Hz_1==flist_Hz_2).all())
        except:
            print ""
            print "frequency ranges"
            print "flist_Hz_1[0] = ", flist_Hz_1[0]
            print "flist_Hz_1[-1] = ", flist_Hz_1[-1]
            print "flist_Hz_2[0] = ", flist_Hz_2[0]
            print "flist_Hz_2[0] = ", flist_Hz_2[-1]
            raise ValueError('frequency arrays not equal. Check frequency ranges above, they should be the same')


        #now we have constructed a unique frequency list
        #between both input waveforms we can safely assign
        #one frequency list 'flist' for the rest of the calculation.
        flist = flist_Hz_1

        # get h1 and h2 only over fmin, fmax
        h1 = ph1.htilde[mask1]
        h2 = ph2.htilde[mask2]

        try:
            assert(len(h1) == len(h2))
        except:
            print ""
            raise ValueError('length of strains not equal')

        n = len(h1)

        #setup psd
        if psd_fun is None:
            psd = np.ones(n)
        else:
            psd = psd_fun(flist)
            #TODO: weird ratio thing...
            # psd = psd_fun((flist[-1] - flist[0])/4) / psd #EasyMatch
            # psd = psd_fun(100) / psd # michael
            # psd = 4 *(flist[-1] - flist[0]) / psd #pycbc
            psd = 1./psd #simplest choice...
        h1abs = np.abs(h1)
        h2abs = np.abs(h2)
        norm1 = np.dot(h1abs, h1abs*psd)
        norm2 = np.dot(h2abs, h2abs*psd)
        integrand = h1 * h2.conj() * psd
        zpf = 5
        # padding = np.zeros(n*zpf)
        # integrand_zp = np.concatenate([padding, integrand, padding])
        #NOTE: when the length of the integrand_zp is some powers of two then the fft can be much longer than expected
        #TODO: develop robust pad_to_pow_2 and test with different sample rates, some powers of 2 give slow results.s
        integrand_zp = pad_to_pow_2(integrand, zpf)
        csnr = np.asarray(np.fft.fft(integrand_zp)) # numpy.fft = Mma iFFT with our conventions
        return np.max(np.abs(csnr)) / np.sqrt(norm1*norm2)


def pad_to_pow_2(array, zpf):
    """given input array and zpf (zero pad factor)
    return concat(a, array, b)
    where a and b are zero padded by len(array) * zpf
    and the whole array is now the next power of 2 in length
    """
    n        = len(array)
    desired_n = n * zpf
    # print ""
    # print "old len = ", n
    # print "desired_n len = ", desired_n
    # print "current ex = ", np.log2(n)
    exponent = np.ceil(np.log2(desired_n))
    # print "exponent = ", exponent
    # if exponent % 2 == 1:
        # exponent += 1
    # print "exponent = ", exponent
    nextpow2 = int(2**exponent)
    # print "nextpow2 = ", nextpow2
    to_add   = int(np.abs(nextpow2 - n))
    # print "to_add = ", to_add
    # print "to_add/2 = ", to_add/2
    leftpad  = np.zeros(to_add/2 + 1)
    # print "len(leftpad) = ", len(leftpad)
    rightpad = np.zeros(to_add/2)
    # print "len(rightpad) = ", len(rightpad)
    # print "len(array) = ", len(array)
    ret      = np.concatenate([leftpad, array, rightpad])
    # print "new len = ", len(ret)
    # print "new exp = ", (np.log2(len(ret)))
    # print "needed len = ", 2**np.ceil(np.log2(len(ret)))
    return ret
