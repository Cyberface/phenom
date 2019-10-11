from phenom.utils.utils import setmask, findindex, pad_to_pow_2
from phenom import MakeWaveformSeries
import numpy as np



#TODO: Refactor match code
#TODO: turn the masking and selection of the frequency
#series into a method which can be access outside.
#This is very useful for debugging purposes.

class Match(object):
    """docstring for Match"""
    def __init__(self):
        # print "init Match class"
        pass

    ################################################
    # TODO: the input to this function, ph1 and ph2
    # should be members of a `FrequencySeries` class
    # with attributes flist, deltaF etc.
    def match(self, ph1, ph2, fmin=0, fmax=0, psd_fun=None, zpf=5):
        """
        input:
            ph1 : instance of phenom.MakeWaveformSeries
            ph2 : instance of phenom.MakeWaveformSeries

        computes the overlap between normalised
        templates (h_plus only!),
        optimised over time and phase.
        Templates must be in frequency domain.
        ph1, ph2 = instances of PhenomD class.
        Input data must be in Hz

        zpf : default 5 : zero pad factor
            automatically pads to power of 2
        """

        if isinstance(ph1, MakeWaveformSeries) and isinstance(ph2, MakeWaveformSeries):
            pass
        else:
            raise TypeError('ph1 and ph2 need to be instances of phenom.MakeWaveformSeries')


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
            print("")
            print("frequency spacing 1 = {0}".format(ph1.df))
            print("frequency spacing 2 = {0}".format(ph2.df))
            print("frequency spacing should be the same")
            print("")
            print("number of points 1 = {0}".format(len(ph1.flist_Hz)))
            print("number of points 2 = {0}".format(len(ph2.flist_Hz)))
            print("number of points should be the same")
            print("")
            print("frequency ranges")
            print("flist_Hz_1[0] = ", flist_Hz_1[0])
            print("flist_Hz_1[-1] = ", flist_Hz_1[-1])
            print("flist_Hz_2[0] = ", flist_Hz_2[0])
            print("flist_Hz_2[-1] = ", flist_Hz_2[-1])
            raise ValueError('lengths not equal. Check frequency ranges above, they should be the same. If they are then check the sample rate for both inputs. They should be the same.')
        try:
            tol = 6
            np.testing.assert_almost_equal(flist_Hz_1, flist_Hz_2, tol)
        except:
            raise ValueError('flist_Hz_1 and flist_Hz_2 frequency arrays not equal at the level of (number of decimal places) tolerance = {0}'.format(tol))


        #now we have constructed a unique frequency list
        #between both input waveforms we can safely assign
        #one frequency list 'flist' for the rest of the calculation.
        flist = flist_Hz_1

        # get h1 and h2 only over fmin, fmax
        try:
            h1 = ph1.hptilde[mask1]
        except:
            raise AttributeError('h1 does not have attribute hptilde')
        try:
            h2 = ph2.hptilde[mask2]
        except:
            raise AttributeError('h2 does not have attribute hptilde')

        try:
            assert(len(h1) == len(h2))
        except:
            print("")
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
        integrand_zp = pad_to_pow_2(integrand, zpf)
        csnr = np.asarray(np.fft.fft(integrand_zp)) # numpy.fft = Mma iFFT with our conventions
        return np.max(np.abs(csnr)) / np.sqrt(norm1*norm2)
