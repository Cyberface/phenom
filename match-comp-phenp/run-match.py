import phenom

import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import numpy as np
from phenom.utils.utils import Constants, HztoMf

m1 = 160.
# m2 = 8.31
m2 = 110.
# chi1x = 0.5   # chi1x = 0.8
chi1y = 0.   # chi1y = 0.
chi1z = 0.   # chi1z = 0.2
chi2x = 0.   # chi2x = 0.
chi2y = 0.   # chi2y = 0.
chi2z = 0.   # chi2z = 0.4
f_min = 10.
fRef = 10.
# delta_f = 1/64.
delta_f = 1/32.

inclination = np.pi / 8.


input_params = {}
input_params.update({'m1' : m1})
# input_params.update({'m2' : m2})
# input_params.update({'chi1x' : chi1x})
input_params.update({'chi1y' : chi1y})
input_params.update({'chi1z' : chi1z})
input_params.update({'chi2x' : chi2x})
input_params.update({'chi2y' : chi2y})
input_params.update({'chi2z' : chi2z})
input_params.update({'f_min' : f_min})
input_params.update({'fRef' : fRef})
input_params.update({'inclination' : inclination})
input_params.update({'delta_f' : delta_f})

phenp = phenom.PhenomP(**input_params)



fl_Hz = np.arange(phenp.p['f_min'], phenp.p['f_max'], phenp.p['delta_f'])
#NOTE: CONVERTED TO Mf
fl_Hz = HztoMf(fl_Hz, phenp.p['Mtot'])

php1pars = input_params
php1pars.update({'m2' : 100.})
php1pars.update({'chi1x' : 0.5})
php2pars = input_params
php2pars.update({'m2' : 100.})
php2pars.update({'chi1x' : 0.5})
php1 = phenom.PhenomP(**php1pars)
php2 = phenom.PhenomP(**php2pars)


match = phenom.Match()
print "[without psd] match = ", match.match(php1, php2, fmin=10, fmax=0, psd_fun=None)

psd_fun = phenom.read_and_interp_psd_from_txt('/Users/sebastian/git/phenom/phenom' + '/psd/data/ZERO_DET_high_P.txt')
print "[with psd] match = ", match.match(php1, php2, fmin=10, fmax=0, psd_fun=psd_fun)



#

import lal
import lalsimulation as lalsim




# _lshp, _lshc = lalsim.SimInspiralChooseFDWaveform(0, delta_f, m1*lal.MSUN_SI, m2*lal.MSUN_SI,
#                 chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, f_min, 0, fRef, 1e6*phenom.Constants.PC_SI, inclination,
#                 0, 0, None, None, -1, -1, lalsim.IMRPhenomPv2)

# _lshp, _lshc = lalsim.SimInspiralChooseFDWaveform(0, delta_f, m1*lal.MSUN_SI, m2*lal.MSUN_SI,
#                 chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, f_min, 0, fRef, 1e6*phenom.Constants.PC_SI, inclination,
#                 0, 0, None, None, -1, -1, lalsim.IMRPhenomPv2)

_lshp, _lshc = lalsim.SimInspiralChooseFDWaveform(0, delta_f, m1*lal.MSUN_SI, 100.*lal.MSUN_SI,
                0.1, chi1y, 0.8, chi2x, chi2y, chi2z, f_min, 0, fRef, 1e6*phenom.Constants.PC_SI, inclination,
                0, 0, None, None, -1, -1, lalsim.IMRPhenomPv2)


lsf = np.arange(_lshp.data.length) * _lshp.deltaF

#NOTE: CONVERTED TO Mf
lsf = HztoMf(lsf, m1+m2)


lshp = _lshp.data.data
lshc = _lshc.data.data

from scipy.interpolate import interp1d
lshp_i = interp1d(lsf, lshp)
# lshc_i = interp1d(lsf, lshc)

print "phenom frequency domain = ", fl_Hz[0], fl_Hz[-1]
print "lal frequency domain = ", lsf[0], lsf[-1]



lshp_flHz = lshp_i(fl_Hz[:-1])

phenp_i = interp1d(fl_Hz, phenp.hp)
phenp_sample = phenp_i(fl_Hz[:-1])


class make_waveform_series(object):
    def __init__(self, flist, htilde):
        self.flist_Hz = flist
        # self.htilde = htilde
        self.hp = htilde


lal_waveform = make_waveform_series(fl_Hz[:-1], lshp_flHz)
phenp_waveform = make_waveform_series(fl_Hz[:-1], phenp_sample)


plt.figure()
plt.plot(phenp_waveform.flist_Hz, np.abs(phenp_waveform.hp), label='phenp')
plt.plot(lal_waveform.flist_Hz, np.abs(lal_waveform.hp), label='lal')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.show()

plt.figure()
plt.plot(phenp_waveform.flist_Hz, np.unwrap(np.angle(phenp_waveform.hp)), label='phenp')
plt.plot(lal_waveform.flist_Hz, np.unwrap(np.angle(lal_waveform.hp)), label='lal')
plt.legend(loc='best')
plt.show()


print lal_waveform.hp


match = phenom.Match()
print "\n\n\n"
print "match: phenom vs lal"
print "[without psd] match = {0:.15f}".format(match.match(phenp_waveform, lal_waveform, fmin=10, fmax=0, psd_fun=None))

# print "[with psd] match = {0:.15f}".format(match.match(phenp_waveform, lal_waveform, fmin=10, fmax=0, psd_fun=psd_fun))
