import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt

import phenom
import numpy as np

import lal
import lalsimulation as lalsim

import myfft as fft
from helpers import *

# t, hp, hc = CallTDWaveform(approx="SEOBNRv2")
# t, hp, hc = CallTDWaveform(approx="IMRPhenomD")
# t, hp, hc = CallTDWaveform(approx="IMRPhenomPv2", chi1x=0.5, iota=np.pi/3., eta=0.16)
# t, hp, hc = CallTDWaveform(approx="SEOBNRv3", chi1x=0.5, iota=np.pi/3., eta=0.16)

m1=80.4782639
m2=16.384655
M, eta = phenom.M_eta_m1_m2(m1, m2)


t, hp, hc = CallTDWaveform(approx="IMRPhenomPv2",
                            M=M, eta=eta,
                            chi1x=0.5)


# plt.figure()
# plt.plot(t - t[peakindex(hp)], np.real(hp), label='IMRPhenomD')
# plt.legend(loc='best')
# plt.show()

ptaper_lower = phenom.planck_taper(t, t[0], t[0] + 1.)
hp *= ptaper_lower

f, hptilde = fft.fft(t, hp)


# plt.figure()
# mask = f > 0
# plt.plot( f[mask], np.abs(hptilde[mask]), label='TDtoFD' )
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# ifft and compare

###
#

phenompv3 = phenom.Waveform(approximant="IMRPhenomPv3")
from copy import copy
phenpv3_1 = copy(phenompv3)

phenpv3_1.input_params['m1']=80.4782639
phenpv3_1.input_params['m2']=16.384655
# phenpv3_1.input_params['chi1x']=0.062809065
# phenpv3_1.input_params['chi1y']=0.528722703
# phenpv3_1.input_params['chi1z']=-0.77006942
# phenpv3_1.input_params['chi2x']=-0.102698207
# phenpv3_1.input_params['chi2y']=-0.0977499112
# phenpv3_1.input_params['chi2z']=-0.0815029368
# phenpv3_1.input_params['inclination']=2.85646439
phenpv3_1.input_params['chi1x']=0.5
phenpv3_1.input_params['chi1y']=0.
phenpv3_1.input_params['chi1z']=0.
phenpv3_1.input_params['chi2x']=0.
phenpv3_1.input_params['chi2y']=0.
phenpv3_1.input_params['chi2z']=0.
phenpv3_1.input_params['inclination']=0.
phenpv3_1.input_params['f_min']=10.
phenpv3_1.input_params['delta_f']=1.0/256.

#phenomp_v3 waveform generator
phenpv3_1.phenompv3(phenpv3_1.input_params)

f = phenpv3_1.flist_Hz
hptilde = phenpv3_1.hptilde
ptaper_lower = phenom.planck_taper(f, f[0], f[0] + 1.)
hptilde *= ptaper_lower


# plt.figure()
# mask = f > 0
# plt.plot( f[mask], np.abs(hptilde[mask]), label='phenpv3' )
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

#
###

tnew, hpnew = fft.myifft( f, hptilde, f[0], 5. )

plt.figure()
plt.plot(tnew - tnew[peakindex(hpnew)], np.real(hpnew), label='FDtoTD')
# plt.plot(t - t[peakindex(hp)], np.real(hp), label='LALSIM')
plt.legend(loc='best')
plt.xlim(-2,2)
plt.show()
