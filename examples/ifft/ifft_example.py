
#to profile
#python -m cProfile -s tottime ifft_example.py


import phenom

import matplotlib
# matplotlib.use('MacOSX')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from phenom.utils.utils import Constants, HztoMf, pad_to_pow_2

import lal
import lalsimulation as lalsim

from ifft import *


#
# m1 = 90./2.
# m2 = 30./2.
# chi1x = 0.9
# chi1y = 0.
# chi1z = 0.
# chi2x = 0.
# chi2y = 0.
# chi2z = 0.
# delta_f = 1./32.
# f_min = 7.
# fRef = f_min
# inclination = np.pi/3. * 0.

# m1=80.4782639
# m2=16.384655
# chi1x=0.062809065
# chi1y=0.528722703
# chi1z=-0.77006942
# chi2x=-0.102698207
# chi2y=-0.0977499112
# chi2z=-0.0815029368
# delta_f=1.0/128.
# # delta_f=1.0/256.
# # delta_f=1.0/512.
# # delta_f=1.0/1024.
# f_min=10.
# fRef = 10.
# inclination=2.85646439

m1=52.384655
m2=16.384655
chi1x=0.9
chi1y=0.
chi1z=0.
chi2x=0.
chi2y=0.
chi2z=0.
delta_f=1.0/128.
# delta_f=1.0/256.
# delta_f=1.0/512.
# delta_f=1.0/1024.
f_min = 5.
fRef = 10.
inclination=2.85646439


ph_phpLAL = phenom.Waveform(approximant='IMRPhenomPv2_LAL',
                            m1=m1,
                            m2=m2,
                            chi1x=chi1x,
                            chi1y=chi1y,
                            chi1z=chi1z,
                            chi2x=chi2x,
                            chi2y=chi2y,
                            chi2z=chi2z,
                            delta_f=delta_f,
                            f_min=f_min,
                            fRef=fRef,
                            inclination=inclination)

plt.figure()
plt.plot(ph_phpLAL.flist_Hz, np.absolute(ph_phpLAL.hptilde), label='phenom.phenp_LAL')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('./FD_amplitude.png')
plt.close()



htilde = ph_phpLAL.hptilde + 1.j * ph_phpLAL.hctilde

times, h = invfft( ph_phpLAL.flist_Hz, htilde, f0=5. )

plt.figure( figsize=( 14, 8 ) )
plt.plot(times, np.real(h))
plt.plot(times, np.abs(h))
# plt.xlim(-4, 1)
plt.xlim(-40, 1)
plt.savefig('./TD_amp.png', dpi=100)
plt.close()

plt.figure( figsize=( 14, 8 ) )
plt.plot(times, np.real(h))
plt.plot(times, np.abs(h))
plt.xlim(-0.2, 0.1)
plt.savefig('./TD_amp_zoom.png', dpi=100)
plt.close()


phenompv3 = phenom.Waveform(approximant="IMRPhenomPv3")
from copy import copy
phenpv3_1 = copy(phenompv3)

phenpv3_1.input_params['m1']=m1
phenpv3_1.input_params['m2']=m2
phenpv3_1.input_params['chi1x']=chi1x
phenpv3_1.input_params['chi1y']=chi1y
phenpv3_1.input_params['chi1z']=chi1z
phenpv3_1.input_params['chi2x']=chi2x
phenpv3_1.input_params['chi2y']=chi2y
phenpv3_1.input_params['chi2z']=chi2z
phenpv3_1.input_params['inclination']=inclination
phenpv3_1.input_params['f_min']=f_min
phenpv3_1.input_params['fRef']=fRef
phenpv3_1.input_params['delta_f']=delta_f

print("starting phenompv3 generator")

#phenomp_v3 waveform generator
phenpv3_1.phenompv3(phenpv3_1.input_params)

plt.figure()
plt.plot(phenpv3_1.flist_Hz, np.absolute(phenpv3_1.hptilde), label='phenom.v3')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('./FD_amplitude_phenpv3.png')
plt.close()

htilde = phenpv3_1.hptilde + 1.j * phenpv3_1.hctilde

print("computing phenompv3 ifft")

times, h = invfft( phenpv3_1.flist_Hz, htilde, f0=5. )

plt.figure( figsize=( 14, 8 ) )
plt.plot(times, np.real(h))
plt.plot(times, np.abs(h))
# plt.xlim(-4, 1)
plt.xlim(-40, 1)
plt.savefig('./TD_amp_phenpv3.png', dpi=100)
plt.close()

plt.figure( figsize=( 14, 8 ) )
plt.plot(times, np.real(h))
plt.plot(times, np.abs(h))
plt.xlim(-0.2, 0.1)
plt.savefig('./TD_amp_zoom_phenpv3.png', dpi=100)
plt.close()
