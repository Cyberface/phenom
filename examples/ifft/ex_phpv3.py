import phenom

import matplotlib
# matplotlib.use('MacOSX')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from phenom.utils.utils import Constants, HztoMf


import phenom


m1 = 35.
m2 = 30.
chi1x = 0.9
chi1y = 0.
chi1z = 0.
chi2x = 0.
chi2y = 0.
chi2z = 0.
delta_f = 1./8.
f_min = 30.
fRef = f_min
inclination = np.pi/3.


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
