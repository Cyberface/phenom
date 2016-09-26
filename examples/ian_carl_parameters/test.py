import phenom

import matplotlib
# matplotlib.use('MacOSX')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from phenom.utils.utils import Constants, HztoMf


import phenom

import lal
import lalsimulation as lalsim

delta_f=1/16.
f_min=30.
f_ref = 30.

phenompv3 = phenom.Waveform(approximant="IMRPhenomPv3")
from copy import copy
phenpv3_1 = copy(phenompv3)


carl_params = {}
carl_params['a'] = {}
carl_params['b'] = {}
carl_params['c'] = {}
carl_params['d'] = {}



carl_params['a']['m1'] = 80.4782639
carl_params['a']['m2'] = 16.384655
carl_params['a']['chi1x'] = 0.062809065
carl_params['a']['chi1y'] = 0.528722703
carl_params['a']['chi1z'] = -0.77006942
carl_params['a']['chi2x'] = -0.102698207
carl_params['a']['chi2y'] = -0.0977499112
carl_params['a']['chi2z'] = -0.0815029368
carl_params['a']['thetaJ'] = 2.85646439



# carl_params['a']['Mtot'] = 65.054
# carl_params['a']['eta'] = 0.15
# carl_params['a']['chi1'] = -0.773
# carl_params['a']['chi2'] = 0.054
# carl_params['a']['chip'] = -0.161
# carl_params['a']['thetaJ'] = -0.44
# carl_params['a']['alpha0'] = -0.039

# carl_params['a']['Mtot'] = 62.748
# carl_params['a']['eta'] = 0.144
# carl_params['a']['chi1'] = -0.772
# carl_params['a']['chi2'] = -0.153
# carl_params['a']['chip'] = -0.134
# carl_params['a']['thetaJ'] = 1.084
# carl_params['a']['alpha0'] = 2.773


# carl_params['a']['m1'], carl_params['a']['m2'] = phenom.utils.m1_m2_M_eta( carl_params['a']['Mtot'], carl_params['a']['eta'] )

# carl_params['a']['chi1x'] = np.cos(carl_params['a']['alpha0']) * carl_params['a']['chip']
# carl_params['a']['chi1y'] = np.sin(carl_params['a']['alpha0']) * carl_params['a']['chip']

# carl_params['b']['Mtot'] =
# carl_params['b']['eta'] =
# carl_params['b']['chi1'] =
# carl_params['b']['chi2'] =
# carl_params['b']['chip'] =
# carl_params['b']['thetaJ'] =
# carl_params['b']['alpha0'] =

phenpv3_1.input_params['m1']=carl_params['a']['m1']
phenpv3_1.input_params['m2']=carl_params['a']['m2']
# phenpv3_1.input_params['chi1z']=carl_params['a']['chi1']
# phenpv3_1.input_params['chi2z']=carl_params['a']['chi2']
# phenpv3_1.input_params['chi1x']=carl_params['a']['chi1x']
# phenpv3_1.input_params['chi1y']=carl_params['a']['chi1y']
phenpv3_1.input_params['chi1x']=carl_params['a']['chi1x']
phenpv3_1.input_params['chi1y']=carl_params['a']['chi1y']
phenpv3_1.input_params['chi1z']=carl_params['a']['chi1z']
phenpv3_1.input_params['chi2x']=carl_params['a']['chi2x']
phenpv3_1.input_params['chi2y']=carl_params['a']['chi2y']
phenpv3_1.input_params['chi2z']=carl_params['a']['chi2z']

phenpv3_1.input_params['inclination']=carl_params['a']['thetaJ']
phenpv3_1.input_params['f_min']=f_min
phenpv3_1.input_params['delta_f']=delta_f
phenpv3_1.input_params['fRef']=f_ref
#phenomp_v3 waveform generator
phenpv3_1.phenompv3(phenpv3_1.input_params)

# print dir(phenpv3_1)

ph_phpLAL = phenom.Waveform(approximant='IMRPhenomPv2_LAL',
                            m1=carl_params['a']['m1'],
                            m2=carl_params['a']['m2'],
                            chi1x=carl_params['a']['chi1x'],
                            chi1y=carl_params['a']['chi1y'],
                            chi1z=carl_params['a']['chi1z'],
                            chi2x=carl_params['a']['chi2x'],
                            chi2y=carl_params['a']['chi2y'],
                            chi2z=carl_params['a']['chi2z'],
                            delta_f=delta_f,
                            f_min=f_min,
                            fRef=f_ref,
                            inclination=carl_params['a']['thetaJ'])

plt.figure()
plt.plot(phenpv3_1.flist_Hz, np.absolute(phenpv3_1.hptilde), label='phenom.phenpv3')
plt.plot(ph_phpLAL.flist_Hz, np.absolute(ph_phpLAL.hptilde), label='phenom.phenp_LAL')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('./test.png')
