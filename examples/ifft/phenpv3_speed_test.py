import phenom

import numpy as np
from phenom.utils.utils import Constants, HztoMf

import lal
import lalsimulation as lalsim



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

m1=80.4782639
m2=16.384655
chi1x=0.062809065
chi1y=0.528722703
chi1z=-0.77006942
chi2x=-0.102698207
chi2y=-0.0977499112
chi2z=-0.0815029368
delta_f=1.0/256.
f_min=5.
fRef = 10.
inclination=2.85646439

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
for i in range(10):
    phenpv3_1.phenompv3(phenpv3_1.input_params)

# phenpv3_1.phenompv3(phenpv3_1.input_params)
