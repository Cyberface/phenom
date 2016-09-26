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

# ph_phpLAL = phenom.Waveform(approximant='IMRPhenomPv2_LAL',
#                             m1=m1,
#                             m2=m2,
#                             chi1x=chi1x,
#                             chi1y=chi1y,
#                             chi1z=chi1z,
#                             chi2x=chi2x,
#                             chi2y=chi2y,
#                             chi2z=chi2z,
#                             delta_f=delta_f,
#                             f_min=f_min,
#                             fRef=fRef,
#                             inclination=inclination)

for i in range(10):
    phenom.Waveform(approximant='IMRPhenomPv2_LAL',
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
