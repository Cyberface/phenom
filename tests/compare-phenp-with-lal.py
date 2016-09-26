import lal
import lalsimulation as lalsim

m1 = 3.
m2 = 9.
chi1x = 0.4
chi1y = 0.02
chi1z = 0.1
chi2x = 0.6
chi2y = 0.111
chi2z = 0.2
f_min = 10.
fRef = 0.
inclination = 0.9


print "\n\n"
print "LAL"
print "\n\n"

_lshp, _lshc = lalsim.SimInspiralChooseFDWaveform(0, 1./64.,
                                                m1*lal.MSUN_SI, m2*lal.MSUN_SI,
                                                chi1x, chi1y, chi1z,
                                                chi2x, chi2y, chi2z,
                                                f_min, 0, fRef,
                                                1e6*lal.PC_SI,
                                                inclination,
                                                0, 0,
                                                None, None,
                                                -1, -1,
                                                lalsim.IMRPhenomPv2)

input_params = {}
input_params.update({'m1' : m1})
input_params.update({'m2' : m2})
input_params.update({'chi1x' : chi1x})
input_params.update({'chi1y' : chi1y})
input_params.update({'chi1z' : chi1z})
input_params.update({'chi2x' : chi2x})
input_params.update({'chi2y' : chi2y})
input_params.update({'chi2z' : chi2z})
input_params.update({'f_min' : f_min})
input_params.update({'fRef' : fRef})
input_params.update({'inclination' : inclination})

print "\n\n"
print "phenom"
print "\n\n"


import phenom
phenp = phenom.PhenomP(**input_params)

# lsf = (np.arange(_lshp.data.length) * _lshp.deltaF)
# lshp = (_lshp.data.data)
# lshc = (_lshc.data.data)
