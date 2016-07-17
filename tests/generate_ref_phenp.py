import phenom
import numpy as np

m1 = 16.
m2 = 12.
chi1x = 0.6
chi1y = 0.
chi1z = 0.
chi2x = 0.
chi2y = 0.
chi2z = 0.
f_min = 30.
fRef = 30.
delta_f = 1/8.
inclination = np.pi / 8.

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
input_params.update({'delta_f' : delta_f})

phenp = phenom.PhenomP(**input_params)

amp = np.absolute(phenp.hp)
phase = np.unwrap(np.angle(phenp.hp))


f=phenp.flist_Hz
a=amp
p=phase

np.savetxt('./phenp_ref.dat', np.column_stack((f, a, p)))
