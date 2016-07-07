import phenom

print "initialising instance of PhenomD class"
ph = phenom.PhenomD(m1=10., m2=10., chi1z=0.3, chi2z=0.2, f_min=10, delta_f=1./128.)
# ph = phenom.PhenomD(m1=0.5 * MSUN, m2=0.5 * MSUN, chi1z=0.3, chi2z=0.2)
# ph = phenom.PhenomD(m1=0.5 * MSUN, m2=0.5 * MSUN, chi1z=0., chi2z=0.)
# ph = phenom.PhenomD(m1=10 * MSUN, m2=10 * MSUN, chi1z=0.3, chi2z=0.2, f_min=10)


print "physical parameters of the instance:"
print "ph.p['m1'] = ", ph.p['m1']
print "ph.p['m2'] = ", ph.p['m2']
print "ph.p['chi1z'] = ", ph.p['chi1z']
print "ph.p['chi2z'] = ", ph.p['chi2z']

print ""

print "f_min = ", ph.p['f_min']
print "f_max = ", ph.p['f_max']
print "delta_f = ", ph.p['delta_f']

print "generate htilde with 'ph.IMRPhenomDGenerateFD()'"
ph.IMRPhenomDGenerateFD()

print "get amplitude and phase from htilde with 'ph.getampandphase(htilde)'"
ph.getampandphase(ph.htilde)

print "output frequency series = ph.flist_Hz."
print "to get the frequencies in dimensionless (geometric) frequencies use"
print "phenom.HztoMf(ph.flist_Hz, ph.p['Mtot'])"


print "test plots"

import matplotlib
# matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import numpy as np
from phenom.utils.utils import HztoMf

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
fig.suptitle('amplitde (left), phase (right)')

ax[0].plot(ph.flist_Hz, ph.amp)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_ylabel('$|h(f)|$')
ax[0].set_xlabel('$f(Hz)$')

ax[1].plot(ph.flist_Hz, ph.phase)
ax[1].set_xscale('log')
ax[1].set_ylabel('$\phi(f)$')
ax[1].set_xlabel('$f(Hz)$')

fig.tight_layout()
plt.savefig('./amp-and-phase-example-plot.png')
