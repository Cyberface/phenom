import matplotlib
# matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import numpy as np

import phenom

#can call Waveform class with default values
wf1 = phenom.Waveform()

print "default_args = {0}".format(wf1.default_args)

#can call Waveform call with keyword ags
wf2 = phenom.Waveform(approximant='IMRPhenomD',m1=20., m2=150., chi1z=0.8, chi2z=1, delta_f=1/2., f_min=1)

print "wf2 args = {0}".format(wf2.input_params)

#you can still call the PhenomD class directly and generate strain
wf3 = phenom.PhenomD()
wf3.IMRPhenomDGenerateFD()



plt.plot(wf1.flist_Hz, np.absolute(wf1.htilde), label='wf1')
plt.plot(wf2.flist_Hz, np.absolute(wf2.htilde), label='wf2')
plt.plot(wf3.flist_Hz, np.absolute(wf3.htilde), label='wf3')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('test.png')
plt.close()
