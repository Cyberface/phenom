import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt

import phenom
import numpy as np

import lal
import lalsimulation as lalsim

import tdfd
from helpers import *

t={}; hp={}; hc={};
f={}; hptilde={}; hctilde={};

# t['lal'], hp['lal'], hc['lal'] = CallTDWaveform(approx="SEOBNRv3", chi1x=0., iota=0., eta=0.16, srate=2**10)
t['lal'], hp['lal'], hc['lal'] = CallTDWaveform(approx="IMRPhenomPv2", chi1x=0., iota=0., eta=0.16, srate=2**10)

# window data
start_window = phenom.planck_taper(t['lal'], 0., 2.)
hp['lal-wind'] = start_window * hp['lal']

plt.figure()
plt.plot( t['lal'], hp['lal'], label='lal' )
plt.plot( t['lal'], hp['lal-wind'], label='lal-window' )
plt.legend(loc='best')
plt.show()


f['my'], hptilde['my'] = tdfd.my_fft(t['lal'], hp['lal'])
f['my'], hptilde['my-wind'] = tdfd.my_fft(t['lal'], hp['lal-wind'])

start_window = phenom.planck_taper(f['my'], 10., 15.)
hptilde['my-wind-2'] = start_window * hptilde['my-wind']




mask = f['my'] > 0
plt.figure()
plt.plot( f['my'][mask], np.abs( hptilde['my'][mask] ), label='my - FD' )
plt.plot( f['my'][mask], np.abs( hptilde['my-wind'][mask] ), label='my - FD - wind' )
# plt.plot( f['my'][mask], np.abs( hptilde['my-wind-2'][mask] ), label='my - FD - wind - 2' )
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.show()


# t['my'], hp['my'] = tdfd.my_ifft(f['my'], hptilde['my'])
# t['my'], hp['my'] = tdfd.my_ifft(f['my'], hptilde['my-wind'] )
t['my'], hp['my'] = tdfd.my_ifft(f['my'], hptilde['my-wind-2'] )


plt.figure()
plt.plot( t['lal'] - t['lal'][peakindex(hp['lal'])], hp['lal'], label='lal' )
plt.plot( t['my'] - t['my'][peakindex(hp['my'])], np.real(hp['my']), label='myifft' )
plt.legend(loc='best')
plt.show()
