import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt

import phenom
import numpy as np

import lal
import lalsimulation as lalsim

import myfft as fft
from helpers import *


# I need to test my ifft algorithm.
# to do this I will do the following:
# 1. take phenpv2, a FD approximant.
# use LAL to compute the time domain and frequency domain data for this waveform
# 2. use my ifft algorithm to compute the ifft of the FD data from 1.
# 3. compare my TD with TD(LAL), hopefully they are the same.

# 1. compute TD and FD from LAL for phenpv2
t={}; hp={}; hc={};
f={}; hptilde={}; hctilde={};

t['lal'], hp['lal'], hc['lal'] = CallTDWaveform(approx="IMRPhenomPv2", chi1x=0.5, iota=np.pi/3., eta=0.16, srate=2**10, f_min=5., f_ref=10.)

f['lal'], hptilde['lal'], hctilde['lal'] = CallFDWaveform(approx="IMRPhenomPv2", chi1x=0.5, iota=np.pi/3., eta=0.16, srate=1./0.0017)

# plt.figure()
# plt.plot( t['lal'], np.real( hp['lal'] ), label='lal - TD' )
# plt.legend(loc='best')
# plt.show()


# mask = f['lal'] > 0
# plt.figure()
# plt.plot( f['lal'][mask], np.abs( hptilde['lal'][mask] ), label='lal - FD' )
# plt.xscale('log')
# plt.yscale('log')
# plt.legend(loc='best')
# plt.show()

# 2. use my ifft algorithm to compute the TD data from the lal FD data.

start_window = phenom.planck_taper( f['lal'], 10. , 10. + 30. )
end_window = 1. - phenom.planck_taper( f['lal'], 150., 150. + 300. )
hptilde['lal'] *= start_window * end_window

mask = f['lal'] > 0
plt.figure()
plt.plot( f['lal'][mask], np.abs( hptilde['lal'][mask] ), label='lal - FD' )
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.show()

t['my'], hp['my'] = fft.myifft( f['lal'], hptilde['lal'], f['lal'][0], 5. )

plt.figure()
plt.plot( t['lal'] - t['lal'][peakindex(hp['lal'])], np.real( hp['lal'] ), label='lal - TD' )
plt.plot( t['my'] - t['my'][peakindex(hp['my'])], np.real( hp['my'] ), label='my - TD' )
plt.legend(loc='best')
plt.show()
