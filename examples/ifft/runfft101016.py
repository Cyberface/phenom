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

t['lal'], hp['lal'], hc['lal'] = CallTDWaveform(approx="IMRPhenomPv2", chi1x=0., iota=0., eta=0.16, srate=2**10)

f['lal'], hptilde['lal'], hctilde['lal'] = CallFDWaveform(approx="IMRPhenomPv2", chi1x=0., iota=0., eta=0.16, srate=1./0.01, f_min=5.)

# mask = f['lal'] > 0
# plt.figure()
# plt.plot( f['lal'][mask], np.abs( hptilde['lal'][mask] ), label='lal - FD' )
# # plt.plot( f['lal'], np.abs( hptilde['lal'] ), label='lal - FD' )
# plt.xscale('log')
# plt.yscale('log')
# plt.legend(loc='best')
# plt.show()

#compute my ifft
# t['my'], hp['my'] = tdfd.conditioned_ifft(f['lal'], hptilde['lal'], start_window=10., start_window_width=4., end_window=None, end_window_width=None)
t['my'], hp['my'] = tdfd.conditioned_ifft(f['lal'], hptilde['lal'], start_window=5., start_window_width=4., end_window=None, end_window_width=None)

plt.figure()
plt.plot( t['my'] - t['my'][peakindex(hp['my'])] , np.real(hp['my']), label='my' )

# f['lal'], hptilde['lal'], hctilde['lal'] = CallFDWaveform(approx="IMRPhenomPv2", chi1x=0., iota=0., eta=0.16, srate=1./0.005)
# t['my'], hp['my'] = tdfd.conditioned_ifft(f['lal'], hptilde['lal'], start_window=10., start_window_width=11., end_window=None, end_window_width=None)
# plt.plot( t['my'] - t['my'][peakindex(hp['my'])] , np.real(hp['my']), label='more samples' )

plt.plot( t['lal'] - t['lal'][peakindex(hp['lal'])] , np.real(hp['lal']), label='lal' )


# plt.xlim(-0.05, 0.05)
plt.legend(loc='best')
plt.show()
