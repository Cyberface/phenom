import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt

import phenom
import numpy as np

import lal
import lalsimulation as lalsim

import tdfd
from helpers import *



def planck_taper_new(tlist, t1, t2):
    """tlist: array of times
    t1. for t<=t1 then return 0
    t2. for t>=t2 then return 1
    else return 1./(np.exp((t2-t1)/(t-t1)+(t2-t1)/(t-t2))+1)"""

    N = len(tlist)

    tout = np.zeros( N )

    for i, t in enumerate(tlist):
        if t<=t1:
            tout[i] = 0.
        elif t>=t2:
            tout[i] = 1.
        else:
            FAC = ( (t2-t1) / (t-t1) ) + ( (t2-t1) / (t-t2) )
            tout[i] = 1. / ( np.exp( FAC ) + 1. )
    return np.asarray(tout)

def window_data(t, h, start_window=None, start_window_width=None, end_window=None, end_window_width=None):
    if start_window and start_window_width is not None:
        start_window = planck_taper_new( t, start_window, start_window + start_window_width )
        h *= start_window
    if end_window and end_window_width is not None:
        end_window = 1. - planck_taper_new( t, end_window, end_window + end_window_width )
        h *= end_window
    return h


#generate data
t={}; hp={}; hc={};

t['lal'], hp['lal'], hc['lal'] = CallTDWaveform(approx="SEOBNRv3", chi1x=0., iota=0., eta=0.16, srate=2**10)

hpnew = window_data( t['lal'], hp['lal'], start_window=0.001, start_window_width=1. )



#window by hand

startwindow = planck_taper_new( t['lal'], 0., 1. )
winhp = hp['lal'] * startwindow

plt.figure()
plt.plot( t['lal'], hp['lal'], label='lal' )
# plt.plot( t['lal'], winhp, label='win hand' )
plt.plot( t['lal'], hpnew, label='lal-cond' )
plt.legend(loc='best')
plt.show()

















#
