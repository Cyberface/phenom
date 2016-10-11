import phenom
import numpy as np

import scipy
from scipy.fftpack import fft, fftfreq, fftshift, ifft

def fft(t, h):
    """
    t : in units of seconds
    h : in units of strain
    """
    dt = t[1] - t[0]
    N = len(t)

    htilde = scipy.fftpack.fft(h) * dt
    f = scipy.fftpack.fftfreq(N, dt)

    # mask = ( f > 0 )
    # return f[mask], htilde[mask]
    return f, htilde

def myifft(f, htilde, f0, taper_low_width):
    """
    f : in units of Hz
    f0 : float, start frequency of taper
    taper_low_width : float, width in Hz of taper
    """
    phase = np.unwrap(np.angle(htilde))
    phase_shift = (phase[0] - phase[-1])
    # phase_shift = 0.
    # phase_shift = (phase[-500] - phase[-1])
    phase_shift = phase[-1]


    extra_cycles = 6.0
    f_min = f0 + (f[1]-f[0])
    textra = int( np.ceil( extra_cycles / f_min) )

    # htilde *= np.exp( -1.j * 2. * np.pi * f * textra)
    htilde *= np.exp( -1.j * 2. * np.pi * f * phase_shift)




    win_minus = phenom.planck_taper( f, f0, f0 + taper_low_width )
    htilde *= win_minus

    N = len(f)

    df = f[1] - f[0]
    dt = 1.0 / ( df * N )

    # I am not sure what the factor of two is here...
    # FIXME: Is this factor of 2 correct or should it be somewhere else?
    h = 2 * ifft(htilde) / dt
    maxTime = dt * N
    # print("highest time (maxTime) = {0}".format(maxTime))
    # print("dt = {0}".format(dt))
    t = np.arange( 0., maxTime, dt )

    return t, h
