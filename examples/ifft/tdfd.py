from phenom.utils import pad_to_pow_2, planck_taper
from scipy.fftpack import fft, fftfreq, fftshift, ifft
from numpy import arange, pi, exp

def my_fft(t, h):

    # compute frequencies
    dt = t[1] - t[0]
    N = len(h)
    f = fftfreq( N, dt )

    # compute fft
    htilde = fft( h ) * dt

    return f, htilde

def my_ifft(f, htilde):

    # compute times
    df = f[1] - f[0]
    N = len(htilde)
    dt = 1. / ( N * df )
    Tmax = N * dt

    t = arange( 0., Tmax, dt )

    # phase shift to avoid wrap around

    extra_cycles = 6.
    tshift = extra_cycles / f[1] * dt

    htilde *= exp( -1.j * 2. * pi * df *  tshift )

    # compute ifft
    h = ifft( htilde ) / dt

    return t, h
