import numpy as np

def peakindex(x):
    """
    x : numpy array.
    returns:
    index of maximum of x
    """
    return list(np.abs(x)).index(np.max(np.abs(x)))

def invfft(f, htilde, f0, start_window_width=2., zpf=2):
    """
    f : numpy array. array of frequencies (Hz)
    htilde : numpy array. fourier domain strain
    f0 : flaot. start frequency of lower window (Hz)
    start_window_width : float. width of start window (Hz)
    zpf : int. factor to zero pad to.
    returns
    =======
    t : numpy array. times series of ifft
    h : numpy array. time domain strain (ifft of htilde)
    """
    from phenom.utils import planck_taper
    from scipy.fftpack import ifft

    # to avoid wrap around perform a phase shift
    # in a relatively arbitrary way
    phase = np.unwrap( np.angle( htilde ) )
    phase_shift = phase[-1] - phase[int(np.ceil(len(phase)/4.))]
    phase_shift = phase[-1]




    extra_cycles = 6.0
    f_min = f0
    textra = int( np.ceil( extra_cycles / f_min) )

    htilde *= np.exp( -1.j * 2. * np.pi * f * phase_shift)
    # htilde *= np.exp( -1.j * 2. * np.pi * f * textra)

    start_window = planck_taper( f, f0, f0 + start_window_width )

    htilde *= start_window

    # htilde = pad_to_pow_2(htilde, 1)

    df = f[1] - f[0]
    dt = 1.0 / ( df * len(htilde) )
    maxTime = dt * len(htilde)


    # compute ifft
    h = ifft( htilde ) / dt

    times = np.arange( 0., maxTime, dt )

    #shift times so that peak amplitude is at t = 0
    maxindex = peakindex( h )

    return times - times[maxindex], h
