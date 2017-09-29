from __future__ import division
from numpy import sqrt, pi, absolute, ndarray, asarray, concatenate, zeros, max, dot, exp, arctan2, cos, sin, arange
from numpy.linalg import norm

class Constants:
    # >>> import lal
    # >>> print lal.MTSUN_SI
    # 4.92549102554e-06
    MTSUN_SI = 4.92549102554e-06
    MSUN_SI = 1.9885469549614615e+30
    MRSUN_SI = 1476.6250614046494
    PC_SI = 3.085677581491367e+16
    LAL_PI = 3.141592653589793 #TODO: replace instances of pi in code with LAL_PI


def ceilpow2(n):
    """
    # from pyCBC this one works... mine didnt work.
    convenience function to determine a power-of-2 upper frequency limit"""
    from math import frexp
    signif,exponent = frexp(n)
    if (signif < 0):
        return 1
    if (signif == 0.5):
        exponent -= 1
    return (1) << exponent

def pad_to_pow_2(arr, zpf):
    """
    arr : array or list to pad
    zpf : zero padding factor : integer
        if 0 then returns original length of arr, rounded to nearst power of 2
    """
    n = len(arr)

    initial_pad_left = zpf * n
    initial_pad_right = zpf * n
    initial_finial_length = initial_pad_left + n + initial_pad_right

    next_power_2 = ceilpow2(initial_finial_length)

    to_add = absolute(next_power_2 - initial_finial_length)

    #if even, add to both sides equally
    #if odd, add extra one to left side
    if (to_add % 2) == 0:
        to_add_left = int(to_add / 2)
        to_add_right = int(to_add / 2)
    else:
        to_add_left = int(to_add / 2 + 1)
        to_add_right = int(to_add / 2)

    left = zeros(initial_pad_left + to_add_left)
    right = zeros(initial_pad_right + to_add_right)

    return concatenate([left, arr, right])

def setmask(arr, x1=None, x2=None):
    """setmask(arr, x1, x2)
    arr = 1D arr
    x1 = lower value
    x2 = upper value
    by default it returns a mask that
    is the full range of arr
    returns
    =======
    mask, x1, x2
    if input x1 and x2 are out of bounds
    then it sets x1 and x2 to the boundry of arr"""
    #setup up masks
    if (x1 == None) | (x1 <= arr[0]):
        #if no starting value for fit given use lowest
        #or if starting value too low, defaulting to inital value in data
        x1 = arr[0]
    if (x2 == None) | (x2 >= arr[-1]):
        #if no ending value for fit given use highest
        #or if starting value too high, defaulting to ending value in data
        x2 = arr[-1]
    #data to be fit is masked
    #data will only be fit over
    #the masked values
    mask = ( arr >= x1 ) & ( arr <= x2 )
    return mask, x1, x2

def findindex(arr, val):
    """
    given an array and a value, return
    the index of array that is closet to val.
    input:
        arr : numpy array
        val : float or int, value to find
    """
    if isinstance(arr, ndarray) is True:
        try:
            idx = absolute(arr - val).argmin()
            return arr[idx]
        except:
            raise ValueError('input arr must be either numpy array, list or tuple')

def MftoHz(Mf, M):
    """MftoHz(Mf, M)
    """
    return Mf / (Constants.MTSUN_SI*M)

def HztoMf(Hz, M):
    """HztoMf(Hz, M)
    """
    return Hz * (Constants.MTSUN_SI*M)

def StoM(S, Mtot):
    """StoM(S, Mtot)
    """
    return S / (Constants.MTSUN_SI*Mtot)

def MtoS(M, Mtot):
    """MtoS(Hz, Mtot)
    """
    return M * (Constants.MTSUN_SI*Mtot)

def Mc_m1_m2(m1, m2):
    """
    Computes the chirp mass (Mc) from the component masses
    input: m1, m2
    output: Mc
    """
    Mc = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
    return Mc

def eta_Mc_M(Mc, M):
    """
    Computes the symmetric-mass-ratio from the Chirp Mass and
    total mass
    input: Mc, M
    output: eta
    """
    return (Mc/M)**(5./3.)

def Mc_M_eta(M, eta):
    """
    Computes the chirp mass (Mc) from the total mass and
    symmetric-mass-ratio
    input: M, eta
    output: Mc
    """
    return M * eta**(3./5.)

def Mc_eta_m1_m2(m1, m2):
    """
    Computes the symmetric mass-ratio (eta)
    and chirp mass (Mc) from the component masses
    input: m1, m2
    output: Mc, eta
    """
    Mc = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
    eta = m1*m2/(m1+m2)**2
    return Mc, eta

def M_eta_m1_m2(m1, m2):
    """
    Computes the symmetric mass-ratio (eta)
    and total mass (M) from the component masses
    input: m1, m2
    output: M, eta
    """
    M = m1+m2
    eta = m1*m2/(m1+m2)**2
    return M, eta

def eta_from_q(q):
    """
    converts mass-ratio to symmetric mass-ratio
    input: q
    output: eta
    """
    return q/(1.+q)**2

def q_from_eta(eta):
    """
    Assumes m1 >= m2
    converts symmetric-mass-ratio to mass-ratio
    input: eta
    output: q
    """
    Seta = sqrt(1. - 4. * eta)
    return (1. + Seta - 2. * eta)/(2. * eta)

def m1_m2_M_eta(M, eta):
    """
    Assumes m1 >= m2
    Computes the component masses m1 and m2
    from the total mass and symmetric mass-ratio.
    input: M, eta
    output: m1, m2
    """
    Seta = sqrt(1. - 4. * eta)
    m1 = 1./2. * (M + Seta * M)
    m2 = 1./2. * (M - Seta * M)
    return m1, m2

def m1_m2_M_q(M, q):
    """
    Assumes m1 >= m2
    Computes the component masses m1 and m2
    from the total mass and mass-ratio.
    input: M, q
    output: m1, m2
    """
    m1 = M*q/(1.+q)
    m2 = M/(1.+q)
    return m1, m2

def chieffPH(m1, m2, s1, s2):
    """
    Computes the effective spin from PhenomB/C
    input: m1, m2, s1z, s2z
    output: chi effective
    """
    return (m1*s1 + m2*s2) / (m1 + m2)

def chipn(eta, chi1z, chi2z):
    """
    Computes the effective spin parameter from PhenomD
    input: eta, s1z, s2z
    output: chi effective"""
    # Convention m1 >= m2 and chi1 is the spin on m1
    delta = sqrt(1.0 - 4.0 * eta)
    chi_s = (chi1z + chi2z) / 2.0;
    chi_a = (chi1z - chi2z) / 2.0;
    return chi_s * (1.0 - eta * 76.0/113.0) + delta * chi_a;

def chip_fun(m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, lnhatx=0., lnhaty=0., lnhatz=1.):
    """
    Computes the chip effective precession parameter defined
    in PhysRevD.91.024043 (Below Eq.3.2).
    This function assumes that the spins are defined in a
    frame where the zhat unit vector corresponds to the
    Lnhat unit vector. As such the default values for lnhat(x,y,z)=(0.,0.,1.)

    We also assume that m1>=m2 and q>1.
    Input spins are dimensionless

    Returns:
        chip (dimensionless precession parameter)
            dimensionless spin components along the lnhat vector
        chi1_l
        chi2_l
    """
    # enforce m1 >= m2 and chi1 is on m1
    if m1<m2: # swap spins and masses
        # chi1z, chi2z = chi2z, chi1z
        chi1x, chi1y, chi1z, chi2x, chi2y, chi2z = float(chi2x), float(chi2y), float(chi2z), float(chi1x), float(chi1y), float(chi1z)
        m1, m2 = float(m2), float(m1)

    m1_2 = m1**2.
    m2_2 = m2**2.

    chi1norm = norm([chi1x, chi1y, chi1z])
    chi2norm = norm([chi2x, chi2y, chi2z])
    try:
        assert(chi1norm<=1.)
    except:
        raise AssertionError("chi1norm = {0}. chi1norm should be less than 1.: chi1x = {1}, chi1y = {2}, chi1z = {3}".format(chi1norm, chi1x, chi1y, chi1z))
    try:
        assert(chi2norm<=1.)
    except:
        raise AssertionError("chi2norm = {0}. chi2norm should be less than 1.: chi2x = {1}, chi2y = {2}, chi2z = {3}".format(chi2norm, chi2x, chi2y, chi2z))


    lnnorm = norm([lnhatx, lnhaty, lnhatz])

    tol = 1e-6
    try:
        assert(absolute(1.-lnnorm)<1e-6)
    except:
        raise AssertionError("lnnorm = {0}. lnnorm should be unit at tol = {1}: lnhatx = {2}, lnhaty = {3}, lnhatz = {4}".format(1.-lnnorm, tol, lnhatx, lnhaty, lnhatz))

    #compute the aligned spin component. The component of the spins along lnhat
    chi1_l = dot( [chi1x, chi1y, chi1z], [lnhatx, lnhaty, lnhatz] )
    chi2_l = dot( [chi2x, chi2y, chi2z], [lnhatx, lnhaty, lnhatz] )

    #compute component of spins perpendicular to lnhat
    chi1_perp_x = chi1x - chi1_l * lnhatx
    chi1_perp_y = chi1y - chi1_l * lnhaty
    chi1_perp_z = chi1z - chi1_l * lnhatz

    chi2_perp_x = chi2x - chi2_l * lnhatx
    chi2_perp_y = chi2y - chi2_l * lnhaty
    chi2_perp_z = chi2z - chi2_l * lnhatz

    #magnitude of in-plane dimensionless spins
    chi1_perp = norm( [chi1_perp_x, chi1_perp_y, chi1_perp_z] )
    chi2_perp = norm( [chi2_perp_x, chi2_perp_y, chi2_perp_z] )

    A1 = 2 + (3*m2) / (2*m1)
    A2 = 2 + (3*m1) / (2*m2)

    #magnitude of in-plane Dimensionfull spins
    S1_perp = chi1_perp * m1_2
    S2_perp = chi2_perp * m2_2

    ASp1 = A1 * S1_perp
    ASp2 = A2 * S2_perp

    num = max([ASp1, ASp2])
    den = A1*m1_2

    chip = num / den

    return chip, chi1_l, chi2_l


def amp0Func(eta):
    """
    amplitude scaling factor defined by eq. 17 in 1508.07253
    """
    return (sqrt(2.0/3.0)*sqrt(eta)) / pi**(1.0/6.0)

def pow_2_of(number):
    """
    squares input number
    helper function from lalsimulation/src/LALSimIMRPhenomD_internals.h
    """
    return number*number

def pow_3_of(number):
    """
    cube input number
    helper function from lalsimulation/src/LALSimIMRPhenomD_internals.h
    """
    return number*number*number

def pow_4_of(number):
    """
    fourth power of number
    helper function from lalsimulation/src/LALSimIMRPhenomD_internals.h
    """
    pow2 = pow_2_of(number)
    return pow2 * pow2

class UsefulPowers(object):
    """init_useful_powers from phenomD LAL code"""
    def __init__(self, number):
        self.sixth        = number ** (1.0/6.0)
        self.third        = self.sixth * self.sixth
        self.two_thirds   = number / self.third
        self.four_thirds  = number * (self.third)
        self.five_thirds  = self.four_thirds * (self.third)
        self.two          = number * number
        self.seven_thirds = self.third * self.two
        self.eight_thirds = self.two_thirds * self.two
        self.seven_sixths = number * self.sixth

def WignerdCoefficients(v, SL, eta, Sp):
    """
    v   : Cubic root of (Pi * Frequency (geometric))
    SL  : Dimensionfull aligned spin
    eta : Symmetric mass-ratio
    Sp  : Dimensionfull spin component in the orbital plane
    """
    from phenom.pn.pn import PhenomPL2PN

    #We define the shorthand s := Sp / (L + SL)
    L = PhenomPL2PN(v, eta)
    #We ignore the sign of L + SL below.
    s = Sp / (L + SL)  #s := Sp / (L + SL)
    s2 = s*s
    cos_beta = 1.0 / sqrt(1.0 + s2)
    cos_beta_half = + sqrt( (1.0 + cos_beta) / 2.0 )  #cos(beta/2)
    sin_beta_half = + sqrt( (1.0 - cos_beta) / 2.0 )  #sin(beta/2)

    return cos_beta_half, sin_beta_half

def planck_taper(tlist, t1, t2):
    """tlist: array of times
    t1. for t<=t1 then return 0
    t2. for t>=t2 then return 1
    else return 1./(np.exp((t2-t1)/(t-t1)+(t2-t1)/(t-t2))+1)"""
    tout = []
    for t in tlist:
        if t<=t1:
            tout.append(0.)
        elif t>=t2:
            tout.append(1.)
        else:
            tout.append(1./(exp((t2-t1)/(t-t1)+(t2-t1)/(t-t2))+1))
    return asarray(tout)

def velocity_to_frequency(v, M):
    return v**(3.0) / (M * Constants.MTSUN_SI * Constants.LAL_PI)

def frequency_to_velocity(f, M):
    return (Constants.LAL_PI * M * Constants.MTSUN_SI * f)**(1.0/3.0)

def f_SchwarzISCO(M):
    """
    Innermost stable circular orbit (ISCO) for a test particle 
    orbiting a Schwarzschild black hole
    Parameters
    ----------
    M : float or numpy.array
        Total mass in solar mass units
    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    return velocity_to_frequency((1.0/6.0)**(0.5), M)


