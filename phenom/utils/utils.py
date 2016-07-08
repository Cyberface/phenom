from numpy import sqrt, pi, absolute, ndarray, asarray

class Constants:
    # >>> import lal
    # >>> print lal.MTSUN_SI
    # 4.92549102554e-06
    MTSUN_SI = 4.92549102554e-06
    MSUN_SI = 1.9885469549614615e+30
    MRSUN_SI = 1476.6250614046494
    PC_SI = 3.085677581491367e+16

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

def _find_nearest(arr,val):
    """
    helper function for 'findindex'.
    given an array and a value, return
    the index of array that is closet to val.
    input:
        arr : numpy array
        val : float or int, value to find
    """
    if isinstance(arr, ndarray) is not True:
        arr = asarray(arr)
    idx = absolute(arr - val).argmin()
    return arr[idx]

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
            arr = arr.tolist()
            try:
                return arr.index(val)
            except:
                return _find_nearest(asarray(arr), val)
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
