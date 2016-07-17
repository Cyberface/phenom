from __future__ import division

from phenom.utils.utils import m1_m2_M_q, eta_from_q
from numpy import pi, log, arccos, sign, sqrt

def PhenomPAlpha(omega, q, chi1x, chi1z, order=-1):
    """PhenomPAlpha(omega, q, chi1x, chi1z, order=-1)
    NOTE the total mass is assumed to be one and the results
    are in units of the total mass (which is one)
    omega = orbital angular frequency
    q = mass-ratio
    chi1x = in-plane spin component on larger body
    chi1z = aligned spin component on larger body
    order = PN order [default = -1 = all terms]
        other possible values:
        "ONE", "TWO", "THREE", "FOUR", "LOG"
        returns ONLY that partular term, NOT the summation.
    """

    mtot = 1.

    m1, m2 = m1_m2_M_q(mtot, q)
    # nu = symmetric mass-ratio (also called eta)
    nu = eta_from_q(q)
    dm = m1-m2

    #1/omega
    _TERM_ONE_ = + (-35./192.+ (5.*dm)/(64.*m1))/omega

    #1/omega**(2/3)
    _TERM_TWO_ = + ( \
    + (15.*dm*m1*chi1z)/(128.*mtot**2.*nu) \
    - (35.*m1**2.*chi1z)/(128.*mtot**2.*nu) \
    )/omega**(2./3.)


    #1/omega**(1/3)
    _TERM_THREE_ = \
    + (-5515./3072. \
    + (4555.*dm)/(7168.*m1) \
    - (515.*nu)/384. \
    - (15.*dm**2.*nu)/(256.*m1**2.) \
    + (175.*dm*nu)/(256.*m1) \
    + (15.*dm*m1**3.*chi1x**2.)/(128.*mtot**4.*nu**2.) \
    - (35.*m1**4.*chi1x**2.)/(128.*mtot**4.*nu**2.) \
    )/omega**(1./3.)

    #omega**(1/3)
    _TERM_FOUR_ = \
    + ( \
    40121485./9289728. - (27895885.*dm)/(21676032.*m1) \
    + (39695.*nu)/86016. \
    + (1615.*dm**2.*nu)/(28672.*m1**2.) \
    + (265.*dm*nu)/(14336.*m1) \
    + (955.*nu**2.)/576. \
    - (15.*dm**3.*nu**2.)/(1024.*m1**3.) \
    + (35.*dm**2.*nu**2.)/(256.*m1**2.) \
    - (2725.*dm*nu**2.)/(3072.*m1) \
    + (485.*dm*m1**3.*chi1x**2.)/(14336.*mtot**4.*nu**2.) \
    + (475.*m1**4.*chi1x**2.)/(6144.*mtot**4.*nu**2.) \
    + (15.*dm**2.*m1**2.*chi1x**2.)/(256.*mtot**4.*nu) \
    - (145.*dm*m1**3.*chi1x**2.)/(512.*mtot**4.*nu) \
    + (575.*m1**4.*chi1x**2.)/(1536.*mtot**4.*nu) \
    + (15.*dm*m1**7.*chi1x**4.)/(512.*mtot**8.*nu**4.)
    - (35.*m1**8.*chi1x**4.)/(512.*mtot**8.*nu**4.) \
    + (15.*dm*m1*pi*chi1z)/(16.*mtot**2.*nu) \
    - (35.*m1**2.*pi*chi1z)/(16.*mtot**2.*nu) \
    + (375.*dm**2.*m1**2.*chi1z**2.)/(256.*mtot**4.*nu) \
    - (1815.*dm*m1**3.*chi1z**2.)/(256.*mtot**4.*nu) \
    + (1645.*m1**4.*chi1z**2.)/(192.*mtot**4.*nu) \
    - (15.*dm*m1**7.*chi1x**2.*chi1z**2.)/(128.*mtot**8.*nu**4.) \
    + (35.*m1**8.*chi1x**2.*chi1z**2.)/(128.*mtot**8.*nu**4.) \
    )*omega**(1./3.)

    #log(omega)
    _LOG_TERMS_ = \
    - (35.*pi*log(omega))/48. \
    + (5.*dm*pi*log(omega))/(16.*m1) \
    + (5.*dm**2.*chi1z*log(omega))/(16.*mtot**2.) \
    - (5.*dm*m1*chi1z*log(omega))/(3.*mtot**2.) \
    + (2545.*m1**2.*chi1z*log(omega))/(1152.*mtot**2.) \
    - (2035.*dm*m1*chi1z*log(omega))/(21504.*mtot**2.*nu) \
    + (2995.*m1**2.*chi1z*log(omega))/(9216.*mtot**2.*nu) \
    + (5.*dm*m1**5.*chi1x**2.*chi1z*log(omega))/(128.*mtot**6.*nu**3.) \
    - (35.*m1**6.*chi1x**2.*chi1z*log(omega))/(384.*mtot**6.*nu**3.)

    if order == -1:
        RET = _TERM_ONE_ + _TERM_TWO_ + _TERM_THREE_ + _TERM_FOUR_ + _LOG_TERMS_
    elif order == "ONE":
        RET = _TERM_ONE_
    elif order == "TWO":
        RET = _TERM_TWO_
    elif order == "THREE":
        RET = _TERM_THREE_
    elif order == "FOUR":
        RET = _TERM_FOUR_
    elif order == "LOG":
        RET = _LOG_TERMS_
    elif order == "SUMONE_NOLOG":
        RET = _TERM_ONE_
    elif order == "SUMONE":
        RET = _TERM_ONE_ + _LOG_TERMS_
    elif order == "SUMTWO_NOLOG":
        RET = _TERM_ONE_ + _TERM_TWO_
    elif order == "SUMTWO":
        RET = _TERM_ONE_ + _TERM_TWO_ + _LOG_TERMS_
    elif order == "SUMTHREE":
        RET = _TERM_ONE_ + _TERM_TWO_ + _TERM_THREE_ + _LOG_TERMS_
    elif order == "SUMTHREE_NOLOG":
        RET = _TERM_ONE_ + _TERM_TWO_ + _TERM_THREE_
    else:
        "Only order = -1 implemented."

    return RET

def PhenomPL2PN(v, eta):
    """PhenomPL2PN(v, eta)
    v = orbital velocity : Cubic root of (Pi * Frequency (geometric))
    eta = symmetric mass-ratio

    Simple 2PN version of the orbital angular momentum L,
    without any spin terms expressed as a function of v.
    For IMRPhenomP(v2).

    Reference:
    - Boh&eacute; et al, 1212.5520v2 Eq 4.7 first line
    """
    return (eta/v) * (1. +  (3./2. + eta/6. )*v**2. + (3.373 - 19.*eta/8. - eta**2./24.)*v**4.)

def PhenomPBeta(omega, q, chi1x, chi1z):
    """PhenomPBeta(omega, q, chi1x, chi1z)
    (also called iota by PhenomP devs)
    NOTE the total mass is assumed to be one and the results
    are in units of the total mass (which is one)
    omega = orbital angular velocity
    q = mass-ratio
    chi1x = in-plane spin component on larger body
    chi1z = aligned spin component on larger body
    """
    mtot = 1.

    m1, m2 = m1_m2_M_q(mtot, q)
    # eta = symmetric mass-ratio (also called nu)
    eta = eta_from_q(q)
    # Sx (also called Sperp)
    Sx = m1**2. * chi1x
    # Sz (also called Spara)
    Sz = m1**2. * chi1z

    #compute PN orbital velocity
    v = omega**(1./3.)

    #s := Sp / (L + SL)
    s = Sx / (PhenomPL2PN(v, eta) + Sz)
    s2 = s*s

    beta = arccos( sign(s) * 1./sqrt(1+s2) )

    return beta

def PhenomPEpsilon(omega, q, chi1x, chi1z, order=-1):
    """PhenomPEpsilon(omega, q, chi1x, chi1z, order=-1)
    NOTE the total mass is assumed to be one and the results
    are in units of the total mass (which is one)
    omega = orbital angular frequency
    q = mass-ratio
    chi1x = in-plane spin component on larger body
    chi1z = aligned spin component on larger body
    order = PN order [default = -1 = all terms]
        other possible values:
        "ONE", "TWO", "THREE", "FOUR", "LOG"
        returns ONLY that partular term, NOT the summation.
    """
    mtot = 1.

    m1, m2 = m1_m2_M_q(mtot, q)
    # nu = symmetric mass-ratio (also called eta)
    nu = eta_from_q(q)
    dm = m1-m2

    #1/omega
    _TERM_ONE_ = (-35/192 + (5*dm)/(64*m1)) / omega

    #1/omega**(2/3)
    _TERM_TWO_ = ( \
            + (15*dm*m1*chi1z)/(128*mtot**2*nu) \
            - (35*m1**2*chi1z)/(128*mtot**2*nu) \
            ) / omega**(2/3)

    #1/omega**(1/3)
    _TERM_THREE_ = \
    +  (-5515/3072 + (4555*dm)/(7168*m1) \
    - (515*nu)/384 - (15*dm**2*nu)/(256*m1**2) \
    + (175*dm*nu)/(256*m1) \
    ) /omega**(1/3)

    #omega**(1/3)
    _TERM_FOUR_ = \
    ( \
    40121485/9289728 - (27895885*dm)/(21676032*m1) \
    + (39695*nu)/86016 + (1615*dm**2*nu)/(28672*m1**2) \
    + (265*dm*nu)/(14336*m1) + (955*nu**2)/576 \
    - (15*dm**3*nu**2)/(1024*m1**3) \
    + (35*dm**2*nu**2)/(256*m1**2) \
    - (2725*dm*nu**2)/(3072*m1) \
    + (15*dm*m1*pi*chi1z)/(16*mtot**2*nu) \
    - (35*m1**2*pi*chi1z)/(16*mtot**2*nu) \
    + (375*dm**2*m1**2*chi1z**2)/(256*mtot**4*nu) \
    - (1815*dm*m1**3*chi1z**2)/(256*mtot**4*nu) \
    + (1645*m1**4*chi1z**2)/(192*mtot**4*nu) \
    ) * omega**(1/3)

    #log(omega)
    _LOG_TERMS_ = \
    - (35*pi*log(omega))/48 + (5*dm*pi*log(omega))/(16*m1) \
    + (5*dm**2*chi1z*log(omega))/(16*mtot**2) \
    - (5*dm*m1*chi1z*log(omega))/(3*mtot**2) \
    + (2545*m1**2*chi1z*log(omega))/(1152*mtot**2) \
    - (2035*dm*m1*chi1z*log(omega))/(21504*mtot**2*nu) \
    + (2995*m1**2*chi1z*log(omega))/(9216*mtot**2*nu)

    if order == -1:
        RET = _TERM_ONE_ + _TERM_TWO_ + _TERM_THREE_ + _TERM_FOUR_ + _LOG_TERMS_
    elif order == "ONE":
        RET = _TERM_ONE_
    elif order == "TWO":
        RET = _TERM_TWO_
    elif order == "THREE":
        RET = _TERM_THREE_
    elif order == "FOUR":
        RET = _TERM_FOUR_
    elif order == "LOG":
        RET = _LOG_TERMS_
    elif order == "SUMONE_NOLOG":
        RET = _TERM_ONE_
    elif order == "SUMONE":
        RET = _TERM_ONE_ + _LOG_TERMS_
    elif order == "SUMTWO_NOLOG":
        RET = _TERM_ONE_ + _TERM_TWO_
    elif order == "SUMTWO":
        RET = _TERM_ONE_ + _TERM_TWO_ + _LOG_TERMS_
    elif order == "SUMTHREE":
        RET = _TERM_ONE_ + _TERM_TWO_ + _TERM_THREE_ + _LOG_TERMS_
    elif order == "SUMTHREE_NOLOG":
        RET = _TERM_ONE_ + _TERM_TWO_ + _TERM_THREE_
    else:
        "Only order = -1 implemented."

    return RET
