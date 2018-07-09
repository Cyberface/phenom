"""
functions that allow you to generate
data for the phenomp version 2 and version 3
precession angles in a common API
"""

import numpy as np
from phenom import HztoMf, PhenomPAlpha, PhenomPEpsilon, PhenomPBeta, chip_fun, chieffPH

try:
    import lal
except ImportError:
    raise ImportError('could not import lal')

try:
    import lalsimulation as lalsim
except ImportError:
    raise ImportError('could not import lalsimulation')

"""
helper functions to evaluate phenomPv2 and phenomPv3
precession angles
"""

def convert_from_cartesian_to_polar(x, y, z):
    """
    function to convert from 3d cartesian components to polar angles and vector magnitude.
    https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
    """
    mag = np.sqrt(x * x + y * y + z * z)
    if (np.abs(mag - 0.) < 1e-9):
        polar = 0.
        azimuthal = 0.
    else:
        polar = np.arccos(z / mag)
        azimuthal = np.arctan2(y, x)
    return mag, polar, azimuthal

def evaluate_phenomPv3_angles(flist, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, fref, ExpansionOrder=-1, PN="3"):
    """
    flist : REAL8Sequence - orbital frequency (Hz)
    # pv3 angles need orbital frequency (Hz)
    input: m1, m2 in SI units
    spins are dimensionless
    evaluated for frequencies flist (Hz)
    fref : gw reference frequency
    returns alpha, beta and epsilon angles

    ExpansionOrder=-1 (can be 1,2,3,4,5,-1) the expansion order to use in phiz and zeta
    PN=3 (can be "0","2" or "3") either 3PN or not
    """
    phiz_of_f = lal.CreateREAL8Sequence(len(flist.data))  # phiz or alpha
    zeta_of_f = lal.CreateREAL8Sequence(len(flist.data))  # zeta or epsilon
    costhetaL_of_f = lal.CreateREAL8Sequence(len(flist.data))  # costhetaL or beta

    costhetaL = 1.
    phiL = 0.

    chi1, theta1, phi1 = convert_from_cartesian_to_polar(s1x, s1y, s1z)
    costheta1 = np.cos(theta1)
    chi2, theta2, phi2 = convert_from_cartesian_to_polar(s2x, s2y, s2z)
    costheta2 = np.cos(theta2)

    if PN == "3":
        lalsim.ComputeAngles3PN(phiz_of_f, zeta_of_f, costhetaL_of_f, flist, m1, m2, costhetaL,
                                phiL, costheta1, phi1, chi1, costheta2, phi2, chi2, fref, ExpansionOrder)
    elif PN == "2":
        lalsim.ComputeAngles2PNNonSpinning(phiz_of_f, zeta_of_f, costhetaL_of_f, flist, m1, m2, costhetaL,
                             phiL, costheta1, phi1, chi1, costheta2, phi2, chi2, fref, ExpansionOrder)
    elif PN == "0":
        lalsim.ComputeAngles(phiz_of_f, zeta_of_f, costhetaL_of_f, flist, m1, m2, costhetaL,
                             phiL, costheta1, phi1, chi1, costheta2, phi2, chi2, fref, ExpansionOrder)
    else:
        raise(ValueError("Only PN = '0', '2' or '3' supported "))

    return phiz_of_f.data, np.arccos(costhetaL_of_f.data), zeta_of_f.data

def evaluate_phenomPv2_angles(flist, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, fref, order=-1):
    """
    flist = Real8Sequence of Orbital frequency
    input: m1, m2 in SI units
    spins are dimensionless
    fref : gw reference frequency
    order : PN order for alpha and epsilon
    returns alpha, beta and epsilon angles
    """

    q = m1 / m2

    if q < 1.:
        raise(ValueError("mass-ratio < 1. m1 must be the larger black hole"))

    # the angles are functions of dimensionless orbital angular GW frequency
    Momega = 2. * np.pi * HztoMf(flist.data, (m1 + m2) / lal.MSUN_SI)

    # divide by 2 to go from GW to orbital frequency
    Momega_ref = 2. * np.pi * HztoMf(fref/2, (m1 + m2) / lal.MSUN_SI)

    # compute chip and chieff
    chip, chi1l, chi2l = chip_fun(m1/ lal.MSUN_SI, m2/ lal.MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z)
    chieff = chieffPH(m1/ lal.MSUN_SI, m2/ lal.MSUN_SI, chi1l, chi2l)

    alpha_ref = PhenomPAlpha(Momega_ref, q, chip, chi1l, order=order)
    epsilon_ref = PhenomPEpsilon(Momega_ref, q, chip, chi1l, order=order)

    alpha = lal.CreateREAL8Sequence(len(flist.data))
    epsilon = lal.CreateREAL8Sequence(len(flist.data))
    beta = lal.CreateREAL8Sequence(len(flist.data))

    for i, f in enumerate(Momega):
        alpha.data[i] = PhenomPAlpha(f, q, chip, chi1l, order=order) - alpha_ref
        epsilon.data[i] = PhenomPEpsilon(f, q, chip, chi1l, order=order) - epsilon_ref
        #beta.data[i] = PhenomPBeta(f, q, chip, chieff)
        beta.data[i] = PhenomPBeta(f, q, chip, chi1l)

    return alpha.data, beta.data, epsilon.data

