from numpy import sqrt

# Final Spin and Radiated Energy formulas described in 1508.07250

def FinalSpin0815_s(eta, s):
    """
    Formula to predict the final spin. Equation 3.6 arXiv:1508.07250
    s defined around Equation 3.6.
    """
    eta2 = eta*eta
    eta3 = eta2*eta
    eta4 = eta3*eta
    s2 = s*s
    s3 = s2*s
    s4 = s3*s

    return 3.4641016151377544*eta - 4.399247300629289*eta2 + \
    9.397292189321194*eta3 - 13.180949901606242*eta4 + \
    (1 - 0.0850917821418767*eta - 5.837029316602263*eta2)*s + \
    (0.1014665242971878*eta - 2.0967746996832157*eta2)*s2 + \
    (-1.3546806617824356*eta + 4.108962025369336*eta2)*s3 + \
    (-0.8676969352555539*eta + 2.064046835273906*eta2)*s4

def FinalSpin0815(eta, chi1, chi2):
    """
    Wrapper function for FinalSpin0815_s.
    """
    # Convention m1 >= m2
    Seta = sqrt(1.0 - 4.0*eta)
    m1 = 0.5 * (1.0 + Seta)
    m2 = 0.5 * (1.0 - Seta)
    m1s = m1*m1
    m2s = m2*m2
    # s defined around Equation 3.6 arXiv:1508.07250
    s = (m1s * chi1 + m2s * chi2)
    return FinalSpin0815_s(eta, s)

def EradRational0815_s(eta, s):
    """
    Formula to predict the total radiated energy. Equation 3.7 and 3.8 arXiv:1508.07250
    Input parameter s defined around Equation 3.7 and 3.8.
    """
    eta2 = eta*eta
    eta3 = eta2*eta
    eta4 = eta3*eta

    return ((0.055974469826360077*eta + 0.5809510763115132*eta2 - \
    0.9606726679372312*eta3 + 3.352411249771192*eta4)* \
    (1. + (-0.0030302335878845507 - 2.0066110851351073*eta \
    + 7.7050567802399215*eta2)*s))/(1. + (-0.6714403054720589 - \
    1.4756929437702908*eta + 7.304676214885011*eta2)*s)

def EradRational0815(eta, chi1, chi2):
    """
    Wrapper function for EradRational0815_s.
    """
    # Convention m1 >= m2
    Seta = sqrt(1.0 - 4.0*eta)
    m1 = 0.5 * (1.0 + Seta)
    m2 = 0.5 * (1.0 - Seta)
    m1s = m1*m1
    m2s = m2*m2
    # arXiv:1508.07250
    s = (m1s * chi1 + m2s * chi2) / (m1s + m2s)

    return EradRational0815_s(eta, s)

def fring(eta, chi1, chi2, finspin):
    """
    fring is the real part of the ringdown frequency
    1508.07250 figure 9
    """
    from scipy import interpolate
    from phenom.utils.QNMdata.QNMphenomd import QNMData

    if (finspin > 1.0):
        print("PhenomD fring function: final spin > 1.0 not supported\n")

    ifring = interpolate.interp1d(QNMData['QNMData_a'], QNMData['QNMData_fring'])

    return ifring(finspin) / (1.0 - EradRational0815(eta, chi1, chi2))

def fdamp(eta, chi1, chi2, finspin):
    """
    fdamp is the complex part of the ringdown frequency
    1508.07250 figure 9
    """
    from scipy import interpolate
    from phenom.utils.QNMdata.QNMphenomd import QNMData

    if (finspin > 1.0):
        print("PhenomD fdamp function: final spin > 1.0 not supported\n")

    ifdamp = interpolate.interp1d(QNMData['QNMData_a'], QNMData['QNMData_fdamp'])

    return ifdamp(finspin) / (1.0 - EradRational0815(eta, chi1, chi2))


def FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH(m1, m2, chip, chi1z, chi2z):
    #TODO: Need to write a function that takes in chi1x,...chi2z, etc and returns chip etc.
    """
    Wrapper for final-spin formula based on:
    - IMRPhenomD's FinalSpin0815() for aligned spins.

    We use their convention m1>m2
    and put <b>all in-plane spin on the larger BH</b>.

    In the aligned limit return the FinalSpin0815 value.

    m1     #/**< Mass of companion 1 (solar masses) */
    m2     #/**< Mass of companion 2 (solar masses) */
    chip   #/**< chip parameter, Eq. 3.4 PhysRevD.91.024043 */
    chi1z  #/**< Aligned spin of BH 1 */
    chi2z  #/**< Aligned spin of BH 2 */
    """
    from math import copysign

    m1 = float(m1)
    m2 = float(m2)
    chip = float(chip)
    chi1z = float(chi1z)
    chi2z = float(chi2z)

    M = m1 + m2
    eta = m1*m2/(M*M)

    if (m1 >= m2):
        q_factor = m1/M
        af_parallel = FinalSpin0815(eta, chi1z, chi2z)
    else:
        q_factor = m2/M
        af_parallel = FinalSpin0815(eta, chi2z, chi1z)

    Sperp = chip * q_factor*q_factor

    af = copysign(1.0, af_parallel) * sqrt(Sperp*Sperp + af_parallel*af_parallel)
    return af
