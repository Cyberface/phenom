import phenom
import numpy as np

import lal
import lalsimulation as lalsim

def peakindex(x):
    return list(np.abs(x)).index(np.max(np.abs(x)))

def q_from_eta(eta):
    """
    Assumes m1 >= m2
    converts symmetric-mass-ratio to mass-ratio
    input: eta
    output: q
    """
    Seta = np.sqrt(1. - 4. * eta)
    return (1. + Seta - 2. * eta)/(2. * eta)

def m1_m2_M_eta(M, eta):
    """
    Assumes m1 >= m2
    Computes the component masses m1 and m2
    from the total mass and symmetric mass-ratio.
    input: M, eta
    output: m1, m2
    """
    Seta = np.sqrt(1. - 4. * eta)
    m1 = 1./2. * (M + Seta * M)
    m2 = 1./2. * (M - Seta * M)
    return m1, m2

def CallTDWaveform(approx="IMRPhenomD", M=100., eta=0.25,
    chi1x=0., chi1y=0., chi1z=0.,
    chi2x=0., chi2y=0., chi2z=0.,
    f_min=10, srate=2**14, f_ref=0.0, iota=0):
    """assuming m1>=m2"""
    deltaT=1./srate
    q = q_from_eta(eta)
    m1, m2 = m1_m2_M_eta(M, eta)
    m1_SI = m1 * lal.MSUN_SI
    m2_SI = m2 * lal.MSUN_SI
    # print 'chi_eff = ', (m1*chi1 + m2*chi2)/M
    # f_max_Hz = f_max / (M * lal.MTSUN_SI)
    phiRef = 0.0
    S1x = chi1x
    S1y = chi1y
    S1z = chi1z
    S2x = chi2x
    S2y = chi2y
    S2z = chi2z
    r = 1e6 * lal.PC_SI
    z = 0.0
    i = iota
    lambda1 = 0.0
    lambda2 = 0.0
    waveFlags = None
    nonGRparams = None
    amplitudeO = -1
    phaseO = -1
    # approximant = lalsim.GetApproximantFromString("IMRPhenomPv2")
    approximant = lalsim.GetApproximantFromString(approx)
    # print approximant
    hp, hc = lalsim.SimInspiralChooseTDWaveform(phiRef,
                                            deltaT,
                                            m1_SI, m2_SI,
                                            S1x, S1y, S1z, S2x, S2y, S2z,
                                            f_min, f_ref,
                                            r,
                                            i,
                                            lambda1, lambda2, waveFlags, nonGRparams,
                                            amplitudeO, phaseO,
                                            approximant)
    t = np.arange(hp.data.length) * hp.deltaT

    #convert to units of total mass (dimensionless)
#     t = phenom.StoM(t, m1 + m2) # keeping time in seconds

    return t, hp.data.data, hc.data.data

def CallFDWaveform(approx="IMRPhenomD", M=100., eta=0.25,
    chi1x=0., chi1y=0., chi1z=0.,
    chi2x=0., chi2y=0., chi2z=0.,
    f_min=10, f_max=0, srate=2**14, f_ref=0.0, iota=0):
    """assuming m1>=m2"""
    deltaF=1./srate
    q = q_from_eta(eta)
    m1, m2 = m1_m2_M_eta(M, eta)
    m1_SI = m1 * lal.MSUN_SI
    m2_SI = m2 * lal.MSUN_SI
    # print 'chi_eff = ', (m1*chi1 + m2*chi2)/M
    # f_max_Hz = f_max / (M * lal.MTSUN_SI)
    phiRef = 0.0
    S1x = chi1x
    S1y = chi1y
    S1z = chi1z
    S2x = chi2x
    S2y = chi2y
    S2z = chi2z
    r = 1e6 * lal.PC_SI
    z = 0.0
    i = iota
    lambda1 = 0.0
    lambda2 = 0.0
    waveFlags = None
    nonGRparams = None
    amplitudeO = -1
    phaseO = -1
    # approximant = lalsim.GetApproximantFromString("IMRPhenomPv2")
    approximant = lalsim.GetApproximantFromString(approx)
    hp, hc = lalsim.SimInspiralChooseFDWaveform(phiRef,
                                            deltaF,
                                            m1_SI, m2_SI,
                                            S1x, S1y, S1z, S2x, S2y, S2z,
                                            f_min, f_max, f_ref,
                                            r,
                                            i,
                                            lambda1, lambda2, waveFlags, nonGRparams,
                                            amplitudeO, phaseO,
                                            approximant)
    f = np.arange(hp.data.length) * hp.deltaF

    #convert to units of total mass (dimensionless)
#     f = phenom.HztoMf(f, m1 + m2) # keeping frequency in Hz

    return f, hp.data.data, hc.data.data
