import phenom

import matplotlib
# matplotlib.use('MacOSX')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from phenom.utils.utils import Constants, HztoMf, pad_to_pow_2

import lal
import lalsimulation as lalsim

from ifft import *


m1=52.384655
m2=16.384655
chi1x=0.5
chi1y=0.
chi1z=0.
chi2x=0.
chi2y=0.
chi2z=0.
# delta_f=1.0/128.
delta_f=1.0/256.
# delta_f=1.0/512.
# delta_f=1.0/1024.
f_min = 5.
fRef = 10.
inclination=2.85646439

Mtot = m1 + m2
q = m1/m2
eta = m1*m2/(m1+m2)**2.
srate=2**11.

ph_phpLAL = phenom.Waveform(approximant='IMRPhenomPv2_LAL',
                            m1=m1,
                            m2=m2,
                            chi1x=chi1x,
                            chi1y=chi1y,
                            chi1z=chi1z,
                            chi2x=chi2x,
                            chi2y=chi2y,
                            chi2z=chi2z,
                            delta_f=delta_f,
                            f_min=f_min,
                            fRef=fRef,
                            inclination=inclination)

htilde = ph_phpLAL.hptilde + 1.j * ph_phpLAL.hctilde

times, h = invfft( ph_phpLAL.flist_Hz, htilde, f0=4. )

plt.figure( figsize=( 14, 8 ) )
plt.plot(times, np.real(h))
plt.plot(times, np.abs(h))
# plt.xlim(-4, 1)
plt.xlim(-40, 1)
plt.savefig('./TD_amp.png', dpi=100)
plt.close()

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

def CallTDWaveform(approx, M, eta, chi1z, chi2z, chi1x, f_min=10, srate=2**14, f_ref=0.0, iota=0):
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
    S1y = 0.0
    S1z = chi1z
    S2x = 0.0
    S2y = 0.0
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
    # t = phenom.StoM(t, m1 + m2)


    maxindex = peakindex( hp.data.data )

    return t - t[maxindex], hp.data.data, hc.data.data

t={}
hp={}
hc={}

t['v3'], hp['v3'], hc['v3'] = CallTDWaveform("SEOBNRv3", Mtot, eta, chi1z, chi2z, chi1x, f_min=f_min, srate=srate, iota=inclination)
t['phenpv2'], hp['phenpv2'], hc['phenpv2'] = CallTDWaveform("IMRPhenomPv2", Mtot, eta, chi1z, chi2z, chi1x, f_min=f_min, srate=srate, iota=inclination)


plt.figure()
plt.plot(t['v3'], np.real(hp['v3']))
plt.plot(t['v3'], np.abs(hp['v3'] + 1.j * hc['v3']))
plt.xlim(-100, 1)
plt.savefig('./TD_amp_SEOBNRv3.png', dpi=100)
plt.close()



plt.figure( figsize=( 14, 8 ) )
plt.plot(times, np.real(h)*500, label='phenpv2-py')
plt.plot(times, np.abs(h)*500, label='phenpv2-py')
plt.plot(t['phenpv2'], np.real(hp['phenpv2']), label='phenpv2')
plt.plot(t['phenpv2'], np.abs(hp['phenpv2'] + 1.j * hc['phenpv2']), label='phenpv2-amp')
plt.xlim(-4, 1)
# plt.xlim(-40, 1)
plt.legend(loc='best')
plt.savefig('./TD_amp_compare_phenom.png', dpi=100)
plt.close()









plt.figure( figsize=( 14, 8 ) )
# plt.plot(times, np.real(h), label='phpv2')
# plt.plot(times, np.abs(h), label='phpv2-amp')
plt.plot(t['phenpv2'], np.real(hp['phenpv2']), label='phenpv2')
plt.plot(t['phenpv2'], np.abs(hp['phenpv2'] + 1.j * hc['phenpv2']), label='phenpv2-amp')
plt.plot(t['v3'], np.real(hp['v3']), label='EOB')
plt.plot(t['v3'], np.abs(hp['v3'] + 1.j * hc['v3']), label='EOB-amp')
plt.xlim(-10, 1)
plt.legend(loc='best')
plt.savefig('./compare_EOB_phenP_TD.png', dpi=100)
plt.close()



















#
