
# source /Users/sebastian/work/git/phenomp-upgrade/lscsoft/activate-lal.sh
import sys
import os

# sys.path.append(
#     '/Users/sebastian/work/git/phenomp-upgrade/lscsoft/lalsimtools-analysis/time-domain-plot/')

# import angles

from phenom.testing import prec_angle_helper as pah

import lal
import lalsimulation as lalsim

os.environ["LAL_WAVEFORM_TOOLS_PATH"] = "/Users/sebastian/work/git/lal_waveform_tools/python/"
import lalsimtools
from lalsimtools import skymaxmatch as smm

from phenom.utils import f_SchwarzISCO, M_eta_m1_m2
from phenom.utils.remnant import FinalSpin0815, fring
from phenom import MftoHz, HztoMf
import phenom

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import matplotlib.gridspec as gridspec
import numpy as np

from pycbc.pnutils import hybrid_meco_frequency


def get_prec_angles(pars, fring_Hz, phenomp_version="pv2", ExpOrder=5, PNOrbL="3", pv2Order=-1):


    my_f_SchISCO = f_SchwarzISCO((pars['m1'] + pars['m2']) / lal.MSUN_SI)

    # f_orb_max = int(1.2 * fring_Hz / 2)
    f_orb_max = int(my_f_SchISCO / 2.)
    f_gw_min = pars['f_min']
    f_gw_ref = pars['f_ref']
    # f_gw_ref = f_gw_min * 2
    f_orb_min = f_gw_min/2
    f_orb_df = 0.01
    Npoints = int((f_orb_max - f_orb_min) / f_orb_df)

    fGW_hz = lal.CreateREAL8Sequence(Npoints)
    fGW_hz.data = np.arange(f_orb_min, f_orb_max, f_orb_df)
    # convert from orbital to GW
    fGW_hz.data *= 2

    f_orb_hz = lal.CreateREAL8Sequence(Npoints)
    f_orb_hz.data = np.arange(f_orb_min, f_orb_max, f_orb_df)

    # pv2 and pv3 angles need orbital frequency (Hz)

    if phenomp_version == "pv2":

        alpha, betav2, epsilon = pah.evaluate_phenomPv2_angles(
            f_orb_hz, pars['m1'], pars['m2'],
            pars['S1x'], pars['S1y'], pars['S1z'],
            pars['S2x'], pars['S2y'], pars['S2z'],
            f_gw_ref, pv2Order)

        alpha = alpha - alpha[0]
        beta = betav2
        epsilon = epsilon - epsilon[0]

    elif phenomp_version == "pv3":

        phiz, beta, zeta = pah.evaluate_phenomPv3_angles(
            f_orb_hz, pars['m1'], pars['m2'],
            pars['S1x'], pars['S1y'], pars['S1z'],
            pars['S2x'], pars['S2y'], pars['S2z'], f_gw_ref,
            ExpOrder, PNOrbL)

        alpha = phiz - phiz[0]
        beta = beta
        epsilon = zeta - zeta[0]

    prec_angles = {}

    # return prec_angles as a function of GW frequ

    # The plots are a function of the GW frequency
    prec_angles['Mf'] = HztoMf(fGW_hz.data, (pars['m1']+pars['m2'])/lal.MSUN_SI)
    prec_angles['f'] = fGW_hz.data
    prec_angles['f_orb_hz'] = f_orb_hz.data
    prec_angles['alpha'] = alpha
    prec_angles['beta'] = beta
    prec_angles['epsilon'] = epsilon

    return prec_angles


def mainfunc(case):
    if case == 1:
        mass1, mass2 = phenom.m1_m2_M_q(10, 1)
        spin1x, spin1y, spin1z = 0.9, 0., 0.
        spin2x, spin2y, spin2z = 0., 0., 0.
    else:
        print("unknown case - exiting")
        import sys; sys.exit()

    SIGNAL_DISTANCE = 1e3 * 1e6 * lal.PC_SI
    SIGNAL_INCLINATION = np.pi / 2.
    SIGNAL_PHASE = np.pi
    longAscNodes = 0.
    sample_rate = 2**10.
    # sample_rate = 2**12.
    # sample_rate = 2**13.
    deltaT = 1. / sample_rate
    params = None
    # f_lower = 8.
    # f_ref = 20.
    # f_lower = 10.
    # f_ref = 10.

    f_lower = 3.
    f_ref = 3.


    # f_lower = 20.
    # f_ref = 20.

    mtot, eta = M_eta_m1_m2(mass1, mass2)

    f_SchISCO = f_SchwarzISCO(mtot)
    f_hybrid_meco = hybrid_meco_frequency(mass1, mass2, spin1z, spin2z)

    fin_spin = FinalSpin0815(eta, spin1z, spin2z)
    fring_Mf = fring(eta, spin1z, spin2z, fin_spin)
    fring_Hz = MftoHz(fring_Mf, mtot)

    # generator_function="get_td_waveform"
    generator_function = "SimInspiralTD"

    pv2 = "IMRPhenomPv2"
    pv3 = "IMRPhenomPv3"

    signal_pars = {
        "m1": mass1 * lal.MSUN_SI,
        "m2": mass2 * lal.MSUN_SI,
        "S1x": spin1x,
        "S1y": spin1y,
        "S1z": spin1z,
        "S2x": spin2x,
        "S2y": spin2y,
        "S2z": spin2z,
        "distance": SIGNAL_DISTANCE,
        "inclination": SIGNAL_INCLINATION,
        "phiRef": SIGNAL_PHASE,
        "longAscNodes": longAscNodes,
        "eccentricity": 0.,
        "meanPerAno": 0.,
        "deltaT": deltaT,
        "f_min": f_lower,
        "f_ref": f_ref,
        "LALparams": params,
        "approximant": ""
    }

    pv2_pars = signal_pars.copy()
    pv2_pars.update({
        "approximant": lalsim.GetApproximantFromString(pv2),
    })

    pv3_pars = signal_pars.copy()
    pv3_pars.update({
        "approximant": lalsim.GetApproximantFromString(pv3),
        "f_min":0.6 * f_lower
    })

    #
    # get pv2 and pv3 precession angles
    pv2angles = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv2",
                                ExpOrder=5, PNOrbL="3")

    pv2angles_leadingorder = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv2",
                                ExpOrder=5, PNOrbL="3", pv2Order="ONE")

    # pv3angles_Exp5_PN3 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                      ExpOrder=5, PNOrbL="3")
    # pv3angles_Exp5_PN2 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                      ExpOrder=5, PNOrbL="2")
    # pv3angles_Exp5_PN0 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                      ExpOrder=5, PNOrbL="0")

    pv3angles_Exp5_PN2 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
                                        ExpOrder=5, PNOrbL="2")
    # pv3angles_Exp4_PN2 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                      ExpOrder=4, PNOrbL="2")
    # pv3angles_Exp3_PN2 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                      ExpOrder=3, PNOrbL="2")
    # pv3angles_Exp2_PN2 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                      ExpOrder=2, PNOrbL="2")
    # pv3angles_Exp1_PN2 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                      ExpOrder=1, PNOrbL="2")
    pv3angles_ExpAll_PN2 = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
                                        ExpOrder=-1, PNOrbL="2")

    # get beta angle also at 0PN and 3PN
    # tmp = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                ExpOrder=-1, PNOrbL="0")
    # pv3_beta_PN0 = {}
    # pv3_beta_PN0['f'] = tmp['f']
    # pv3_beta_PN0['beta'] = tmp['beta']
    # tmp = get_prec_angles(pv2_pars, fring_Hz, phenomp_version="pv3",
    #                                ExpOrder=-1, PNOrbL="3")
    # pv3_beta_PN3 = {}
    # pv3_beta_PN3['f'] = tmp['f']
    # pv3_beta_PN3['beta'] = tmp['beta']



    plt.figure()
    # plt.plot(pv2angles['Mf'], pv2angles['alpha'], label="PhenomPv2")
    plt.plot(pv2angles['f_orb_hz'], pv2angles['alpha'], label="PhenomPv2")
    plt.plot(pv3angles_Exp5_PN2['f_orb_hz'], pv3angles_Exp5_PN2['alpha'], label="PhenomPv3")

    # plt.ylim(-10, 5)
    # plt.ylim(-10,110)
    # plt.ylim(0, 12)
    # plt.ylim(0, 60)
    # plt.xlim(0,0.02)
    plt.xlabel("orbital freq (Hz)")
    plt.ylabel(r"$\alpha(f)-\alpha(f_0)$")
    plt.legend()
    plt.tight_layout()
    plt.savefig('alpha.png')


if __name__ == "__main__":


    mainfunc(1)