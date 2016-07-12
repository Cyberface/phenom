from __future__ import division

from phenom.utils.utils import M_eta_m1_m2, chipn, Constants, chieffPH, WignerdCoefficients
from phenom.utils import swsh
from phenom.waveform import PhenomD
from phenom.pn.pn import PhenomPAlpha, PhenomPBeta, PhenomPEpsilon

from numpy import exp, arange, zeros, array, conj

class PhenomP():


    #TODO:
    #1. spherical harmonics

    # highest frequency the model will go to. Inherited from PhenomD at the moment
    # Dimensionless frequency (Mf) at which define the end of the waveform
    Mf_CUT = 0.2


    """docstring for PhenomP"""
    def __init__(self,  m1=10., m2=10., chi1x=0.,  chi1y=0., chi1z=0., chi2x=0., chi2y=0., chi2z=0., f_min=20., f_max=0., delta_f=1/64., distance=1e6 * Constants.PC_SI, fRef=0., phiRef=0.):
        """
        input:
        m1 (Msun)
        m2 (Msun)
        chi1z (dimensionless)
        chi2z (dimensionless)
        f_min (Hz)
        f_max (Hz)
        delta_f (sample rate Hz)
        distance (m) : Default 1e6 * Constants.PC_SI = 1 mega parsec
        fRef (reference frequency Hz)
        phiRef (orbital phase at fRef)"""

        # NOTE PhenomP in LAL assumes that m2>m1. We use the opposite convention here!

        # enforce m1 >= m2 and chi1 is on m1
        if m1<m2: # swap spins and masses
            chi1z, chi2z = chi2z, chi1z
            m1, m2 = m2, m1

        self.p                = {}
        self.p['m1']          = float(m1)
        self.p['m2']          = float(m2)
        self.p['chi1x']       = float(chi1x)
        self.p['chi1y']       = float(chi1y)
        self.p['chi1z']       = float(chi1z)
        self.p['chi2x']       = float(chi2x)
        self.p['chi2y']       = float(chi2y)
        self.p['chi2z']       = float(chi2z)
        self.p['f_min']       = float(f_min)
        self.p['f_max']       = float(f_max)
        self.p['delta_f']     = float(delta_f)
        self.p['distance']    = float(distance)
        self.p['fRef']        = float(fRef)
        self.p['phiRef']      = float(phiRef)

        self.p['Mtot'], self.p['eta'] = M_eta_m1_m2(self.p['m1'], self.p['m2'])


        # Dimensionfull aligned spin.
        self.p['SL'] = self.p['chi1z']*self.p['m1']**2. + self.p['chi2z']*self.p['m2']**2.
        #Dimensionfull spin component in the orbital plane. S_perp = S_2_perp
        # self.p['Sperp'] = self.p['chip']*(self.p['m1']**2.)
        #NOTE: Currently code only works if m1>m2 and chip=chi1x
        self.p['Sperp'] = self.p['chi1x']*(self.p['m1']**2.)
        #NOTE: In LAL PhenomP uses q=m2/m1, q>=1, we use opposite.
        self.p['q'] = self.p['m1'] / self.p['m2']

        self.p['M_sec'] = self.p['Mtot'] * Constants.MTSUN_SI # Conversion factor Hz -> dimensionless frequency

        self.piM = Constants.LAL_PI * self.p['M_sec']

        #chipn is used in the aligned model 'phenomD'
        self.p['chipn']  = chipn(self.p['eta'], self.p['chi1z'], self.p['chi2z'])
        #chieff is used in the twisting part
        self.p['chieff'] = chieffPH(self.p['m1'], self.p['m2'], self.p['chi1z'], self.p['chi2z'])

        if self.p['f_max'] == 0. : self.p['f_max'] = self.Mf_CUT / self.p['M_sec'] # converted from Mf to Hz
        if self.p['fRef'] == 0.  : self.p['fRef'] = self.p['f_min']
        self.p['MfRef'] = self.p['fRef'] * self.p['M_sec']
        # orbital angular frequency = pi * f_GW = omega_GW / 2
        self.p['omega_Ref'] = self.p['MfRef'] * Constants.LAL_PI

        #Compute Ylm's only once and pass them to PhenomPCoreOneFrequency() below.
        Y2m = {}
        thetaJ = 0.2345 #TODO: FIX ME THIS IS JUST A TMP!
        ytheta  = thetaJ
        yphi    = 0.
        Y2m['Y2m2'] = swsh.SWSH(2, -2, ytheta, yphi).val
        Y2m['Y2m1'] = swsh.SWSH(2, -1, ytheta, yphi).val
        Y2m['Y20']  = swsh.SWSH(2,  0, ytheta, yphi).val
        Y2m['Y21']  = swsh.SWSH(2,  1, ytheta, yphi).val
        Y2m['Y22']  = swsh.SWSH(2,  2, ytheta, yphi).val

        #compute amplitude and phase from aligned spin model
        #the amplitude is unscaled.
        aligned_amp, aligned_phase = self._generate_aligned_spin_approx(self.p)

        print "testing"
        print aligned_amp[0], aligned_phase[0]

        #'amp_scale = amp0 in LALSimIMRPhenomP.c'
        self.p['amp_scale'] = self.p['Mtot'] * Constants.MRSUN_SI * self.p['Mtot'] * Constants.MTSUN_SI / self.p['distance']

        # aligned spin strain
        self.hP = self.p['amp_scale'] * aligned_amp * exp(-1.j * aligned_phase)

        alpha_at_omega_Ref = self._alpha_precessing_angle(self.p['omega_Ref'], self.p)
        epsilon_at_omega_Ref = self._epsilon_precessing_angle(self.p['omega_Ref'], self.p)
        print "testing"
        print alpha_at_omega_Ref
        print epsilon_at_omega_Ref

        #TODO: Code up the spherical harmonics into a class

        omega_flist_hz = arange(self.p['f_min'], self.p['f_max'], self.p['delta_f']) * Constants.LAL_PI
        self.alpha = zeros(len(omega_flist_hz))
        self.epsilon = zeros(len(omega_flist_hz))
        for i in range(len(omega_flist_hz)):
            omega = omega_flist_hz[i] * self.p['M_sec']
            self.alpha[i] = self._alpha_precessing_angle(omega, self.p)
            self.epsilon[i] = self._epsilon_precessing_angle(omega, self.p)

        self.alpha -= alpha_at_omega_Ref #TODO When I have the 'compute phenomP model parameters' function will need to subtract off alpha0 too!!
        self.epsilon -= epsilon_at_omega_Ref


        cBetah, sBetah = WignerdCoefficients(omega**(1./3.), self.p['SL'], self.p['eta'], self.p['Sperp'])

        cBetah2 = cBetah*cBetah
        cBetah3 = cBetah2*cBetah
        cBetah4 = cBetah3*cBetah
        sBetah2 = sBetah*sBetah
        sBetah3 = sBetah2*sBetah
        sBetah4 = sBetah3*sBetah

        """Compute Wigner d coefficients
        The expressions below agree with refX [Goldstein?] and Mathematica
        d2  = Table[WignerD[{2, mp, 2}, 0, -\[Beta], 0], {mp, -2, 2}]
        dm2 = Table[WignerD[{2, mp, -2}, 0, -\[Beta], 0], {mp, -2, 2}]
        """
        sqrt_6 = 2.44948974278317788
        d2   = array([sBetah4, 2*cBetah*sBetah3, sqrt_6*sBetah2*cBetah2, 2*cBetah3*sBetah, cBetah4])
        #Exploit symmetry d^2_{-2,-m} = (-1)^m d^2_{2,m}
        dm2  = array([d2[4], -d2[3], d2[2], -d2[1], d2[0]])

        Y2mA = array([Y2m['Y2m2'], Y2m['Y2m1'], Y2m['Y20'], Y2m['Y21'], Y2m['Y22']])
        hp_sum = 0.
        hc_sum = 0.
        #Sum up contributions to \tilde h+ and \tilde hx
        #Precompute powers of e^{i m alpha}
        cexp_i_alpha = exp(+1.j*self.alpha)
        cexp_2i_alpha = cexp_i_alpha*cexp_i_alpha
        cexp_mi_alpha = 1.0/cexp_i_alpha
        cexp_m2i_alpha = cexp_mi_alpha*cexp_mi_alpha
        cexp_im_alpha = array([cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha])

        for m in [-2, -1, 0, 1, 2]:
            T2m   = cexp_im_alpha[-m+2] * dm2[m+2] *      Y2mA[m+2]  # = cexp(-I*m*alpha) * dm2[m+2] *      Y2mA[m+2]
            Tm2m  = cexp_im_alpha[m+2]  * d2[m+2]  * conj(Y2mA[m+2]) # = cexp(+I*m*alpha) * d2[m+2]  * conj(Y2mA[m+2])
            hp_sum +=     T2m + Tm2m
            hc_sum += +1.j*(T2m - Tm2m)

        eps_phase_hP = exp(-2*1.j*self.epsilon) * self.hP / 2.0
        self.hp = eps_phase_hP * hp_sum
        self.hc = eps_phase_hP * hc_sum

    def _generate_aligned_spin_approx(self, p):
        """
        returns amplitude and phase form aligned spin model.
        Currently only implemented PhenomD
        """

        #QUESTION: Should I set fRef and phiRef to zero here when
        #generating phenomD for phenomP to use?
        #I feel like I should....

        #generate phenomD
        ph = PhenomD(m1=p['m1'], m2=p['m2'],
                    chi1z=p['chi1z'], chi2z=p['chi2z'],
                    f_min=p['f_min'], f_max=p['f_max'],
                    delta_f=p['delta_f'],
                    distance=p['distance'],
                    # fRef=0., phiRef=0.)
                    fRef=p['fRef'], phiRef=p['phiRef'],
                    finspin_func="FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH",
                    chi1x=p['chi1x'])
        print "finspin = {0}".format(ph.model_pars['finspin'])
        #could change this to generate at a single frequency point.
        ph.IMRPhenomDGenerateFD()
        ph.getampandphase(ph.htilde)
        #this means we now have ph.amp and ph.phase.
        #this returns the physically scaled amplitude to the input distance
        #we don't want that for phenomP
        #so we undo the scaling
        ph.amp /= ph.model_pars['amp0']
        return ph.amp, ph.phase

    def _alpha_precessing_angle(self, omega, p):
        """This function uses omega (the orbital angular frequency) as its argument!
        orbital angular frequency = pi * f_GW = omega_GW / 2
        """
        q = p['q']
        chi1x = p['chi1x']
        chi1z = p['chi1z']
        return PhenomPAlpha(omega, q, chi1x, chi1z, order=-1)

    def _epsilon_precessing_angle(self, omega, p):
        """This function uses omega (the orbital angular frequency) as its argument!
        orbital angular frequency = pi * f_GW = omega_GW / 2
        """
        q = p['q']
        chi1x = p['chi1x']
        chi1z = p['chi1z']
        return PhenomPEpsilon(omega, q, chi1x, chi1z, order=-1)








#
