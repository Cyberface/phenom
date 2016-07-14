from __future__ import division

from phenom.utils.utils import M_eta_m1_m2, chipn, Constants, chieffPH, WignerdCoefficients
from phenom.utils import swsh
from phenom.waveform import PhenomD
from phenom.pn.pn import PhenomPAlpha, PhenomPBeta, PhenomPEpsilon, PhenomPL2PN

from numpy import exp, arange, zeros, array, conj, sin, cos, dot, max, arccos, arctan2
from numpy.linalg import norm

class PhenomP():


    #TODO:
    #1. spherical harmonics

    # highest frequency the model will go to. Inherited from PhenomD at the moment
    # Dimensionless frequency (Mf) at which define the end of the waveform
    Mf_CUT = 0.2


    """docstring for PhenomP"""
    def __init__(self,  m1=10., m2=10., chi1x=0.,  chi1y=0., chi1z=0., chi2x=0., chi2y=0., chi2z=0., f_min=20., f_max=0., delta_f=1/64., distance=1e6 * Constants.PC_SI, fRef=0., phiRef=0., inclination=0.):
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
        phiRef (orbital phase at fRef)
        inclination (rad) angle between orbital anglar momentm and z-axis ?"""

        # NOTE PhenomP in LAL assumes that m2>m1. We use the opposite convention here!

        # enforce m1 >= m2 and chi1 is on m1
        if m1<m2: # swap spins and masses
            # chi1z, chi2z = chi2z, chi1z
            chi1x, chi1y, chi1z, chi2x, chi2y, chi2z = chi2x, chi2y, chi2z, chi1x, chi1y, chi1z
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
        self.p['inclination'] = float(inclination)

        self.p['Mtot'], self.p['eta'] = M_eta_m1_m2(self.p['m1'], self.p['m2'])

        #NOTE: In LAL PhenomP uses q=m2/m1, q>=1, we use opposite.
        self.p['q'] = self.p['m1'] / self.p['m2']

        self.p['M_sec'] = self.p['Mtot'] * Constants.MTSUN_SI # Conversion factor Hz -> dimensionless frequency

        if self.p['f_max'] == 0. : self.p['f_max'] = self.Mf_CUT / self.p['M_sec'] # converted from Mf to Hz
        if self.p['fRef'] == 0.  : self.p['fRef'] = self.p['f_min']
        self.p['MfRef'] = self.p['fRef'] * self.p['M_sec']
        # orbital angular frequency = pi * f_GW = omega_GW / 2
        self.p['omega_Ref'] = self.p['MfRef'] * Constants.LAL_PI

        #TODO: somehow need to add inclination angle as an input
        #But not exactly sure what the inclination angle is...

        self.p.update(self._PhenomPCalculateModelParameters(self.p))

        #Sperp and SL are used in the computation of the Wigner-d matrix
        #Dimensionless spin component in the orbital plane. S_perp = S_2_perp
        #NOTE: This assumes the m1 + m2 = 1. and hence then /Mtot**2
        self.p['Sperp'] = self.p['chip']*(self.p['m1'] / self.p['Mtot'])**2.
        #Dimensionfull aligned spin
        #This also assumes that m1 + m2 = 1. and hence then /Mtot**2
        self.p['SL'] = (self.p['chi1_l']*self.p['m1']**2. + self.p['chi2_l']*self.p['m2']**2.) / self.p['Mtot']**2.

        print "self.p['SL'] = ", self.p['SL']
        print "self.p['Sperp'] = ", self.p['Sperp']

        # import sys
        # sys.exit(1)

        #chipn is used in the aligned model 'phenomD'
        # self.p['chipn']  = chipn(self.p['eta'], self.p['chi1z'], self.p['chi2z'])
        self.p['chipn']  = chipn(self.p['eta'], self.p['chi1_l'], self.p['chi2_l'])
        #chieff is used in the twisting part
        # self.p['chieff'] = chieffPH(self.p['m1'], self.p['m2'], self.p['chi1z'], self.p['chi2z'])
        self.p['chieff'] = chieffPH(self.p['m1'], self.p['m2'], self.p['chi1_l'], self.p['chi2_l'])

        # dimensionless aligned spin of the largest BH
        self.p['chil'] = (self.p['Mtot'] / self.p['m1']) * self.p['chieff']

        print "self.p['chipn'] = ", self.p['chipn']
        print "self.p['chieff'] = ", self.p['chieff']
        print "self.p['chil'] = ", self.p['chil']

        #Compute Ylm's only once and pass them to PhenomPCoreOneFrequency() below.
        Y2m = {}
        ytheta  = self.p['thetaJ']
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
        print "self.p['omega_Ref'] = ", self.p['omega_Ref']
        print "alpha_at_omega_Ref = ", alpha_at_omega_Ref
        print "epsilon_at_omega_Ref = ", epsilon_at_omega_Ref

        #TODO: Code up the spherical harmonics into a class

        omega_flist_hz = arange(self.p['f_min'], self.p['f_max'], self.p['delta_f']) * Constants.LAL_PI
        self.alpha = zeros(len(omega_flist_hz))
        self.epsilon = zeros(len(omega_flist_hz))
        for i in range(len(omega_flist_hz)):
            omega = omega_flist_hz[i] * self.p['M_sec']
            self.alpha[i] = self._alpha_precessing_angle(omega, self.p)
            self.epsilon[i] = self._epsilon_precessing_angle(omega, self.p)

        self.alpha -= alpha_at_omega_Ref - self.p['alpha0']
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
                    chi1z=p['chi1_l'], chi2z=p['chi2_l'],
                    f_min=p['f_min'], f_max=p['f_max'],
                    delta_f=p['delta_f'],
                    distance=p['distance'],
                    # fRef=0., phiRef=0.)
                    fRef=p['fRef'], phiRef=p['phiRef'],
                    finspin_func="FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH",
                    chip=p['chip'])
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
        chi1x = p['chip'] #TODO This needs to be chip
        chi1z = p['chil'] #TODO dimensionless aligned spin of the largest BH
        return PhenomPAlpha(omega, q, chi1x, chi1z, order=-1)

    def _epsilon_precessing_angle(self, omega, p):
        """This function uses omega (the orbital angular frequency) as its argument!
        orbital angular frequency = pi * f_GW = omega_GW / 2
        """
        #TODO Change doc string in pn.py to be correct values.
        #TODO Change it in v3utils too!
        q = p['q']
        chi1x = p['chip'] #TODO This needs to be chip
        chi1z = p['chil'] #TODO dimensionless aligned spin of the largest BH
        return PhenomPEpsilon(omega, q, chi1x, chi1z, order=-1)


    def _PhenomPCalculateModelParameters(self, p):
        """
        'SimIMRPhenomPCalculateModelParameters' in 'LALSimIMRPhenomP.c'
        convention m1 >= m2, q = m1/m2 >= 1
        m1,            Mass of companion 1 (Msun)
        m2,            Mass of companion 2 (Msun)
        f_ref,         Reference GW frequency (Hz)
        inc,           inclination angle (rad), angle between L and z (z \equiv N , the line of sight in phenomP)
        s1x,           Initial value of s1x: dimensionless spin of BH 1
        s1y,           Initial value of s1y: dimensionless spin of BH 1
        s1z,           Initial value of s1z: dimensionless spin of BH 1
        s2x,           Initial value of s2x: dimensionless spin of BH 2
        s2y,           Initial value of s2y: dimensionless spin of BH 2
        s2z,           Initial value of s2z: dimensionless spin of BH 2
        """

        print "p['m1'] = ", p['m1']
        print "p['m2'] = ", p['m2']
        if p['m1'] < p['m2']:
            raise ValueError('m1 = {0}, m2 = {1}. Convention error, this function needs m1 > m2'.format(p['m1'], p['m2']))

        #check that the spin magnitude is <=1
        if norm([p['chi1x'], p['chi1y'], p['chi1z']]) > 1.:
            raise ValueError('chi1 has a magnitude > 1')
        if norm([p['chi2x'], p['chi2y'], p['chi2z']]) > 1.:
            raise ValueError('chi2 has a magnitude > 1')

        m1_2 = p['m1']**2.
        m2_2 = p['m2']**2.

        #we start out in the Lhat = zhat frame
        #and define the spin w.r.t this frame.
        #Then, we incline the orbital frame w.r.t to the z-axis
        #by the angle inc.
        #This is done by a rotation about the y-axis, so the y-components do not change
        #in LAL this step is done in XLALSimInspiralInitialConditionsPrecessingApproxs in LALSimInspiralSpinTaylor.c
        #But it's simple so I just do it in this function.

        print("spins before rotation by {0} = ".format(p['inclination']))
        print("chi1x = {0}, chi1y = {1}, chi1z = {2}".format(p['chi1x'], p['chi1y'], p['chi1z']))
        print("chi2x = {0}, chi2y = {1}, chi2z = {2}".format(p['chi2x'], p['chi2y'], p['chi2z']))


        p['chi1x'], p['chi1z'] = self.ROTATEY(p['inclination'], p['chi1x'], p['chi1z'])
        p['chi2x'], p['chi2z'] = self.ROTATEY(p['inclination'], p['chi2x'], p['chi2z'])

        print("spins after rotation by {0} = ".format(p['inclination']))
        print("chi1x = {0}, chi1y = {1}, chi1z = {2}".format(p['chi1x'], p['chi1y'], p['chi1z']))
        print("chi2x = {0}, chi2y = {1}, chi2z = {2}".format(p['chi2x'], p['chi2y'], p['chi2z']))



        #from this we construct the orbital angular momentum
        #Again, this is a rotation about the y-axis.
        lnhatx = sin(p['inclination'])
        lnhaty = 0.
        lnhatz = cos(p['inclination'])

        print("lnhatx = {0}, lnhaty = {1}, lnhatz = {2}".format(lnhatx, lnhaty, lnhatz))

        #compute the aligned spin component. The component of the spins along lnhat
        chi1_l = dot( [p['chi1x'], p['chi1y'], p['chi1z']], [lnhatx, lnhaty, lnhatz] )
        chi2_l = dot( [p['chi2x'], p['chi2y'], p['chi2z']], [lnhatx, lnhaty, lnhatz] )

        #compute component of spins perpendicular to lnhat
        chi1_perp_x = p['chi1x'] - chi1_l * lnhatx
        chi1_perp_y = p['chi1y'] - chi1_l * lnhaty
        chi1_perp_z = p['chi1z'] - chi1_l * lnhatz

        chi2_perp_x = p['chi2x'] - chi2_l * lnhatx
        chi2_perp_y = p['chi2y'] - chi2_l * lnhaty
        chi2_perp_z = p['chi2z'] - chi2_l * lnhatz

        #magnitude of in-plane dimensionless spins
        chi1_perp = norm( [chi1_perp_x, chi1_perp_y, chi1_perp_z] )
        chi2_perp = norm( [chi2_perp_x, chi2_perp_y, chi2_perp_z] )

        print("chi1_perp = {0}".format(chi1_perp))
        print("chi2_perp = {0}".format(chi2_perp))

        print("m1_2 = {0}".format(m1_2))
        print("m2_2 = {0}".format(m2_2))

        #magnitude of in-plane Dimensionfull spins
        S1_perp = chi1_perp * m1_2
        S2_perp = chi2_perp * m2_2
        print("S1_perp = ", S1_perp)
        print("S2_perp = ", S2_perp)

        #Below Eq.3.2  PhysRevD.91.024043
        A1 = 2 + (3*p['m2']) / (2*p['m1'])
        A2 = 2 + (3*p['m1']) / (2*p['m2'])
        print("A1 = ", A1)
        print("A2 = ", A2)
        ASp1 = A1*S1_perp
        ASp2 = A2*S2_perp
        print("ASp1 = ", ASp1)
        print("ASp2 = ", ASp2)

        num = max([ASp1, ASp2])
        den = A1*m1_2
        print("num = ", num)
        print("den = ", den)
        #chip = max(A1 Sp1, A2 Sp2) / (A_i m_i^2) The denomenator references the larger BH
        chip = num / den

        #compute L, J0 and orientation angles
        piM = Constants.LAL_PI * p['M_sec']
        v_ref = (piM * p['fRef'])**(1./3.)

        #Use 2PN approximation for initial L
        #magnitude of L
        L0 = p['Mtot']**2. * PhenomPL2PN(v_ref, p['eta'])

        #compute initial J
        #NOTE: we the spins need to be dimensionfull
        Jx0 = L0 * lnhatx + p['chi1x']*m1_2 + p['chi2x']*m2_2
        Jy0 = L0 * lnhaty + p['chi1y']*m1_2 + p['chi2y']*m2_2
        Jz0 = L0 * lnhatz + p['chi1z']*m1_2 + p['chi2z']*m2_2
        J0 = norm( [ Jx0, Jy0, Jz0 ] )

        #Compute thetaJ, the angle between J0 and line of sight (z-direction)
        if (J0 < 1e-10):
            print("Warning: |J0| < 1e-10. Setting thetaJ = 0.\n")
            thetaJ = 0.
        else:
            thetaJ = arccos(Jz0 / J0)

        #phiJ, We only use this angle internally since it is degenerate with alpha0.
        #NOTE:
        #in C code
        #if (Jx0 < DBL_MIN && Jy0 < DBL_MIN)
        #I think the replacement is the same
        if (Jx0 <= 0. and Jy0 <= 0.):
            phiJ = 0.
        else:
            phiJ = arctan2(Jy0, Jx0) #Angle of J0 in the plane of the sky
        #NOTE: Compared to the similar code in SpinTaylorF2 we have defined phiJ as the angle between the positive
        #(rather than the negative) x-axis and the projection of J0, since this is a more natural definition of the angle.
        #We have also renamed the angle from psiJ to phiJ.

        #Rotate Lnhat back to frame where J is along z and the line of
        #sight in the Oxz plane with >0 projection in x, to figure out initial alpha
        #The rotation matrix is
        #{
        #{-cos(thetaJ)*cos(phiJ), -cos(thetaJ)*sin(phiJ), sin(thetaJ)},
        #{sin(phiJ), -cos(phiJ), 0},
        #{cos(phiJ)*sin(thetaJ), sin(thetaJ)*sin(phiJ),cos(thetaJ)}
        #}

        rotLx = -lnhatx*cos(thetaJ)*cos(phiJ) - lnhaty*cos(thetaJ)*sin(phiJ) + lnhatz*sin(thetaJ)
        rotLy = lnhatx*sin(phiJ) - lnhaty*cos(phiJ)
        if (rotLx == 0.0 and rotLy == 0.0):
            alpha0 = 0.0
        else:
            alpha0 = arctan2(rotLy, rotLx)

        print {"chi1_l" : chi1_l, "chi2_l" : chi2_l, "chip": chip, "thetaJ" : thetaJ, "alpha0" : alpha0}

        return {"chi1_l" : chi1_l, "chi2_l" : chi2_l, "chip": chip, "thetaJ" : thetaJ, "alpha0" : alpha0}

    def ROTATEY(self, angle, vx, vz):
        """LALSimInspiralSpinTaylor.c
        Rotate components of a vector aboyt y axis
        y-component doesn't change so we don't include it.
        It's just the usual rotation matrix."""
    	tmp1 = vx*cos(angle) + vz*sin(angle)
    	tmp2 = - vx*sin(angle) + vz*cos(angle)
    	vx = tmp1
    	vz = tmp2
        return vx, vz






#
