from phenom.waveform.waveform import Waveform
from phenom.utils.utils import M_eta_m1_m2, chipn, UsefulPowers, __MTSUN_SI__, __MSUN_SI__, pow_2_of, pow_3_of
from phenom.utils.remnant import fring, fdamp, FinalSpin0815
from numpy import sqrt, pi, arange, zeros, exp, fabs

class PhenomDInternalsAmplitude(object):
    """docstring for PhenomDInternalsAmplitude"""
    def __init__(self):
        pass
    # ///////////////////////////// Amplitude: Inspiral functions /////////////////////////

    # // Phenom coefficients rho1, ..., rho3 from direct fit
    # // AmpInsDFFitCoeffChiPNFunc[eta, chiPN]

    def rho1_fun(self, p):
        """
        rho_1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 3931.8979897196696 - 17395.758706812805*eta \
        + (3132.375545898835 + 343965.86092361377*eta - 1.2162565819981997e6*eta2)*xi \
        + (-70698.00600428853 + 1.383907177859705e6*eta - 3.9662761890979446e6*eta2)*xi2 \
        + (-60017.52423652596 + 803515.1181825735*eta - 2.091710365941658e6*eta2)*xi3

    def rho2_fun(self, p):
        """
        rho_2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return -40105.47653771657 + 112253.0169706701*eta \
        + (23561.696065836168 - 3.476180699403351e6*eta + 1.137593670849482e7*eta2)*xi \
        + (754313.1127166454 - 1.308476044625268e7*eta + 3.6444584853928134e7*eta2)*xi2 \
        + (596226.612472288 - 7.4277901143564405e6*eta + 1.8928977514040343e7*eta2)*xi3

    def rho3_fun(self, p):
        """
        rho_3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 83208.35471266537 - 191237.7264145924*eta + \
        (-210916.2454782992 + 8.71797508352568e6*eta - 2.6914942420669552e7*eta2)*xi \
        + (-1.9889806527362722e6 + 3.0888029960154563e7*eta - 8.390870279256162e7*eta2)*xi2 \
        + (-1.4535031953446497e6 + 1.7063528990822166e7*eta - 4.2748659731120914e7*eta2)*xi3

    def init_amp_ins_prefactors(self, p, powers_of_pi):
        """helper function for AmpInsAnsatz
        p : dict
        powers_of_pi : instant of UsefulPowers class
        output:
        output a dictionary called prefactors"""
        from phenom.utils.utils import amp0Func

        prefactors = {}

        eta = p['eta']

        prefactors['amp0'] = amp0Func(p['eta'])

        chi1 = p['chi1z']
        chi2 = p['chi2z']
        # phenom coefficients
        rho1 = p['rho1']
        rho2 = p['rho2']
        rho3 = p['rho3']

        chi12 = chi1*chi1
        chi22 = chi2*chi2
        eta2 = eta*eta
        eta3 = eta*eta2

        Pi = pi
        Pi2 = powers_of_pi.two
        Seta = sqrt(1.0 - 4.0*eta)

        prefactors['two_thirds'] = ((-969 + 1804*eta)*powers_of_pi.two_thirds)/672.
        prefactors['one'] = ((chi1*(81*(1 + Seta) - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Pi)/48.
        prefactors['four_thirds'] = (	(-27312085.0 - 10287648*chi22 - 10287648*chi12*(1 + Seta) + 10287648*chi22*Seta
        							 + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta
        							 + 35371056*eta2
        							 )
        						* powers_of_pi.four_thirds) / 8.128512e6
        prefactors['five_thirds'] = (powers_of_pi.five_thirds * (chi2*(-285197*(-1 + Seta) + 4*(-91902 + 1579*Seta)*eta - 35632*eta2)
        														+ chi1*(285197*(1 + Seta) - 4*(91902 + 1579*Seta)*eta - 35632*eta2)
        														+ 42840*(-1.0 + 4*eta)*Pi
        														)
        							) / 32256.
        prefactors['two'] = - (Pi2*(-336*(-3248849057.0 + 2943675504*chi12 - 3339284256*chi1*chi2 + 2943675504*chi22)*eta2
        						  - 324322727232*eta3
        						  - 7*(-177520268561 + 107414046432*chi22 + 107414046432*chi12*(1 + Seta)
        								- 107414046432*chi22*Seta + 11087290368*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi
        								)
        						  + 12*eta*(-545384828789 - 176491177632*chi1*chi2 + 202603761360*chi22
        									+ 77616*chi12*(2610335 + 995766*Seta) - 77287373856*chi22*Seta
        									+ 5841690624*(chi1 + chi2)*Pi + 21384760320*Pi2
        									)
        							)
        					)/6.0085960704e10
        prefactors['seven_thirds']= rho1
        prefactors['eight_thirds'] = rho2
        prefactors['three'] = rho3

        return prefactors

    def AmpInsAnsatz(self, Mf, powers_of_Mf, prefactors):
        # The Newtonian term in LAL is fine and we should use exactly the same (either hardcoded or call).
        # We just use the Mathematica expression for convenience.
        """
        input:
            Mf : float
                dimensionless frequency
            powers_of_Mf : instance of UsefulPowers class
            prefactors : dict
                output from init_amp_ins_prefactors function
        Inspiral amplitude plus rho phenom coefficents. rho coefficients computed
        in rho1_fun, rho2_fun, rho3_fun functions.
        Amplitude is a re-expansion. See 1508.07253 and Equation 29, 30 and Appendix B arXiv:1508.07253 for details
        """
        Mf2 = powers_of_Mf.two
        Mf3 = Mf*Mf2

        return 1 + powers_of_Mf.two_thirds * prefactors['two_thirds'] \
        		+ Mf * prefactors['one'] + powers_of_Mf.four_thirds * prefactors['four_thirds'] \
        		+ powers_of_Mf.five_thirds * prefactors['five_thirds'] + Mf2 * prefactors['two'] \
        		+ powers_of_Mf.seven_thirds * prefactors['seven_thirds'] \
                + powers_of_Mf.eight_thirds * prefactors['eight_thirds'] \
        		+ Mf3 * prefactors['three']

    def DAmpInsAnsatz(self, Mf, p, powers_of_pi, powers_of_Mf):
        """
        input frequency : Mf
        Take the AmpInsAnsatz expression and compute the first derivative
        with respect to frequency to get the expression below.
        """
        eta = p['eta']
        chi1 = p['chi1z']
        chi2 = p['chi2z']
        rho1 = p['rho1']
        rho2 = p['rho2']
        rho3 = p['rho3']

        chi12 = chi1*chi1
        chi22 = chi2*chi2
        eta2 = eta*eta
        eta3 = eta*eta2
        Mf2 = Mf*Mf
        Pi = pi
        Pi2 = powers_of_pi.two
        Seta = sqrt(1.0 - 4.0*eta)

        return ((-969 + 1804*eta)*powers_of_pi.two_thirds)/(1008.*powers_of_Mf.third) \
        + ((chi1*(81*(1 + Seta) - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Pi)/48. \
        + ((-27312085 - 10287648*chi22 - 10287648*chi12*(1 + Seta) \
        + 10287648*chi22*Seta + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta \
        + 35371056*eta2)*powers_of_Mf.third*powers_of_pi.four_thirds)/6.096384e6 \
        + (5*powers_of_Mf.two_thirds*powers_of_pi.five_thirds*(chi2*(-285197*(-1 + Seta) \
        + 4*(-91902 + 1579*Seta)*eta - 35632*eta2) + chi1*(285197*(1 + Seta) \
        - 4*(91902 + 1579*Seta)*eta - 35632*eta2) + 42840*(-1 + 4*eta)*Pi))/96768. \
        - (Mf*Pi2*(-336*(-3248849057.0 + 2943675504*chi12 - 3339284256*chi1*chi2 + 2943675504*chi22)*eta2 - 324322727232*eta3 \
        - 7*(-177520268561 + 107414046432*chi22 + 107414046432*chi12*(1 + Seta) - 107414046432*chi22*Seta \
        + 11087290368*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi) \
        + 12*eta*(-545384828789.0 - 176491177632*chi1*chi2 + 202603761360*chi22 + 77616*chi12*(2610335 + 995766*Seta) \
        - 77287373856*chi22*Seta + 5841690624*(chi1 + chi2)*Pi + 21384760320*Pi2)))/3.0042980352e10 \
        + (7.0/3.0)*pow(Mf,4.0/3.0)*rho1 + (8.0/3.0)*powers_of_Mf.five_thirds*rho2 + 3*Mf2*rho3

# /////////////////////////// Amplitude: Merger-Ringdown functions ///////////////////////
#
# // Phenom coefficients gamma1, ..., gamma3
# // AmpMRDAnsatzFunc[]

    def gamma1_fun(self, p):
        """
        gamma 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 0.006927402739328343 + 0.03020474290328911*eta \
        + (0.006308024337706171 - 0.12074130661131138*eta + 0.26271598905781324*eta2)*xi \
        + (0.0034151773647198794 - 0.10779338611188374*eta + 0.27098966966891747*eta2)*xi2 \
        + (0.0007374185938559283 - 0.02749621038376281*eta + 0.0733150789135702*eta2)*xi3


    def gamma2_fun(self, p):
        """
        gamma 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 1.010344404799477 + 0.0008993122007234548*eta \
        + (0.283949116804459 - 4.049752962958005*eta + 13.207828172665366*eta2)*xi \
        + (0.10396278486805426 - 7.025059158961947*eta + 24.784892370130475*eta2)*xi2 \
        + (0.03093202475605892 - 2.6924023896851663*eta + 9.609374464684983*eta2)*xi3


    def gamma3_fun(self, p):
        """
        gamma 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 1.3081615607036106 - 0.005537729694807678*eta \
        + (-0.06782917938621007 - 0.6689834970767117*eta + 3.403147966134083*eta2)*xi \
        + (-0.05296577374411866 - 0.9923793203111362*eta + 4.820681208409587*eta2)*xi2 \
        + (-0.006134139870393713 - 0.38429253308696365*eta + 1.7561754421985984*eta2)*xi3

    def AmpMRDAnsatz(self, Mf, p):
        """
        input frequency : Mf
        Ansatz for the merger-ringdown amplitude. Equation 19 arXiv:1508.07253
        """
        fRD = p['fRD']
        fDM = p['fDM']
        gamma1 = p['gamma1']
        gamma2 = p['gamma2']
        gamma3 = p['gamma3']
        fDMgamma3 = fDM*gamma3
        fminfRD = Mf - fRD
        return exp( -(fminfRD)*gamma2 / (fDMgamma3) ) \
        * (fDMgamma3*gamma1) / (pow_2_of(fminfRD) + pow_2_of(fDMgamma3))


    def DAmpMRDAnsatz(self, Mf, p):
        """
        input frequency : Mf
        first frequency derivative of AmpMRDAnsatz
        """
        fRD = p['fRD']
        fDM = p['fDM']
        gamma1 = p['gamma1']
        gamma2 = p['gamma2']
        gamma3 = p['gamma3']

        fDMgamma3 = fDM * gamma3
        pow2_fDMgamma3 = pow_2_of(fDMgamma3)
        fminfRD = Mf - fRD
        expfactor = exp(((fminfRD)*gamma2)/(fDMgamma3))
        pow2pluspow2 = pow_2_of(fminfRD) + pow2_fDMgamma3

        return (-2*fDM*(fminfRD)*gamma3*gamma1) / ( expfactor * pow_2_of(pow2pluspow2)) \
        - (gamma2*gamma1) / ( expfactor * (pow2pluspow2))


    def fmaxCalc(self, p):
        """
        Equation 20 arXiv:1508.07253 (called f_peak in paper)
        analytic location of maximum of AmpMRDAnsatz
        """
        fRD    = p['fRD']
        fDM    = p['fDM']
        gamma2 = p['gamma2']
        gamma3 = p['gamma3']

        # // NOTE: There's a problem with this expression from the paper becoming imaginary if gamma2>=1
        # // Fix: if gamma2 >= 1 then set the square root term to zero.
        if (gamma2 <= 1):
            return fabs(fRD + (fDM*(-1 + sqrt(1 - pow_2_of(gamma2)))*gamma3)/gamma2)
        else:
            return fabs(fRD + (fDM*(-1)*gamma3)/gamma2)

    # ///////////////////////////// Amplitude: Intermediate functions ////////////////////////

    # // Phenom coefficients delta0, ..., delta4 determined from collocation method
    # // (constraining 3 values and 2 derivatives)
    # // AmpIntAnsatzFunc[]

    def AmpIntAnsatz(self, Mf, p):
        """
        Ansatz for the intermediate amplitude. Equation 21 arXiv:1508.07253
        """
        Mf2 = Mf*Mf
        Mf3 = Mf*Mf2
        Mf4 = Mf*Mf3
        return p['delta0'] + p['delta1']*Mf + p['delta2']*Mf2 + p['delta3']*Mf3 + p['delta4']*Mf4

    def AmpIntColFitCoeff(self, p):
        """
        The function name stands for 'Amplitude Intermediate Collocation Fit Coefficient'
        This is the 'v2' value in Table 5 of arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 0.8149838730507785 + 2.5747553517454658*eta \
        + (1.1610198035496786 - 2.3627771785551537*eta + 6.771038707057573*eta2)*xi \
        + (0.7570782938606834 - 2.7256896890432474*eta + 7.1140380397149965*eta2)*xi2 \
        + (0.1766934149293479 - 0.7978690983168183*eta + 2.1162391502005153*eta2)*xi3

    def delta0_fun(self, p, d):
        """
        p : dict
        d : dict of delta
        The following functions (delta{0,1,2,3,4}_fun) were derived
        in mathematica according to
        the constraints detailed in arXiv:1508.07253,
        section 'Region IIa - intermediate'.
        These are not given in the paper.
        Can be rederived by solving Equation 21 for the constraints
        given in Equations 22-26 in arXiv:1508.07253
        """
        f1 = p['f1']
        f2 = p['f2']
        f3 = p['f3']
        v1 = p['v1']
        v2 = p['v2']
        v3 = p['v3']
        d1 = p['d1']
        d2 = p['d2']

        f12 = d['f12']
        f13 = d['f13']
        f14 = d['f14']
        f15 = d['f15']
        f22 = d['f22']
        f23 = d['f23']
        f24 = d['f24']
        f32 = d['f32']
        f33 = d['f33']
        f34 = d['f34']
        f35 = d['f35']

        return -((d2*f15*f22*f3 - 2*d2*f14*f23*f3 + d2*f13*f24*f3 - d2*f15*f2*f32 + d2*f14*f22*f32 \
        - d1*f13*f23*f32 + d2*f13*f23*f32 + d1*f12*f24*f32 - d2*f12*f24*f32 + d2*f14*f2*f33 \
        + 2*d1*f13*f22*f33 - 2*d2*f13*f22*f33 - d1*f12*f23*f33 + d2*f12*f23*f33 - d1*f1*f24*f33 \
        - d1*f13*f2*f34 - d1*f12*f22*f34 + 2*d1*f1*f23*f34 + d1*f12*f2*f35 - d1*f1*f22*f35 \
        + 4*f12*f23*f32*v1 - 3*f1*f24*f32*v1 - 8*f12*f22*f33*v1 + 4*f1*f23*f33*v1 + f24*f33*v1 \
        + 4*f12*f2*f34*v1 + f1*f22*f34*v1 - 2*f23*f34*v1 - 2*f1*f2*f35*v1 + f22*f35*v1 - f15*f32*v2 \
        + 3*f14*f33*v2 - 3*f13*f34*v2 + f12*f35*v2 - f15*f22*v3 + 2*f14*f23*v3 - f13*f24*v3 \
        + 2*f15*f2*f3*v3 - f14*f22*f3*v3 - 4*f13*f23*f3*v3 + 3*f12*f24*f3*v3 - 4*f14*f2*f32*v3 \
        + 8*f13*f22*f32*v3 - 4*f12*f23*f32*v3) / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(f3-f2)))

    def delta1_fun(self, p, d):
        f1 = p['f1']
        f2 = p['f2']
        f3 = p['f3']
        v1 = p['v1']
        v2 = p['v2']
        v3 = p['v3']
        d1 = p['d1']
        d2 = p['d2']

        f12 = d['f12']
        f13 = d['f13']
        f14 = d['f14']
        f15 = d['f15']
        f22 = d['f22']
        f23 = d['f23']
        f24 = d['f24']
        f32 = d['f32']
        f33 = d['f33']
        f34 = d['f34']
        f35 = d['f35']

        return -((-(d2*f15*f22) + 2*d2*f14*f23 - d2*f13*f24 - d2*f14*f22*f3 + 2*d1*f13*f23*f3 \
        + 2*d2*f13*f23*f3 - 2*d1*f12*f24*f3 - d2*f12*f24*f3 + d2*f15*f32 - 3*d1*f13*f22*f32 \
        - d2*f13*f22*f32 + 2*d1*f12*f23*f32 - 2*d2*f12*f23*f32 + d1*f1*f24*f32 + 2*d2*f1*f24*f32 \
        - d2*f14*f33 + d1*f12*f22*f33 + 3*d2*f12*f22*f33 - 2*d1*f1*f23*f33 - 2*d2*f1*f23*f33 \
        + d1*f24*f33 + d1*f13*f34 + d1*f1*f22*f34 - 2*d1*f23*f34 - d1*f12*f35 + d1*f22*f35 \
        - 8*f12*f23*f3*v1 + 6*f1*f24*f3*v1 + 12*f12*f22*f32*v1 - 8*f1*f23*f32*v1 - 4*f12*f34*v1 \
        + 2*f1*f35*v1 + 2*f15*f3*v2 - 4*f14*f32*v2 + 4*f12*f34*v2 - 2*f1*f35*v2 - 2*f15*f3*v3 \
        + 8*f12*f23*f3*v3 - 6*f1*f24*f3*v3 + 4*f14*f32*v3 - 12*f12*f22*f32*v3 + 8*f1*f23*f32*v3) \
        / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(-f2 + f3)))

    def delta2_fun(self, p, d):
        f1 = p['f1']
        f2 = p['f2']
        f3 = p['f3']
        v1 = p['v1']
        v2 = p['v2']
        v3 = p['v3']
        d1 = p['d1']
        d2 = p['d2']

        f12 = d['f12']
        f13 = d['f13']
        f14 = d['f14']
        f15 = d['f15']
        f23 = d['f23']
        f24 = d['f24']
        f32 = d['f32']
        f33 = d['f33']
        f34 = d['f34']
        f35 = d['f35']

        return -((d2*f15*f2 - d1*f13*f23 - 3*d2*f13*f23 + d1*f12*f24 + 2*d2*f12*f24 - d2*f15*f3 \
        + d2*f14*f2*f3 - d1*f12*f23*f3 + d2*f12*f23*f3 + d1*f1*f24*f3 - d2*f1*f24*f3 - d2*f14*f32 \
        + 3*d1*f13*f2*f32 + d2*f13*f2*f32 - d1*f1*f23*f32 + d2*f1*f23*f32 - 2*d1*f24*f32 - d2*f24*f32 \
        - 2*d1*f13*f33 + 2*d2*f13*f33 - d1*f12*f2*f33 - 3*d2*f12*f2*f33 + 3*d1*f23*f33 + d2*f23*f33 \
        + d1*f12*f34 - d1*f1*f2*f34 + d1*f1*f35 - d1*f2*f35 + 4*f12*f23*v1 - 3*f1*f24*v1 + 4*f1*f23*f3*v1 \
        - 3*f24*f3*v1 - 12*f12*f2*f32*v1 + 4*f23*f32*v1 + 8*f12*f33*v1 - f1*f34*v1 - f35*v1 - f15*v2 \
        - f14*f3*v2 + 8*f13*f32*v2 - 8*f12*f33*v2 + f1*f34*v2 + f35*v2 + f15*v3 - 4*f12*f23*v3 + 3*f1*f24*v3 \
        + f14*f3*v3 - 4*f1*f23*f3*v3 + 3*f24*f3*v3 - 8*f13*f32*v3 + 12*f12*f2*f32*v3 - 4*f23*f32*v3) \
        / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(-f2 + f3)))

    def delta3_fun(self, p, d):
        f1 = p['f1']
        f2 = p['f2']
        f3 = p['f3']
        v1 = p['v1']
        v2 = p['v2']
        v3 = p['v3']
        d1 = p['d1']
        d2 = p['d2']

        f12 = d['f12']
        f13 = d['f13']
        f14 = d['f14']
        f22 = d['f22']
        f24 = d['f24']
        f32 = d['f32']
        f33 = d['f33']
        f34 = d['f34']

        return -((-2*d2*f14*f2 + d1*f13*f22 + 3*d2*f13*f22 - d1*f1*f24 - d2*f1*f24 + 2*d2*f14*f3 \
        - 2*d1*f13*f2*f3 - 2*d2*f13*f2*f3 + d1*f12*f22*f3 - d2*f12*f22*f3 + d1*f24*f3 + d2*f24*f3 \
        + d1*f13*f32 - d2*f13*f32 - 2*d1*f12*f2*f32 + 2*d2*f12*f2*f32 + d1*f1*f22*f32 - d2*f1*f22*f32 \
        + d1*f12*f33 - d2*f12*f33 + 2*d1*f1*f2*f33 + 2*d2*f1*f2*f33 - 3*d1*f22*f33 - d2*f22*f33 \
        - 2*d1*f1*f34 + 2*d1*f2*f34 - 4*f12*f22*v1 + 2*f24*v1 + 8*f12*f2*f3*v1 - 4*f1*f22*f3*v1 \
        - 4*f12*f32*v1 + 8*f1*f2*f32*v1 - 4*f22*f32*v1 - 4*f1*f33*v1 + 2*f34*v1 + 2*f14*v2 \
        - 4*f13*f3*v2 + 4*f1*f33*v2 - 2*f34*v2 - 2*f14*v3 + 4*f12*f22*v3 - 2*f24*v3 + 4*f13*f3*v3 \
        - 8*f12*f2*f3*v3 + 4*f1*f22*f3*v3 + 4*f12*f32*v3 - 8*f1*f2*f32*v3 + 4*f22*f32*v3) \
        / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(-f2 + f3)))

    def delta4_fun(self, p, d):
        f1 = p['f1']
        f2 = p['f2']
        f3 = p['f3']
        v1 = p['v1']
        v2 = p['v2']
        v3 = p['v3']
        d1 = p['d1']
        d2 = p['d2']

        f12 = d['f12']
        f13 = d['f13']
        f22 = d['f22']
        f23 = d['f23']
        f32 = d['f32']
        f33 = d['f33']

        return -((d2*f13*f2 - d1*f12*f22 - 2*d2*f12*f22 + d1*f1*f23 + d2*f1*f23 - d2*f13*f3 + 2*d1*f12*f2*f3 \
        + d2*f12*f2*f3 - d1*f1*f22*f3 + d2*f1*f22*f3 - d1*f23*f3 - d2*f23*f3 - d1*f12*f32 + d2*f12*f32 \
        - d1*f1*f2*f32 - 2*d2*f1*f2*f32 + 2*d1*f22*f32 + d2*f22*f32 + d1*f1*f33 - d1*f2*f33 + 3*f1*f22*v1 \
        - 2*f23*v1 - 6*f1*f2*f3*v1 + 3*f22*f3*v1 + 3*f1*f32*v1 - f33*v1 - f13*v2 + 3*f12*f3*v2 - 3*f1*f32*v2 \
        + f33*v2 + f13*v3 - 3*f1*f22*v3 + 2*f23*v3 - 3*f12*f3*v3 + 6*f1*f2*f3*v3 - 3*f22*f3*v3) \
        / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(-f2 + f3)))

    def ComputeDeltasFromCollocation(self, p, amp_prefactors, powers_of_pi):
        """
        Output: None, It adds
        f1,f2,f3,v1,v2,v3,d1,d2
        and
        delta0,delta1,delta2,delta3,delta4
        to the 'p' dict
        Calculates delta_i's
        Method described in arXiv:1508.07253 section 'Region IIa - intermediate'
        """
        eta = p['eta']
        chi = p['chipn']
        # // Three evenly spaced collocation points in the interval [f1,f3].
        f1 = self.AMP_fJoin_INS # in Mf
        f3 = p['fmaxCalc'] # in Mf
        dfx = (f3 - f1)/2.0 # in Mf
        f2 = f1 + dfx # in Mf

        powers_of_f1 = UsefulPowers(f1)

        # // v1 is inspiral model evaluated at f1
        # // d1 is derivative of inspiral model evaluated at f1
        v1 = self.AmpInsAnsatz(f1, powers_of_f1, amp_prefactors)
        d1 = self.DAmpInsAnsatz(f1, p, powers_of_pi, powers_of_f1)

        # // v3 is merger-ringdown model evaluated at f3
        # // d2 is derivative of merger-ringdown model evaluated at f3
        v3 = self.AmpMRDAnsatz(f3, p)
        d2 = self.DAmpMRDAnsatz(f3, p)

        # // v2 is the value of the amplitude evaluated at f2
        # // they come from the fit of the collocation points in the intermediate region
        v2 = self.AmpIntColFitCoeff(p)

        # save these to the `global p dictionary'
        self.p['f1'] = f1
        self.p['f2'] = f2
        self.p['f3'] = f3
        self.p['v1'] = v1
        self.p['v2'] = v2
        self.p['v3'] = v3
        self.p['d1'] = d1
        self.p['d2'] = d2

        # // Now compute the delta_i's from the collocation coefficients
        # // Precompute common quantities here and pass along to delta functions.
        d = {}
        d['f12'] = f1*f1
        d['f13'] = f1*d['f12']
        d['f14'] = f1*d['f13']
        d['f15'] = f1*d['f14']
        d['f22'] = f2*f2
        d['f23'] = f2*d['f22']
        d['f24'] = f2*d['f23']
        d['f32'] = f3*f3
        d['f33'] = f3*d['f32']
        d['f34'] = f3*d['f33']
        d['f35'] = f3*d['f34']
        self.p['delta0'] = self.delta0_fun(self.p, d)
        self.p['delta1'] = self.delta1_fun(self.p, d)
        self.p['delta2'] = self.delta2_fun(self.p, d)
        self.p['delta3'] = self.delta3_fun(self.p, d)
        self.p['delta4'] = self.delta4_fun(self.p, d)

class PhenomDInternals(PhenomDInternalsAmplitude):
    """docstring for PhenomDInternals"""
    def __init__(self):
        pass

# TODO:0 Figure out how to structure PhenomD class
class PhenomD(Waveform, PhenomDInternals):
    """docstring for PhenomD"""
    # Use the constructor from the Waveform Class
    # def __init__(self, **p):
    #     super(PhenomD, self).__init__(p)
    def __init__(self, m1=10. * __MSUN_SI__, m2=10. * __MSUN_SI__, chi1z=0., chi2z=0., f_min=20., f_max=0., delta_f=1/64., distance=1e6, fRef=0., phiRef=0.):
        """
        input:
        m1 (SI)
        m2 (SI)
        chi1z (dimensionless)
        chi2z (dimensionless)
        f_min (Hz)
        f_max (Hz)
        delta_f (sample rate Hz)
        distance (m)
        fRef (reference frequency Hz)
        phiRef (phase at fRef)"""

        # enforce m1 >= m2 and chi1 is on m1
        if m1<m2: # swap spins and masses
            chi1z, chi2z = chi2z, chi1z
            m1, m2 = m2, m1

        # put inputs into a dictionary
        self.p             = {}
        self.p['m1_SI']    = m1
        self.p['m2_SI']    = m2
        self.p['chi1z']    = chi1z
        self.p['chi2z']    = chi2z
        self.p['f_min']    = f_min
        self.p['f_max']    = f_max
        self.p['delta_f']  = delta_f
        self.p['distance'] = distance
        self.p['fRef']     = fRef
        self.p['phiRef']   = phiRef

        # external: SI; internal: solar masses
        self.p['m1'] = self.p['m1_SI'] / __MSUN_SI__
        self.p['m2'] = self.p['m2_SI'] / __MSUN_SI__

        self.p['Mtot'], self.p['eta'] = M_eta_m1_m2(self.p['m1'], self.p['m2'])
        self.p['chipn'] = chipn(self.p['eta'], self.p['chi1z'], self.p['chi2z'])

        self.M_sec = self.p['Mtot'] * __MTSUN_SI__ # Conversion factor Hz -> dimensionless frequency

        self.define_phenomD_constants()

        # Default PhenomD values
        if self.p['f_max'] == 0.: self.p['f_max'] = self.f_CUT / self.M_sec # converted from Mf to Hz

        #initialise UsefulPowers instances
        self.powers_of_pi = UsefulPowers(pi)

        # compute prediction of ringdown frequency
        self.p['finspin'] = FinalSpin0815(self.p['eta'], self.p['chi1z'], self.p['chi2z'])
        self.p['fRD']     = fring(self.p['eta'], self.p['chi1z'], self.p['chi2z'], self.p['finspin'])
        self.p['fDM']     = fdamp(self.p['eta'], self.p['chi1z'], self.p['chi2z'], self.p['finspin'])


        print "self.rho3_fun(self.p) = ", self.rho3_fun(self.p)

        #TODO: populate amplitude and phase phenom coefficient dictionaries
        #inspiral amplitude coeffs
        self.p['rho1'] = self.rho1_fun(self.p)
        self.p['rho2'] = self.rho2_fun(self.p)
        self.p['rho3'] = self.rho3_fun(self.p)

        #inspiral amplitude prefactors
        self.amp_prefactors = super(PhenomDInternals, self).init_amp_ins_prefactors(self.p, self.powers_of_pi)

        # print self.AmpInsAnsatz(0.08, UsefulPowers(0.08), self.amp_prefactors)
        # print self.DAmpInsAnsatz(0.08, self.p, self.powers_of_pi, UsefulPowers(0.08))
        #merger-ringdown amplitude coeffs
        self.p['gamma1'] = self.gamma1_fun(self.p)
        self.p['gamma2'] = self.gamma2_fun(self.p)
        self.p['gamma3'] = self.gamma3_fun(self.p)

        # fmaxCalc is the peak of the merger-ringdown amplitude function
        self.p['fmaxCalc'] = self.fmaxCalc(self.p)

        #intermediate amplitude coeffs
        self.ComputeDeltasFromCollocation(self.p, self.amp_prefactors, self.powers_of_pi)

        #inspiral phase coeffs
        #inspiral phase prefactors
        #intermediate phase coeffs
        #merger-ringdown phase coeffs


        print "self.AmpMRDAnsatz(0.088, self.p) = ",self.AmpMRDAnsatz(0.088, self.p)
        print "self.DAmpMRDAnsatz(0.088, self.p) = ",self.DAmpMRDAnsatz(0.088, self.p)


        # print super(PhenomD, self).rho3_fun(self.p['eta'], self.p['chipn'])
        # print super(PhenomDInternals, self).rho3_fun(self.p['eta'], self.p['chipn'])
        # super(PhenomD, self).interal()
        # super(PhenomDInternals, self).interal()


        # flist = arange(self.p['f_min'], self.p['f_max'], self.p['delta_f'])
        # amp = zeros(len(flist))
        # phi = zeros(len(flist))
        # for i in range(len(flist)):
            # Mf = flist[i] * self.M_sec
            # self.powers_of_Mf = UsefulPowers(Mf)
            # print Mf
            # amp[i] = self.IMRPhenomDAmplitude(Mf, self.p, self.powers_of_Mf, self.amp_prefactors)
            # phi[i] = self.IMRPhenomDPhase(Mf, self.p, self.pn, self.powers_of_Mf, self.phi_prefactors)


    def define_phenomD_constants(self):
        """
        sets phenomD CONSTANTS
        https://github.com/lscsoft/lalsuite/blob/master/lalsimulation/src/LALSimIMRPhenomD.h
        """
        # /* CONSTANTS */

        # /**
        #  * Dimensionless frequency (Mf) at which define the end of the waveform
        #  */
        self.f_CUT = 0.2

        # /**
        #  * Dimensionless frequency (Mf) at which the inspiral amplitude
        #  * switches to the intermediate amplitude
        #  */
        self.AMP_fJoin_INS = 0.014

        # /**
        #  * Dimensionless frequency (Mf) at which the inspiral phase
        #  * switches to the intermediate phase
        #  */
        self.PHI_fJoin_INS = 0.018

        # /**
        #   * Minimal final spin value below which the waveform might behave pathological
        #   * because the ISCO frequency is too low. For more details, see the review wiki
        #   * page https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/PhenD_LargeNegativeSpins
        #   */
        self.MIN_FINAL_SPIN = -0.717

        # /**
        #   * A large mass ratio causes memory over-runs.
        #   * We test and put the limit an order of magnitude above that of previous waveform models (which were around q=100).
        #   */
        self.MAX_ALLOWED_MASS_RATIO = 5000
        pass

    def IMRPhenomDAmplitude(self, Mf, p, powers_of_Mf, amp_prefactors):
        """// Call ComputeIMRPhenomDAmplitudeCoefficients() first!
        This function computes the IMR amplitude given phenom coefficients.
        Defined in VIII. Full IMR Waveforms arXiv:1508.07253
        """
        # Defined in VIII. Full IMR Waveforms arXiv:1508.07253
        # The inspiral, intermediate and merger-ringdown amplitude parts

        # Transition frequencies
        p['fInsJoin'] = self.AMP_fJoin_INS
        p['fMRDJoin'] = p['fmaxCalc']

        AmpPreFac = amp_prefactors['amp0'] / powers_of_Mf.seven_sixths
        # split the calculation to just 1 of 3 possible mutually exclusive ranges
        # this effectively implements a step function transition function
        if (Mf <= p['fInsJoin']):	# Inspiral range
            AmpIns = AmpPreFac * self.AmpInsAnsatz(Mf, powers_of_Mf, amp_prefactors)
            return AmpIns
        elif (Mf >= p['fMRDJoin']):	# MRD range
            AmpMRD = AmpPreFac * self.AmpMRDAnsatz(Mf, p)
            return AmpMRD
        else:
            #	Intermediate range
            AmpInt = AmpPreFac * self.AmpIntAnsatz(Mf, p)
            return AmpInt

    def IMRPhenomDPhase(self, Mf, p, pn, powers_of_Mf, phi_prefactors):
        pass

    def IMRPhenomDGenerateFD(self, Mf):
        """
        generates frequency domain strain"""

        pass
