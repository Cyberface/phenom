from __future__ import division

import sys
from phenom.waveform.waveform import Waveform
from phenom.utils.utils import M_eta_m1_m2, chipn, UsefulPowers, pow_2_of, pow_3_of, pow_4_of, Constants, HztoMf, StoM
from phenom.utils.remnant import fring, fdamp, FinalSpin0815, FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH
from numpy import sqrt, pi, arange, zeros, exp, fabs, log, arctan, angle, unwrap, absolute

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

    def init_amp_ins_prefactors(self, p, model_pars, powers_of_pi):
        """helper function for AmpInsAnsatz
        p : dict
        powers_of_pi : instant of UsefulPowers class
        output:
        output a dictionary called prefactors"""
        from phenom.utils.utils import amp0Func

        prefactors = {}

        eta = p['eta']

        # eq.17 [1508.07253v2]
        prefactors['PNamp0'] = amp0Func(eta)

        chi1 = p['chi1z']
        chi2 = p['chi2z']
        # phenom coefficients
        rho1 = model_pars['rho1']
        rho2 = model_pars['rho2']
        rho3 = model_pars['rho3']

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

    def DAmpInsAnsatz(self, Mf, p, model_pars, powers_of_pi, powers_of_Mf):
        """
        input frequency : Mf
        Take the AmpInsAnsatz expression and compute the first derivative
        with respect to frequency to get the expression below.
        """
        eta = p['eta']
        chi1 = p['chi1z']
        chi2 = p['chi2z']
        rho1 = model_pars['rho1']
        rho2 = model_pars['rho2']
        rho3 = model_pars['rho3']

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

    def AmpMRDAnsatz(self, Mf, model_pars):
        """
        input frequency : Mf
        Ansatz for the merger-ringdown amplitude. Equation 19 arXiv:1508.07253
        """
        fRD = model_pars['fRD']
        fDM = model_pars['fDM']
        gamma1 = model_pars['gamma1']
        gamma2 = model_pars['gamma2']
        gamma3 = model_pars['gamma3']
        fDMgamma3 = fDM*gamma3
        fminfRD = Mf - fRD
        return exp( -(fminfRD)*gamma2 / (fDMgamma3) ) \
        * (fDMgamma3*gamma1) / (pow_2_of(fminfRD) + pow_2_of(fDMgamma3))


    def DAmpMRDAnsatz(self, Mf, model_pars):
        """
        input frequency : Mf
        first frequency derivative of AmpMRDAnsatz
        """
        fRD = model_pars['fRD']
        fDM = model_pars['fDM']
        gamma1 = model_pars['gamma1']
        gamma2 = model_pars['gamma2']
        gamma3 = model_pars['gamma3']

        fDMgamma3 = fDM * gamma3
        pow2_fDMgamma3 = pow_2_of(fDMgamma3)
        fminfRD = Mf - fRD
        expfactor = exp(((fminfRD)*gamma2)/(fDMgamma3))
        pow2pluspow2 = pow_2_of(fminfRD) + pow2_fDMgamma3

        return (-2*fDM*(fminfRD)*gamma3*gamma1) / ( expfactor * pow_2_of(pow2pluspow2)) \
        - (gamma2*gamma1) / ( expfactor * (pow2pluspow2))


    def fmaxCalc(self, model_pars):
        """
        Equation 20 arXiv:1508.07253 (called f_peak in paper)
        analytic location of maximum of AmpMRDAnsatz
        """
        fRD    = model_pars['fRD']
        fDM    = model_pars['fDM']
        gamma2 = model_pars['gamma2']
        gamma3 = model_pars['gamma3']

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

    def AmpIntAnsatz(self, Mf, model_pars):
        """
        Ansatz for the intermediate amplitude. Equation 21 arXiv:1508.07253
        """
        Mf2 = Mf*Mf
        Mf3 = Mf*Mf2
        Mf4 = Mf*Mf3
        return model_pars['delta0'] + model_pars['delta1']*Mf + model_pars['delta2']*Mf2 + model_pars['delta3']*Mf3 + model_pars['delta4']*Mf4

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

    def delta0_fun(self, delta_params, d):
        """
        delta_params
        d : dict of delta
        The following functions (delta{0,1,2,3,4}_fun) were derived
        in mathematica according to
        the constraints detailed in arXiv:1508.07253,
        section 'Region IIa - intermediate'.
        These are not given in the paper.
        Can be rederived by solving Equation 21 for the constraints
        given in Equations 22-26 in arXiv:1508.07253
        """
        f1 = delta_params['f1']
        f2 = delta_params['f2']
        f3 = delta_params['f3']
        v1 = delta_params['v1']
        v2 = delta_params['v2']
        v3 = delta_params['v3']
        d1 = delta_params['d1']
        d2 = delta_params['d2']

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

    def delta1_fun(self, delta_params, d):
        f1 = delta_params['f1']
        f2 = delta_params['f2']
        f3 = delta_params['f3']
        v1 = delta_params['v1']
        v2 = delta_params['v2']
        v3 = delta_params['v3']
        d1 = delta_params['d1']
        d2 = delta_params['d2']

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

    def delta2_fun(self, delta_params, d):
        f1 = delta_params['f1']
        f2 = delta_params['f2']
        f3 = delta_params['f3']
        v1 = delta_params['v1']
        v2 = delta_params['v2']
        v3 = delta_params['v3']
        d1 = delta_params['d1']
        d2 = delta_params['d2']

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

    def delta3_fun(self, delta_params, d):
        f1 = delta_params['f1']
        f2 = delta_params['f2']
        f3 = delta_params['f3']
        v1 = delta_params['v1']
        v2 = delta_params['v2']
        v3 = delta_params['v3']
        d1 = delta_params['d1']
        d2 = delta_params['d2']

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

    def delta4_fun(self, delta_params, d):
        f1 = delta_params['f1']
        f2 = delta_params['f2']
        f3 = delta_params['f3']
        v1 = delta_params['v1']
        v2 = delta_params['v2']
        v3 = delta_params['v3']
        d1 = delta_params['d1']
        d2 = delta_params['d2']

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

    def ComputeDeltasFromCollocation(self, p, model_pars, powers_of_pi):
        """
        Output: diction of:
        f1,f2,f3,v1,v2,v3,d1,d2
        and
        delta0,delta1,delta2,delta3,delta4.
        In phenomD code this gets added to the model_pars dictionary
        Calculates delta_i's
        Method described in arXiv:1508.07253 section 'Region IIa - intermediate'
        """

        delta_params = {}

        eta = p['eta']
        chi = p['chipn']
        # // Three evenly spaced collocation points in the interval [f1,f3].
        f1 = self.AMP_fJoin_INS # in Mf
        f3 = model_pars['fmaxCalc'] # in Mf
        dfx = (f3 - f1)/2.0 # in Mf
        f2 = f1 + dfx # in Mf

        powers_of_f1 = UsefulPowers(f1)

        # // v1 is inspiral model evaluated at f1
        # // d1 is derivative of inspiral model evaluated at f1
        v1 = self.AmpInsAnsatz(f1, powers_of_f1, model_pars['amp_prefactors'])
        d1 = self.DAmpInsAnsatz(f1, p, model_pars, powers_of_pi, powers_of_f1)

        # // v3 is merger-ringdown model evaluated at f3
        # // d2 is derivative of merger-ringdown model evaluated at f3
        v3 = self.AmpMRDAnsatz(f3, model_pars)
        d2 = self.DAmpMRDAnsatz(f3, model_pars)

        # // v2 is the value of the amplitude evaluated at f2
        # // they come from the fit of the collocation points in the intermediate region
        v2 = self.AmpIntColFitCoeff(p)

        # save these to the delta_params dictionary
        delta_params['f1'] = f1
        delta_params['f2'] = f2
        delta_params['f3'] = f3
        delta_params['v1'] = v1
        delta_params['v2'] = v2
        delta_params['v3'] = v3
        delta_params['d1'] = d1
        delta_params['d2'] = d2

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
        delta_params['delta0'] = self.delta0_fun(delta_params, d)
        delta_params['delta1'] = self.delta1_fun(delta_params, d)
        delta_params['delta2'] = self.delta2_fun(delta_params, d)
        delta_params['delta3'] = self.delta3_fun(delta_params, d)
        delta_params['delta4'] = self.delta4_fun(delta_params, d)
        return delta_params

class PhenomDInternalsPhase(object):
    """docstring for PhenomDInternalsPhase"""
    def __init__(self):
        pass

    # /********************************* Phase functions *********************************/
    # ////////////////////////////// Phase: Ringdown functions ///////////////////////////

    # // alpha_i i=1,2,3,4,5 are the phenomenological intermediate coefficients depending on eta and chiPN
    # // PhiRingdownAnsatz is the ringdown phasing in terms of the alpha_i coefficients

    def alpha1Fit(self, p):
        """
        alpha 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 43.31514709695348 + 638.6332679188081*eta \
        + (-32.85768747216059 + 2415.8938269370315*eta - 5766.875169379177*eta2)*xi \
        + (-61.85459307173841 + 2953.967762459948*eta - 8986.29057591497*eta2)*xi2 \
        + (-21.571435779762044 + 981.2158224673428*eta - 3239.5664895930286*eta2)*xi3

    def alpha2Fit(self, p):
        """
        alpha 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return -0.07020209449091723 - 0.16269798450687084*eta \
        + (-0.1872514685185499 + 1.138313650449945*eta - 2.8334196304430046*eta2)*xi \
        + (-0.17137955686840617 + 1.7197549338119527*eta - 4.539717148261272*eta2)*xi2 \
        + (-0.049983437357548705 + 0.6062072055948309*eta - 1.682769616644546*eta2)*xi3

    def alpha3Fit(self, p):
        """
        alpha 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 9.5988072383479 - 397.05438595557433*eta \
        + (16.202126189517813 - 1574.8286986717037*eta + 3600.3410843831093*eta2)*xi \
        + (27.092429659075467 - 1786.482357315139*eta + 5152.919378666511*eta2)*xi2 \
        + (11.175710130033895 - 577.7999423177481*eta + 1808.730762932043*eta2)*xi3

    def alpha4Fit(self, p):
        """
        alpha 4 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return -0.02989487384493607 + 1.4022106448583738*eta \
        + (-0.07356049468633846 + 0.8337006542278661*eta + 0.2240008282397391*eta2)*xi \
        + (-0.055202870001177226 + 0.5667186343606578*eta + 0.7186931973380503*eta2)*xi2 \
        + (-0.015507437354325743 + 0.15750322779277187*eta + 0.21076815715176228*eta2)*xi3

    def alpha5Fit(self, p):
        """
        alpha 5 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 0.9974408278363099 - 0.007884449714907203*eta \
        + (-0.059046901195591035 + 1.3958712396764088*eta - 4.516631601676276*eta2)*xi \
        + (-0.05585343136869692 + 1.7516580039343603*eta - 5.990208965347804*eta2)*xi2 \
        + (-0.017945336522161195 + 0.5965097794825992*eta - 2.0608879367971804*eta2)*xi3

    def PhiMRDAnsatzInt(self, Mf, model_pars, rholm=1.0, taulm=1.0 ):
        """
        Ansatz for the merger-ringdown phase Equation 14 arXiv:1508.07253
        Rholm was added when IMRPhenomHM (high mode) was added.
        Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
        to 1. and PhenomD is recovered.
        Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
        Again, when Taulm = 1.0 then PhenomD is recovered.
        """
        sqrootf = sqrt(Mf);
        fpow1_5 = Mf * sqrootf;
        # // check if this is any faster: 2 sqrts instead of one pow(x,0.75)
        fpow0_75 = sqrt(fpow1_5); # pow(f,0.75)

        ans = -(model_pars['alpha2']/Mf) \
              + (4.0/3.0) * (model_pars['alpha3'] * fpow0_75) \
              + model_pars['alpha1'] * Mf \
              + model_pars['alpha4'] * rholm * arctan( (Mf - model_pars['alpha5'] * model_pars['fRD']) / (rholm * model_pars['fDM'] * taulm) )

        #
        return ans

    def DPhiMRD(self, Mf, model_pars, eta, rholm=1.0, taulm=1.0):
        """
        First frequency derivative of PhiMRDAnsatzInt
        Rholm was added when IMRPhenomHM (high mode) was added.
        Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
        to 1. and PhenomD is recovered.
        Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
        Again, when Taulm = 1.0 then PhenomD is recovered.
        """
        return (model_pars['alpha1'] + model_pars['alpha2']/pow_2_of(Mf) + model_pars['alpha3']/(Mf**0.25) + \
        model_pars['alpha4']/(model_pars['fDM']*taulm*(1 + pow_2_of(Mf - model_pars['alpha5'] * model_pars['fRD'])/pow_2_of(rholm*model_pars['fDM']*taulm)))) / eta

    # ///////////////////////////// Phase: Intermediate functions /////////////////////////////
    #
    # // beta_i i=1,2,3 are the phenomenological intermediate coefficients depending on eta and chiPN
    # // PhiIntAnsatz is the intermediate phasing in terms of the beta_i coefficients
    #
    #
    # // \[Beta]1Fit = PhiIntFitCoeff\[Chi]PNFunc[\[Eta], \[Chi]PN][[1]]

    def beta1Fit(self, p):
        """
        beta 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 97.89747327985583 - 42.659730877489224*eta \
        + (153.48421037904913 - 1417.0620760768954*eta + 2752.8614143665027*eta2)*xi \
        + (138.7406469558649 - 1433.6585075135881*eta + 2857.7418952430758*eta2)*xi2 \
        + (41.025109467376126 - 423.680737974639*eta + 850.3594335657173*eta2)*xi3

    def beta2Fit(self, p):
        """
        beta 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return -3.282701958759534 - 9.051384468245866*eta \
        + (-12.415449742258042 + 55.4716447709787*eta - 106.05109938966335*eta2)*xi \
        + (-11.953044553690658 + 76.80704618365418*eta - 155.33172948098394*eta2)*xi2 \
        + (-3.4129261592393263 + 25.572377569952536*eta - 54.408036707740465*eta2)*xi3

    def beta3Fit(self, p):
        """
        beta 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return -0.000025156429818799565 + 0.000019750256942201327*eta \
        + (-0.000018370671469295915 + 0.000021886317041311973*eta + 0.00008250240316860033*eta2)*xi \
        + (7.157371250566708e-6 - 0.000055780000112270685*eta + 0.00019142082884072178*eta2)*xi2 \
        + (5.447166261464217e-6 - 0.00003220610095021982*eta + 0.00007974016714984341*eta2)*xi3

    def PhiIntAnsatz(self, Mf, model_pars):
        """
        ansatz for the intermediate phase defined by Equation 16 arXiv:1508.07253
        """
        # // 1./eta in paper omitted and put in when need in the functions:
        # // ComputeIMRPhenDPhaseConnectionCoefficients
        # // IMRPhenDPhase
        return  model_pars['beta1']*Mf - model_pars['beta3']/(3.*pow_3_of(Mf)) + model_pars['beta2']*log(Mf)

    def DPhiIntAnsatz(self, Mf, model_pars, eta):
        """
        First frequency derivative of PhiIntAnsatz
        (this time with 1./eta explicitly factored in)
        """
        return (model_pars['beta1'] + model_pars['beta3']/pow_4_of(Mf) + model_pars['beta2']/Mf) / eta

    def DPhiIntTemp(self, Mf, model_pars, eta, C2Int):
        """
        temporary instance of DPhiIntAnsatz used when computing
        coefficients to make the phase C(1) continuous between regions.
        """
        beta1 = model_pars['beta1']
        beta2 = model_pars['beta2']
        beta3 = model_pars['beta3']

        return C2Int + (beta1 + beta3/pow_4_of(Mf) + beta2/Mf)/eta

    # ///////////////////////////// Phase: Inspiral functions /////////////////////////////
    #
    # // sigma_i i=1,2,3,4 are the phenomenological inspiral coefficients depending on eta and chiPN
    # // PhiInsAnsatzInt is a souped up TF2 phasing which depends on the sigma_i coefficients


    def sigma1Fit(self, p):
        """
        sigma 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 2096.551999295543 + 1463.7493168261553*eta \
        + (1312.5493286098522 + 18307.330017082117*eta - 43534.1440746107*eta2)*xi \
        + (-833.2889543511114 + 32047.31997183187*eta - 108609.45037520859*eta2)*xi2 \
        + (452.25136398112204 + 8353.439546391714*eta - 44531.3250037322*eta2)*xi3

    def sigma2Fit(self, p):
        """
        sigma 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return -10114.056472621156 - 44631.01109458185*eta \
        + (-6541.308761668722 - 266959.23419307504*eta + 686328.3229317984*eta2)*xi \
        + (3405.6372187679685 - 437507.7208209015*eta + 1.6318171307344697e6*eta2)*xi2 \
        + (-7462.648563007646 - 114585.25177153319*eta + 674402.4689098676*eta2)*xi3

    def sigma3Fit(self, p):
        """
        sigma 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 22933.658273436497 + 230960.00814979506*eta \
        + (14961.083974183695 + 1.1940181342318142e6*eta - 3.1042239693052764e6*eta2)*xi \
        + (-3038.166617199259 + 1.8720322849093592e6*eta - 7.309145012085539e6*eta2)*xi2 \
        + (42738.22871475411 + 467502.018616601*eta - 3.064853498512499e6*eta2)*xi3

    def sigma4Fit(self, p):
        """
        sigma 4 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        eta = p['eta']
        chi = p['chipn']
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return -14621.71522218357 - 377812.8579387104*eta \
        + (-9608.682631509726 - 1.7108925257214056e6*eta + 4.332924601416521e6*eta2)*xi \
        + (-22366.683262266528 - 2.5019716386377467e6*eta + 1.0274495902259542e7*eta2)*xi2 \
        + (-85360.30079034246 - 570025.3441737515*eta + 4.396844346849777e6*eta2)*xi3



    def SimInspiralTaylorF2AlignedPhasing(self, p):
        """
        input
            p : dict
        output - pn coefficient dictionary
            pfa : dict
            pfa['N'],
            pfa['v[0]'],
            pfa['v[2]'],
            pfa['v[3]'],
            pfa['v[4]'],
            pfa['v[5]'],
            pfa['vlogv[5]'],
            pfa['v[6]'],
            pfa['v[6]'],
            pfa['vlogv[6]'],
            pfa['v[7]']

        SimInspiralTaylorF2AlignedPhasing
        This code is taken from the top LALSimIMRPhenomD_internals.c
        It implements the, then currrent, state of the PN calculations in
        LAL. After PhenomD was constructed the 3PN spin-spin term
        was added to LAL. This implementation does NOT use this term.
        This implementation + 3PN spin-spin term should equate to LAL (07/06/2016)

         - original text from LALSimIMRPhenomD_internals.c:
        ''This waveform uses the TaylorF2 coefficients for it's inspiral phase augmented
        by higher order phenomenological terms tuned to SEOBv2-Hybrid waveforms.
        Below are lines copied from LALSimInspiralPNCoefficients.c which are the TaylorF2
        phase coefficients we have used.
        We document them here in case changes to that file changes the behaviour
        of this waveform.''
        """
        # define dictionary to store PN coefficients
        from collections import OrderedDict
        # pfa = {}
        pfa = OrderedDict()

        m1 = p['m1']
        m2 = p['m2']
        chi1L = p['chi1z']
        chi2L = p['chi2z']
        mtot = p['Mtot']
        eta = p['eta']
        d = (m1 - m2) / (m1 + m2)
        m1M = m1/mtot
        m2M = m2/mtot

        chi1sq = chi1L * chi1L
        chi2sq = chi2L * chi2L
        chi1dotchi2 = chi1L * chi2L

        LAL_PI = pi
        LAL_GAMMA = 0.5772156649015329
        #hardcode tidal deformation to correspond to BHs i.e. unity
        qm_def1 = 1.
        qm_def2 = 1.

        # // Use the spin-orbit variables from arXiv:1303.7412, Eq. 3.9
        # // We write dSigmaL for their (\delta m/m) * \Sigma_\ell
        # // There's a division by mtotal^2 in both the energy and flux terms
        # // We just absorb the division by mtotal^2 into SL and dSigmaL
        SL = m1M * m1M * chi1L + m2M * m2M * chi2L
        dSigmaL = d * (m2M * chi2L - m1M * chi1L)
        pfaN = 3/(128 * eta)
        # //Non-spin phasing terms - see arXiv:0907.0700, Eq. 3.18
        pfa['v[0]'] = 1
        pfa['v[1]'] = 0
        pfa['v[2]'] = 5*(743/84 + 11 * eta)/9
        pfa['v[3]'] = -16*LAL_PI
        pfa['v[4]'] = 5*(3058.673/7.056 + 5429/7 * eta \
                         + 617 * eta*eta)/72
        pfa['v[5]'] = 5/9 * (7729/84 - 13 * eta) * LAL_PI
        pfa['vlogv[5]'] = 5/3 * (7729/84 - 13 * eta) * LAL_PI
        pfa['v[6]'] = (11583.231236531/4.694215680 \
                         - 640/3 * LAL_PI * LAL_PI - 6848/21*LAL_GAMMA) \
                     + eta * (-15737.765635/3.048192 \
                         + 2255./12. * LAL_PI * LAL_PI) \
                     + eta*eta * 76055/1728 \
                     - eta*eta*eta * 127825/1296
        pfa['v[6]'] += (-6848/21)*log(4.)
        pfa['vlogv[6]'] = -6848/21
        pfa['v[7]'] = LAL_PI * ( 77096675/254016 \
                         + 378515/1512 * eta - 74045/756 * eta*eta)


        # /* Compute 2.0PN SS, QM, and self-spin */
        # // See Eq. (6.24) in arXiv:0810.5336
        # // 9b,c,d in arXiv:astro-ph/0504538
        pn_sigma = eta * (721/48*chi1L*chi2L - 247/48*chi1dotchi2)
        pn_sigma += (720*qm_def1 - 1)/96.0 * m1M * m1M * chi1L * chi1L
        pn_sigma += (720*qm_def2 - 1)/96.0 * m2M * m2M * chi2L * chi2L
        pn_sigma -= (240*qm_def1 - 7)/96.0 * m1M * m1M * chi1sq
        pn_sigma -= (240*qm_def2 - 7)/96.0 * m2M * m2M * chi2sq

        # 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
        # was not available when PhenomD was tuned.
        # so is not used
        # pn_ss3 =  (326.75/1.12 + 557.5/1.8*eta)*eta*chi1L*chi2L
        # pn_ss3 += ((4703.5/8.4+2935/6.*m1M-120.*m1M*m1M)*qm_def1 + (-4108.25/6.72-108.5/1.2*m1M+125.5/3.6*m1M*m1M)) *m1M*m1M * chi1sq
        # pn_ss3 += ((4703.5/8.4+2935./6.*m2M-120.*m2M*m2M)*qm_def2 + (-4108.25/6.72-108.5/1.2*m2M+125.5/3.6*m2M*m2M)) *m2M*m2M * chi2sq

        # // Spin-orbit terms - can be derived from arXiv:1303.7412, Eq. 3.15-16
        pn_gamma = (554345/1134 + 110*eta/9)*SL + (13915/84 - 10*eta/3.)*dSigmaL
        # case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN
        pfa['v[7]'] += (-8980424995/762048 + 6586595*eta/756 \
        - 305*eta*eta/36)*SL - (170978035/48384 - 2876425*eta/672 - 4735*eta*eta/144) * dSigmaL
        # case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        # 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
        # was not available when PhenomD was tuned.
        # so is not used.
        # pfa['v[6]'] += LAL_PI * (3760*SL + 1490*dSigmaL)/3 + pn_ss3
        pfa['v[6]'] += LAL_PI * (3760*SL + 1490*dSigmaL)/3
        # case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        pfa['v[5]'] += -1 * pn_gamma
        pfa['vlogv[5]'] += -3 * pn_gamma
        # case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
        pfa['v[4]'] += -10 * pn_sigma
        # case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
        pfa['v[3]'] += 188*SL/3 + 25*dSigmaL

        # finally we multiply every term by the Newtonian term
        for key in pfa.keys():
            pfa[key] *= pfaN

        return pfa

    def init_phi_ins_prefactors(self, model_pars, powers_of_pi):
        """helper function for phase
        p : dict
        pn : dict of PN coefficients from SimInspiralTaylorF2AlignedPhasing
        powers_of_pi : instant of UsefulPowers class
        output:
        output a dictionary called prefactors"""

        pn = model_pars['pn']

        prefactors = {}

        sigma1 = model_pars['sigma1']
        sigma2 = model_pars['sigma2']
        sigma3 = model_pars['sigma3']
        sigma4 = model_pars['sigma4']
        Pi = pi
        LAL_PI_4 = pi / 4.

        # // PN phasing series
        prefactors['initial_phasing'] = pn['v[5]'] - LAL_PI_4
        prefactors['two_thirds'] = pn['v[7]'] * powers_of_pi.two_thirds
        prefactors['third'] = pn['v[6]'] * powers_of_pi.third
        prefactors['third_with_logv'] = pn['vlogv[6]'] * powers_of_pi.third
        prefactors['logv'] = pn['vlogv[5]']
        prefactors['minus_third'] = pn['v[4]'] / powers_of_pi.third
        prefactors['minus_two_thirds'] = pn['v[3]'] / powers_of_pi.two_thirds
        prefactors['minus_one'] = pn['v[2]'] / Pi
        prefactors['minus_five_thirds'] = pn['v[0]'] / powers_of_pi.five_thirds #// * v^0

        # // higher order terms that were calibrated for PhenomD
        prefactors['one'] = sigma1
        prefactors['four_thirds'] = sigma2 * 3.0/4.0
        prefactors['five_thirds'] = sigma3 * 3.0/5.0
        prefactors['two'] = sigma4 / 2.0

        return prefactors

    def PhiInsAnsatzInt(self, Mf, eta, prefactors, powers_of_Mf, powers_of_pi):
        """
        input
        Mf : dimensionless frequency
        powers_of_Mf : instance of UsefulPowers
        prefactors : dict
                output from init_phi_ins_prefactors function
        p : dict
        Ansatz for the inspiral phase.
        We call the LAL TF2 coefficients here.
        The exact values of the coefficients used are given
        as comments in the top of this file
        Defined by Equation 27 and 28 arXiv:1508.07253
        """

        # // Assemble PN phasing series
        v = powers_of_Mf.third * powers_of_pi.third
        logv = log(v)

        phasing = prefactors['initial_phasing']

        phasing += prefactors['two_thirds']	* powers_of_Mf.two_thirds
        phasing += prefactors['third'] * powers_of_Mf.third
        phasing += prefactors['third_with_logv'] * logv * powers_of_Mf.third
        phasing += prefactors['logv'] * logv
        phasing += prefactors['minus_third'] / powers_of_Mf.third
        phasing += prefactors['minus_two_thirds'] / powers_of_Mf.two_thirds
        phasing += prefactors['minus_one'] / Mf
        phasing += prefactors['minus_five_thirds'] / powers_of_Mf.five_thirds #// * v^0

        # // Now add higher order terms that were calibrated for PhenomD
        phasing += ( prefactors['one'] * Mf + prefactors['four_thirds'] * powers_of_Mf.four_thirds \
        		   + prefactors['five_thirds'] * powers_of_Mf.five_thirds \
        		   + prefactors['two'] * powers_of_Mf.two \
        		 ) / eta

        return phasing


    def DPhiInsAnsatzInt(self, Mf, powers_of_pi, model_pars, eta):
        """
        input
        Mf : dimensionless frequency
        powers_of_pi : instance of UsefulPowers
        prefactors : dict
                output from init_phi_ins_prefactors function
        p : dict
        pn : dict
        First frequency derivative of PhiInsAnsatzInt
        """

        pn = model_pars['pn']

        sigma1 = model_pars['sigma1']
        sigma2 = model_pars['sigma2']
        sigma3 = model_pars['sigma3']
        sigma4 = model_pars['sigma4']
        Pi = pi

        # // Assemble PN phasing series
        v = (Pi*Mf) ** (1./3.)
        logv = log(v)
        v2 = v * v
        v3 = v * v2
        v4 = v * v3
        v5 = v * v4
        v6 = v * v5
        v7 = v * v6
        v8 = v * v7

        # // Apply the correct prefactors to LAL phase coefficients to get the
        # // phase derivative dphi / dMf = dphi/dv * dv/dMf
        Dphasing = 0.0
        Dphasing += +2.0 * pn['v[7]'] * v7
        Dphasing += (pn['v[6]'] + pn['vlogv[6]'] * (1.0 + logv)) * v6
        Dphasing += pn['vlogv[5]'] * v5
        Dphasing += -1.0 * pn['v[4]'] * v4
        Dphasing += -2.0 * pn['v[3]'] * v3
        Dphasing += -3.0 * pn['v[2]'] * v2
        Dphasing += -4.0 * pn['v[1]'] * v
        Dphasing += -5.0 * pn['v[0]']
        Dphasing /= v8 * 3.0/Pi

        # // Now add higher order terms that were calibrated for PhenomD
        Dphasing += ( \
          sigma1 \
        + sigma2 * v / powers_of_pi.third \
        + sigma3 * v2 / powers_of_pi.two_thirds \
        + (sigma4/Pi) * v3 \
        ) / eta

        return Dphasing

    def ComputeIMRPhenDPhaseConnectionCoefficients(self, eta, model_pars, powers_of_pi, rholm=1.0, taulm=1.0):
        """
        p : dict

        returns
            C1Int, C2Int, C1MRD, C2MRD

        This function aligns the three phase parts (inspiral, intermediate and merger-rindown)
        such that they are c^1 continuous at the transition frequencies
        Defined in VIII. Full IMR Waveforms arXiv:1508.07253
        Rholm was added when IMRPhenomHM (high mode) was added.
        Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
        to 1. and PhenomD is recovered.
        Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
        Again, when Taulm = 1.0 then PhenomD is recovered.
        """

        # // Transition frequencies
        # // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
        fInsJoin=self.PHI_fJoin_INS
        fMRDJoin=0.5*model_pars['fRD']

        prefactors = model_pars['phi_prefactors']

        # // Compute C1Int and C2Int coeffs
        # // Equations to solve for to get C(1) continuous join
        # // PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
        # // Joining at fInsJoin
        # // PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
        # // PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int

        DPhiIns = self.DPhiInsAnsatzInt(fInsJoin, powers_of_pi, model_pars, eta)
        DPhiInt = self.DPhiIntAnsatz(fInsJoin, model_pars, eta)
        C2Int = DPhiIns - DPhiInt

        powers_of_fInsJoin = UsefulPowers(fInsJoin)
        C1Int = self.PhiInsAnsatzInt(fInsJoin, eta, prefactors, powers_of_fInsJoin, powers_of_pi) \
        - 1.0/eta * self.PhiIntAnsatz(fInsJoin, model_pars) - C2Int * fInsJoin

        # // Compute C1MRD and C2MRD coeffs
        # // Equations to solve for to get C(1) continuous join
        # // PhiInsInt (f)  =   PhiMRD (f) + C1MRD + C2MRD f
        # // Joining at fMRDJoin
        # // Where \[Phi]InsInt(f) is the \[Phi]Ins+\[Phi]Int joined function
        # // PhiInsInt (fMRDJoin)  =   PhiMRD (fMRDJoin) + C1MRD + C2MRD fMRDJoin
        # // PhiInsInt'(fMRDJoin)  =   PhiMRD'(fMRDJoin) + C2MRD
        # // temporary Intermediate Phase function to Join up the Merger-Ringdown
        PhiIntTempVal = 1.0/eta * self.PhiIntAnsatz(fMRDJoin, model_pars) + C1Int + C2Int*fMRDJoin
        DPhiIntTempVal = self.DPhiIntTemp(fMRDJoin, model_pars, eta, C2Int)
        DPhiMRDVal = self.DPhiMRD(fMRDJoin, model_pars, eta, rholm, taulm)
        C2MRD = DPhiIntTempVal - DPhiMRDVal
        C1MRD = PhiIntTempVal - 1.0/eta * self.PhiMRDAnsatzInt(fMRDJoin, model_pars, rholm, taulm) - C2MRD*fMRDJoin
        return C1Int, C2Int, C1MRD, C2MRD

class PhenomDInternals(PhenomDInternalsAmplitude, PhenomDInternalsPhase):
    """docstring for PhenomDInternals"""
    def __init__(self):
        pass

class PhenomD(PhenomDInternals):
    """docstring for PhenomD"""

    #sets phenomD CONSTANTS
    #https://github.com/lscsoft/lalsuite/blob/master/lalsimulation/src/LALSimIMRPhenomD.h

    # /* CONSTANTS */

    # /**
    #  * Dimensionless frequency (Mf) at which define the end of the waveform
    #  */
    Mf_CUT = 0.2

    # /**
    #  * Dimensionless frequency (Mf) at which the inspiral amplitude
    #  * switches to the intermediate amplitude
    #  */
    AMP_fJoin_INS = 0.014

    # /**
    #  * Dimensionless frequency (Mf) at which the inspiral phase
    #  * switches to the intermediate phase
    #  */
    PHI_fJoin_INS = 0.018

    # /**
    #   * Minimal final spin value below which the waveform might behave pathological
    #   * because the ISCO frequency is too low. For more details, see the review wiki
    #   * page https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/PhenD_LargeNegativeSpins
    #   */
    #NOTE: This is not implemented in this python version
    MIN_FINAL_SPIN = -0.717

    # /**
    #   * A large mass ratio causes memory over-runs.
    #   * We test and put the limit an order of magnitude above that of previous waveform models (which were around q=100).
    #   */
    #NOTE: This is not implemented in this python version
    MAX_ALLOWED_MASS_RATIO = 5000

    def __init__(self,
                 m1=10.,
                 m2=10.,
                 chi1z=0.,
                 chi2z=0.,
                 f_min=20.,
                 f_max=0.,
                 delta_f=1/64.,
                 distance=1e6 * Constants.PC_SI,
                 fRef=0.,
                 phiRef=0.,
                 finspin_func="FinalSpin0815",
                 rholm=1.0,
                 taulm=1.0,
                 **kwargs):
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
        finspin_func='FinalSpin0815' (used because PhenomP changes what is used to calculate the final spin.)
        Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal to 1. and PhenomD is recovered.
        Taulm = fDMlm/fDM22. Ratio of ringdown damping times. Again, when Taulm = 1.0 then PhenomD is recovered.
        **kwargs ( currently the only extra kwarg that is passed is chip or chi1x. )

        Rholm and Taulm were added when IMRPhenomHM (high mode) was added.
        """

        #TODO: Refactor this PhenomD code.
        #so that you can call something like
        #self.generate_phenomD_model_parameters()
        #and it will populate self.model_parameters()
        #containing things like alpha1, alpha2, etc...

        # enforce m1 >= m2 and chi1 is on m1
        if m1<m2: # swap spins and masses
            chi1z, chi2z = chi2z, chi1z
            m1, m2 = m2, m1

        self.p                = {}
        self.p['m1']          = float(m1)
        self.p['m2']          = float(m2)
        self.p['chi1z']       = float(chi1z)
        self.p['chi2z']       = float(chi2z)
        self.p['f_min']       = float(f_min)
        self.p['f_max']       = float(f_max)
        self.p['delta_f']     = float(delta_f)
        self.p['distance']    = float(distance)
        self.p['fRef']        = float(fRef)
        self.p['phiRef']      = float(phiRef)

        #TODO: This will eventually be chip when I write the
        #'transform to model parameters' function
        if "chip" in kwargs:
            self.p['chip'] = float(kwargs['chip'])

        self.p['Mtot'], self.p['eta'] = M_eta_m1_m2(self.p['m1'], self.p['m2'])
        self.p['chipn'] = chipn(self.p['eta'], self.p['chi1z'], self.p['chi2z'])

        #initialise UsefulPowers instances
        self.powers_of_pi = UsefulPowers(pi)

        #model_pars is a dict containing all? the internal variables of phenomD
        #final spin is a variable because phenomP uses a different final spin to
        #phenomD
        self.finspin_func = finspin_func
        self.model_pars = self.compute_model_parameters(self.p, self.finspin_func, rholm, taulm)


        #end of __init__()

    def compute_model_parameters(self, p, finspin_func, rholm=1.0, taulm=1.0):
        """
        compute_model_parameters(p, finspin_func)
        p : dict
        finspin_func : string, name of final spin function. Used in phenomP as it
        uses a different final spin function.
        Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal to 1. and PhenomD is recovered.
        Taulm = fDMlm/fDM22. Ratio of ringdown damping times. Again, when Taulm = 1.0 then PhenomD is recovered.
        Rholm and Taulm were added when IMRPhenomHM (high mode) was added.
        """

        model_pars = {}

        # Compute the amplitude pre-factor
        model_pars['amp0'] = 2. * sqrt(5. / (64.*pi)) * p['Mtot'] * Constants.MRSUN_SI * p['Mtot'] * Constants.MTSUN_SI / p['distance']

        model_pars['M_sec'] = p['Mtot'] * Constants.MTSUN_SI # Conversion factor Hz -> dimensionless frequency

        # Default PhenomD values
        if p['f_max'] == 0. : p['f_max'] = self.Mf_CUT / model_pars['M_sec'] # converted from Mf to Hz
        if p['fRef'] == 0.  : p['fRef'] = p['f_min']

        # compute prediction of ringdown frequency

        if finspin_func == "FinalSpin0815":
            model_pars['finspin'] = FinalSpin0815(p['eta'], p['chi1z'], p['chi2z'])
        elif finspin_func == "FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH":
            model_pars['finspin'] = FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH(p['m1'], p['m2'], p['chip'], p['chi1z'], p['chi2z'])
            if( absolute(model_pars['finspin']) > 1.0 ):
                print("Warning: final spin magnitude {0} > 1. Setting final spin magnitude = 1.".format(model_pars['finspin']))
                model_pars['finspin'] = copysign(1.0, model_pars['finspin'])
        else:
            raise ValueError('finspin_func = {0} not recognised'.format(finspin_func))

        model_pars['fRD']     = fring(p['eta'], p['chi1z'], p['chi2z'], model_pars['finspin'])
        model_pars['fDM']     = fdamp(p['eta'], p['chi1z'], p['chi2z'], model_pars['finspin'])

        #inspiral amplitude coeffs
        model_pars['rho1'] = self.rho1_fun(p)
        model_pars['rho2'] = self.rho2_fun(p)
        model_pars['rho3'] = self.rho3_fun(p)

        #inspiral amplitude prefactors
        #only need to pass rho1, rho2 and rho3 from model_pars
        model_pars['amp_prefactors'] = super(PhenomDInternals, self).init_amp_ins_prefactors(p, model_pars, self.powers_of_pi)

        #merger-ringdown amplitude coeffs
        model_pars['gamma1'] = self.gamma1_fun(p)
        model_pars['gamma2'] = self.gamma2_fun(p)
        model_pars['gamma3'] = self.gamma3_fun(p)

        # fmaxCalc is the peak of the merger-ringdown amplitude function
        model_pars['fmaxCalc'] = self.fmaxCalc(model_pars)

        #intermediate amplitude coeffs
        model_pars.update(self.ComputeDeltasFromCollocation(p, model_pars, self.powers_of_pi))

        # populate phase phenom coefficient dictionaries
        #inspiral PN coefficients
        model_pars['pn'] = self.SimInspiralTaylorF2AlignedPhasing(p)

        #inspiral phase coeffs
        model_pars['sigma1'] = self.sigma1Fit(p)
        model_pars['sigma2'] = self.sigma2Fit(p)
        model_pars['sigma3'] = self.sigma3Fit(p)
        model_pars['sigma4'] = self.sigma4Fit(p)

        #inspiral phase prefactors
        model_pars['phi_prefactors'] = self.init_phi_ins_prefactors(model_pars, self.powers_of_pi)

        #intermediate phase coeffs
        model_pars['beta1'] = self.beta1Fit(p)
        model_pars['beta2'] = self.beta2Fit(p)
        model_pars['beta3'] = self.beta3Fit(p)
        #merger-ringdown phase coeffs
        model_pars['alpha1'] = self.alpha1Fit(p)
        model_pars['alpha2'] = self.alpha2Fit(p)
        model_pars['alpha3'] = self.alpha3Fit(p)
        model_pars['alpha4'] = self.alpha4Fit(p)
        model_pars['alpha5'] = self.alpha5Fit(p)

        #compute phase connection coefficients
        model_pars['C1Int'], model_pars['C2Int'], model_pars['C1MRD'], model_pars['C2MRD'] = self.ComputeIMRPhenDPhaseConnectionCoefficients(p['eta'], model_pars, self.powers_of_pi, rholm, taulm)

        # compute time shift such that the ifft of htilde has it's peak
        # amplitude (near) t=0. (only approximate)
        model_pars['t0'] = self.computetimeshift(model_pars, p['eta'])

        # // incorporating fRef
        model_pars['MfRef'] = model_pars['M_sec'] * p['fRef']
        powers_of_fRef = UsefulPowers(model_pars['MfRef'])
        model_pars['phi_fRef'] = self.IMRPhenomDPhase(model_pars['MfRef'], p['eta'], model_pars, powers_of_fRef, rholm, taulm)
        # // factor of 2 b/c phi0 is orbital phase
        model_pars['phi_precalc'] = 2.*p['phiRef'] + model_pars['phi_fRef']

        return model_pars

    def IMRPhenomDAmplitude(self, Mf, model_pars, powers_of_Mf):
        """Call ComputeIMRPhenomDAmplitudeCoefficients() first!
        This function computes the IMR amplitude given phenom coefficients.
        Defined in VIII. Full IMR Waveforms arXiv:1508.07253
        """
        # Defined in VIII. Full IMR Waveforms arXiv:1508.07253
        # The inspiral, intermediate and merger-ringdown amplitude parts

        # Transition frequencies
        fInsJoin = self.AMP_fJoin_INS
        fMRDJoin = model_pars['fmaxCalc']

        amp_prefactors = model_pars['amp_prefactors']

        AmpPreFac = amp_prefactors['PNamp0'] / powers_of_Mf.seven_sixths
        # split the calculation to just 1 of 3 possible mutually exclusive ranges
        # this effectively implements a step function transition function
        if (Mf <= fInsJoin):	# Inspiral range
            AmpIns = AmpPreFac * self.AmpInsAnsatz(Mf, powers_of_Mf, amp_prefactors)
            return AmpIns
        elif (Mf >= fMRDJoin):	# MRD range
            AmpMRD = AmpPreFac * self.AmpMRDAnsatz(Mf, model_pars)
            return AmpMRD
        else:
            #	Intermediate range
            AmpInt = AmpPreFac * self.AmpIntAnsatz(Mf, model_pars)
            return AmpInt

    def IMRPhenomDPhase(self, Mf, eta, model_pars, powers_of_Mf, rholm=1.0, taulm=1.0 ):
        """
        This function computes the IMR phase given phenom coefficients.
        Defined in VIII. Full IMR Waveforms arXiv:1508.07253
        Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal to 1. and PhenomD is recovered.
        Taulm = fDMlm/fDM22. Ratio of ringdown damping times. Again, when Taulm = 1.0 then PhenomD is recovered.
        Rholm and Taulm were added when IMRPhenomHM (high mode) was added.
        """
        # Defined in VIII. Full IMR Waveforms arXiv:1508.07253
        # The inspiral, intermendiate and merger-ringdown phase parts
        pn = model_pars['pn']

        # split the calculation to just 1 of 3 possible mutually exclusive ranges
        fInsJoin = self.PHI_fJoin_INS
        fMRDJoin = 0.5*model_pars['fRD']
        phi_prefactors = model_pars['phi_prefactors']

        if (Mf <= fInsJoin):	# Inspiral rangee
            PhiIns = self.PhiInsAnsatzInt(Mf, eta, phi_prefactors, powers_of_Mf, self.powers_of_pi)
            return PhiIns
        elif (Mf >= fMRDJoin):	# MRD range
            PhiMRD = 1.0/eta * self.PhiMRDAnsatzInt(Mf, model_pars, rholm, taulm) + model_pars['C1MRD'] + model_pars['C2MRD'] * Mf
            return PhiMRD
        else:
            #	Intermediate range
            PhiInt = 1.0/eta * self.PhiIntAnsatz(Mf, model_pars) + model_pars['C1Int'] + model_pars['C2Int'] * Mf
            return PhiInt

    def IMRPhenomDPhaseDerivFrequencyPoint(self, Mf, eta, model_pars, rholm=1.0, taulm=1.0):
        """
        ONLY TESTED FOR (l,m)=(2,2) MODE

        Implementing https://git.ligo.org/lscsoft/lalsuite/blob/master/lalsimulation/src/LALSimIMRPhenomD.c#L488

        Returns the value of the frequency derivative of the phase.
        Used to compute the chirp time of a phenomD waveform.
        """
        pn = model_pars['pn']

        # split the calculation to just 1 of 3 possible mutually exclusive ranges
        fInsJoin = self.PHI_fJoin_INS
        fMRDJoin = 0.5*model_pars['fRD']
        phi_prefactors = model_pars['phi_prefactors']

        if (Mf <= fInsJoin):	# Inspiral rangee
            DPhiIns = self.DPhiInsAnsatzInt(Mf, self.powers_of_pi, model_pars, eta)
            return DPhiIns
        elif (Mf >= fMRDJoin):	# MRD range
            DPhiMRDval = self.DPhiMRD(Mf, model_pars, eta, rholm, taulm) + model_pars['C2MRD']
            return DPhiMRDval
        else:
            #	Intermediate range
            DPhiInt = self.DPhiIntAnsatz(Mf, model_pars, eta) + model_pars['C2Int']
            return DPhiInt

    def computetimeshift(self, model_pars, eta):
        # //time shift so that peak amplitude is approximately at t=0
        # //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
        t0 = self.DPhiMRD(model_pars['fmaxCalc'], model_pars, eta)
        return t0

    def IMRPhenomDGenerateFDSingleFrequencyPoint(self, Mf):
        """
        standalone phenomD strain generator for a single frequency point
        input:
            Mf : geometric frequency
        output:
            htilde at Mf
        """
        powers_of_Mf = UsefulPowers(Mf)
        amp   = self.IMRPhenomDAmplitude(Mf, self.model_pars, powers_of_Mf)
        phase = self.IMRPhenomDPhase(Mf, self.p['eta'], self.model_pars, powers_of_Mf, rholm=1.0, taulm=1.0)
        # finally supply time and phase shift for fRef and approximately t=0 at amplitude max in time domain.
        phase -= self.model_pars['t0']*(Mf-self.model_pars['MfRef']) + self.model_pars['phi_precalc']
        #amp0 scales to physical distance and correct units for strain
        htilde = self.model_pars['amp0'] * amp * exp(-1.j * phase)
        return htilde

    def IMRPhenomDGenerateFD(self):
        """
        generates frequency domain strain
        calls IMRPhenomDGenerateFDSingleFrequencyPoint() in a loop.
        attributes :
            flist_Hz
            htilde
        """
        self.flist_Hz = arange(self.p['f_min'], self.p['f_max'], self.p['delta_f'])
        self.htilde = zeros(len(self.flist_Hz), dtype=complex)
        for i in range(len(self.flist_Hz)):
            #waveform generation is in terms of Mf so need to convert from Hz to Mf
            Mf = self.flist_Hz[i] * self.model_pars['M_sec']
            self.htilde[i] = self.IMRPhenomDGenerateFDSingleFrequencyPoint(Mf)

    def getampandphase(self, htilde):
        """
        input: htilde
        attributes :
            amp
            phase

        given IMRPhenomDGenerateFD is run and self.htilde is initialised,
        compute the amplitude and phase of the waveform.
        """
        self.amp = absolute(htilde)
        self.phase = unwrap(angle(htilde))

    def ChirpTime(self, fStart, units='s'):
        """
        input:
            fStart (start frequency in Hz)
            units = 's' or 'M' either seconds or in geometic units of the total Mass

        output: the chirp time in seconds
        Return the chirp time in seconds or geometic units for this system
        """

        Mf_PeakFreq = self.model_pars['fmaxCalc']
        Hz_PeakFreq = Mf_PeakFreq / self.model_pars['M_sec']

        Hz_Start = fStart
        Mf_Start = HztoMf(fStart, self.p['Mtot'])

        if (Hz_Start > Hz_PeakFreq):
            print("Error: Peak frequency is greater than start frequency")
            sys.exit(1)

        dphiSt = self.IMRPhenomDPhaseDerivFrequencyPoint(Mf_Start, self.p['eta'], self.model_pars, rholm=1.0, taulm=1.0)
        dphifRD = self.IMRPhenomDPhaseDerivFrequencyPoint(Mf_PeakFreq, self.p['eta'], self.model_pars, rholm=1.0, taulm=1.0)

        dphidiff = dphifRD - dphiSt

        ChirpTimeSec = dphidiff / 2. / Constants.LAL_PI * self.model_pars['M_sec']

        if units == 's':
            ChirpTimeSec = ChirpTimeSec
        elif units == 'M':
            ChirpTimeSec = StoM(ChirpTimeSec, self.p['Mtot'])
        else:
            print("Units must be either 's' or 'M'")
            sys.exit(1)

        return ChirpTimeSec