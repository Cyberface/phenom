from phenom.waveform.waveform import Waveform
from phenom.utils.utils import M_eta_m1_m2, chipn, UsefulPowers
from numpy import sqrt, pi

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

    def init_amp_ins_prefactors(self, p, amp_dict, powers_of_pi):
        """helper function for AmpInsAnsatz
        p : dict
        amp_dict : dict
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
        rho1 = amp_dict['rho1']
        rho2 = amp_dict['rho2']
        rho3 = amp_dict['rho3']

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

    def AmpInsAnsatz(Mf, powers_of_Mf, prefactors):
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


class PhenomDInternals(PhenomDInternalsAmplitude):
    """docstring for PhenomDInternals"""
    def __init__(self):
        pass

# TODO:0 Figure out how to structure PhenomD class
class PhenomD(Waveform, PhenomDInternals):
    """docstring for PhenomD"""
    # Use the constructor from the Waveform Class
    def __init__(self, **p):
        super(PhenomD, self).__init__(p)

        #initialise UsefulPowers instances
        self.powers_of_pi = UsefulPowers(pi)
        # self.powers_of_Mf = UsefulPowers(Mf)

        self.p['Mtot'], self.p['eta'] = M_eta_m1_m2(self.p['m1'], self.p['m2'])
        self.p['chipn'] = chipn(self.p['eta'], self.p['chi1z'], self.p['chi2z'])
        print self.rho3_fun(self.p)

        #TODO: populate amplitude and phase phenom coefficient dictionaries
        self.amp_dict = {}
        self.amp_dict['rho1'] = self.rho1_fun(self.p)
        self.amp_dict['rho2'] = self.rho2_fun(self.p)
        self.amp_dict['rho3'] = self.rho3_fun(self.p)

        self.amp_prefactors = super(PhenomDInternals, self).init_amp_ins_prefactors(self.p, self.amp_dict, self.powers_of_pi)

        # print super(PhenomD, self).rho3_fun(self.p['eta'], self.p['chipn'])
        # print super(PhenomDInternals, self).rho3_fun(self.p['eta'], self.p['chipn'])
        # super(PhenomD, self).interal()
        # super(PhenomDInternals, self).interal()

    def amp(self):
        pass


    def h(self):
        """
        generates frequency domain strain"""
        pass
