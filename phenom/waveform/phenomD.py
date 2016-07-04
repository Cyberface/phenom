from phenom.waveform.waveform import Waveform
from phenom.utils.utils import M_eta_m1_m2, chipn
from numpy import sqrt, pi

class UsefulPowers(object):
    """init_useful_powers from phenomD LAL code"""
    def __init__(self, number):
        self.sixth        = number ** (1.0/6.0)
        self.third        = self.sixth * self.sixth
        self.two_thirds   = number / self.third
        self.four_thirds  = number * (self.third);
        self.five_thirds  = self.four_thirds * (self.third);
        self.two          = number * number;
        self.seven_thirds = self.third * self.two;
        self.eight_thirds = self.two_thirds * self.two;

#initialise powers_of_pi
powers_of_pi = UsefulPowers(pi)

class PhenomDInternals(object):
    """docstring for PhenomDInternals"""
    def __init__(self, m1, m2, chi1z, chi2z):
        self.m1 = m1
        self.m2 = m2
        self.chi1z = chi1z
        self.chi2z = chi2z
        self.M, self.eta = M_eta_m1_m2(self.m1, self.m2)
        self.chi = chipn(self.eta, self.chi1z, self.chi2z)
        # self.arg = arg

class PhenomDInternalsAmplitude(PhenomDInternals):
    """docstring for PhenomDInternalsAmplitude"""
    def __init__(self, arg):
        self.arg = arg

    # ///////////////////////////// Amplitude: Inspiral functions /////////////////////////

    # // Phenom coefficients rho1, ..., rho3 from direct fit
    # // AmpInsDFFitCoeffChiPNFunc[eta, chiPN]

    def rho1_fun(self, eta, chi):
        """
        rho_1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 3931.8979897196696 - 17395.758706812805*eta \
        + (3132.375545898835 + 343965.86092361377*eta - 1.2162565819981997e6*eta2)*xi \
        + (-70698.00600428853 + 1.383907177859705e6*eta - 3.9662761890979446e6*eta2)*xi2 \
        + (-60017.52423652596 + 803515.1181825735*eta - 2.091710365941658e6*eta2)*xi3

    def rho2_fun(self, eta, chi):
        """
        rho_2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return -40105.47653771657 + 112253.0169706701*eta \
        + (23561.696065836168 - 3.476180699403351e6*eta + 1.137593670849482e7*eta2)*xi \
        + (754313.1127166454 - 1.308476044625268e7*eta + 3.6444584853928134e7*eta2)*xi2 \
        + (596226.612472288 - 7.4277901143564405e6*eta + 1.8928977514040343e7*eta2)*xi3

    def rho3_fun(self, eta, chi):
        """
        rho_3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
        """
        xi = -1 + chi
        xi2 = xi*xi
        xi3 = xi2*xi
        eta2 = eta*eta

        return 83208.35471266537 - 191237.7264145924*eta + \
        (-210916.2454782992 + 8.71797508352568e6*eta - 2.6914942420669552e7*eta2)*xi \
        + (-1.9889806527362722e6 + 3.0888029960154563e7*eta - 8.390870279256162e7*eta2)*xi2 \
        + (-1.4535031953446497e6 + 1.7063528990822166e7*eta - 4.2748659731120914e7*eta2)*xi3

class PhenomD(PhenomDInternals, Waveform):
    """docstring for PhenomD"""
    def __init__(self, **p):
        self.p = p
    def amp(self):
        pass
