#TODO: do me. swsh.py = spin-weighted spherical harmonics

from numpy import absolute, sqrt, cos, sin, exp
from phenom.utils.utils import Constants

class SWSH(object):
    """docstring for SWSH = spin-weighted spherical harmonics"""
    def __init__(self, ell=2, m=2, theta=0., phi=0.):
        """
        ell [2]   : int
            mode number ell
        m   [2]   : int
            mode number m
        theta [0] : float
            polar angle (rad)
        phi  [0]  : float
            azimuthal angle (rad)
        From: lal/src/utilities/SphericalHarmonics.c
        """

        self.ell = int(ell)
        self.m = int(m)
        self.theta = float(theta)
        self.phi = float(phi)

        if self.ell < absolute(self.m):
            raise ValueError('ell = {0} < |m|={1}'.format(self.ell, absolute(self.m)))

        self.val = self.getYlm(self.ell, self.m, self.theta, self.phi)


    def getYlm(self, ell, m, theta, phi):
        LAL_PI = Constants.LAL_PI
        if ell == 2:
            if m == -2:
                ret = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta )) * self.cpolar(m, phi)
            elif m == -1:
                ret = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 - cos( theta )) * self.cpolar(m, phi)
            elif m == 0:
                ret = sqrt( 15.0 / ( 32.0 * LAL_PI ) ) * sin( theta )*sin( theta )
            elif m == 1:
                ret = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 + cos( theta )) * self.cpolar(m, phi)
            elif m == 2:
                ret = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta )) * self.cpolar(m, phi)
            else:
                raise ValueError('m = {0} not valid, max absolute value is {1}'.format(m, ell))
            return ret
        else:
            raise ValueError('input ell = {0} not valid. Only ell = 2 supported'.format(ell))

    def cpolar(self, m, phi):
        """
        m     : int
            mode number m
        phi   : float
            azimuthal angle (rad)
        computes exp(1.j * m * phi)"""
        return exp(1.j * m * phi)
