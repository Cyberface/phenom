from __future__ import division

from phenom.utils.utils import M_eta_m1_m2, chipn, Constants, chieffPH

class PhenomP():
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


        self.M_sec = self.p['Mtot'] * Constants.MTSUN_SI # Conversion factor Hz -> dimensionless frequency

        self.piM = Constants.LAL_PI * self.M_sec

        #chipn is used in the aligned model 'phenomD'
        self.p['chipn']  = chipn(self.p['eta'], self.p['chi1z'], self.p['chi2z'])
        #chieff is used in the twisting part
        self.p['chieff'] = chieffPH(self.p['m1'], self.p['m2'], self.p['chi1z'], self.p['chi2z'])



        # ph = phenom.PhenomD(**self.p)
        # ph.IMRPhenomDGenerateFD()
        # ph.getampandphase(ph.htilde)
        #this means we now have ph.amp and ph.phase.










        #
