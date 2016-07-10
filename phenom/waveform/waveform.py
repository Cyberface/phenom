
from phenom.utils.utils import M_eta_m1_m2, chipn, Constants

class Waveform(object):
    """docstring for Waveform
    p : dict
        dictionary of physical parameters:
        'mass1', 'mass2',
        'spin1x', 'spin1y', 'spin1z',
        'spin2x', 'spin2y', 'spin2z',
        'f_min', ,'f_max' 'delta_f',
        'theta_LN',
        'tc', 'phic', 'ra', 'dec', 'polarization',
        'approximant',
        'detectors'
        """

    def __init__(self, approx='phenomD', m1=10., m2=10., chi1z=0., chi2z=0., f_min=20., f_max=0., delta_f=1/64., distance=1e6 * Constants.PC_SI, fRef=0., phiRef=0.):
        """
        input:
        approx (string) : defaul = 'phenomD'
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

        # enforce m1 >= m2 and chi1 is on m1
        if m1<m2: # swap spins and masses
            chi1z, chi2z = chi2z, chi1z
            m1, m2 = m2, m1

        # put inputs into a dictionary
        self.p             = {}
        self.p['approx']   = approx
        self.p['m1']       = float(m1)
        self.p['m2']       = float(m2)
        self.p['chi1z']    = float(chi1z)
        self.p['chi2z']    = float(chi2z)
        self.p['f_min']    = float(f_min)
        self.p['f_max']    = float(f_max)
        self.p['delta_f']  = float(delta_f)
        self.p['distance'] = float(distance)
        self.p['fRef']     = float(fRef)
        self.p['phiRef']   = float(phiRef)

        self.p['Mtot'], self.p['eta'] = M_eta_m1_m2(self.p['m1'], self.p['m2'])
        self.p['chipn'] = chipn(self.p['eta'], self.p['chi1z'], self.p['chi2z'])

        self.waveform_method()

    def waveform_method(self):
        print "I am from Waveform"
        if self.p['approx'] == 'phenomD':
            import phenom
            self.ph = phenom.PhenomD(self.p)
            self.ph.IMRPhenomDGenerateFD()
        elif p['approximant'] == 'PhenomP':
            raise NotImplementedError
        else:
            raise NotImplementedError('approximant = {0} is not implemented.\
 Available approximants = {1}'.format(p['approximant'], self.available_approximants))
        # TODO: Need to add a master function to return hp and hx with appropriate
        # inclination angle and sin, cosine etc.
