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
    # static element
    default_args = {'approximant':'IMRPhenomD',
                    'm1':50.,
                    'm2':50.,
                    'f_min':30.,
                    'chi1x':0.,
                    'chi1y':0.,
                    'chi1z':0.,
                    'chi2x':0.,
                    'chi2y':0.,
                    'chi2z':0.,
                    'inclination':0.,
                    'distance':1e6 * Constants.PC_SI,
                    'f_max':0.,
                    'fRef':0.,
                    'phiRef':0.,
                    'delta_f':1./2.,
                    'delta_t':1./2.}

    available_approximants = ['IMRPhenomD']

    def __init__(self,
        approximant=default_args['approximant'],
        m1=default_args['m1'],
        m2=default_args['m2'],
        chi1x=default_args['chi1x'],
        chi1y=default_args['chi1y'],
        chi1z=default_args['chi1z'],
        chi2x=default_args['chi2x'],
        chi2y=default_args['chi2y'],
        chi2z=default_args['chi2z'],
        inclination=default_args['inclination'],
        distance=default_args['distance'],
        f_min=default_args['f_min'],
        f_max=default_args['f_max'],
        fRef=default_args['fRef'],
        phiRef=default_args['phiRef'],
        delta_f=default_args['delta_f'],
        delta_t=default_args['delta_t']):
        """
        input:
        approx (string) : default = 'IMRPhenomD'
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
        #TODO: automate the names of the fields with values
        self.input_params                = {}
        self.input_params['approximant'] = approximant
        self.input_params['m1']          = float(m1)
        self.input_params['m2']          = float(m2)
        self.input_params['chi1x']       = float(chi1x)
        self.input_params['chi1y']       = float(chi1y)
        self.input_params['chi1z']       = float(chi1z)
        self.input_params['chi2x']       = float(chi2x)
        self.input_params['chi2y']       = float(chi2y)
        self.input_params['chi2z']       = float(chi2z)
        self.input_params['inclination'] = float(inclination)
        self.input_params['distance']    = float(distance)
        self.input_params['f_min']       = float(f_min)
        self.input_params['f_max']       = float(f_max)
        self.input_params['fRef']        = float(fRef)
        self.input_params['phiRef']      = float(phiRef)
        self.input_params['delta_f']     = float(delta_f)
        self.input_params['delta_t']     = float(delta_t)
        #TODO: It would be nicer if only variables that are needed
        #are defined. For examples only time domain approximants need
        #a 'delta_t'. It doesn't really mean anything for frequency domain.

        print "calling _generate_fd_waveform"
        self._generate_fd_waveform(self.input_params)

    def _generate_fd_waveform(self, input_params):
        """
        here is where the main driver waveform functions are called
        from the 'Waveform' class.
        When a new approximant is added, add a conditional here.
        Also add it to the 'available_approximants' list at the top.

        input:
            input_params : dict
        """
        if input_params['approximant'] == 'IMRPhenomD':
            from phenom.waveform import PhenomD
            self.ph = PhenomD(\
                        m1=input_params['m1'], m2=input_params['m2'],\
                        chi1z=input_params['chi1z'], chi2z=input_params['chi2z'],\
                        f_min=input_params['f_min'], f_max=input_params['f_max'],\
                        delta_f=input_params['delta_f'],\
                        distance=input_params['distance'],\
                        fRef=input_params['fRef'], phiRef=input_params['phiRef'])
            self.ph.IMRPhenomDGenerateFD()

        elif input_params['approximant'] == 'IMRPhenomP':
            raise NotImplementedError('PhenomP not implemented')

        else:
            raise NotImplementedError("approximant = {0} is not implemented. Available approximants = {1}".format(input_params['approximant'], self.available_approximants))
        # TODO: Need to add a master function to return hp and hx with appropriate
        # inclination angle and sin, cosine etc.
