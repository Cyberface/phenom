#this file contains the code to interpolate and compute
#the alpha and epsilon precession angles
#from the phenEOB model.

#phenEOB is a phenomenological model of the precession angles
#as predicted by SEOBNRv3
#as a function of frequency

import numpy as np


class InitialisePhenEOBModel(object):
    """This functions sets up the interpolants.
    We only want to do this once when we call phenomP"""
    def __init__(self, model_name):
        from phenom.waveform.phenEOButils import InterpolateAlpha
        self.model_name = model_name
        self.model_paths = self._preset_models(self.model_name)
        self.ia = InterpolateAlpha(self.model_paths['h5paths'], self.model_paths['h5dict_names'])

    def _preset_models(self, model_name):

        # rootdir = '/Users/sebastian/phd/mounts/arcca_mount/'
        rootdir = '/home/spx8sk/'

        if model_name == 'grid5x6step':
            grid5x6step= {}

            grid5x6step['q1h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q1-spin-grid-5x6/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid5x6step['q2h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q2-spin-grid-5x6/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid5x6step['q3h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q3-spin-grid-5x6/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid5x6step['q4h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q4-spin-grid-5x6/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid5x6step['q5h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q5-spin-grid-5x6/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid5x6step['q6h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q6-spin-grid-5x6/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid5x6step['q7h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q7-spin-grid-5x6/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid5x6step['q8h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q8-spin-grid-5x6/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'

            grid5x6step['h5paths'] = [
                grid5x6step['q1h5path0'],
                grid5x6step['q2h5path0'],
                grid5x6step['q3h5path0'],
                grid5x6step['q4h5path0'],
                grid5x6step['q5h5path0'],
                grid5x6step['q6h5path0'],
                grid5x6step['q7h5path0'],
                grid5x6step['q8h5path0']
            ]

            grid5x6step['h5dict_names'] = ['q1', 'q2', 'q3', 'q4', 'q5', 'q6', 'q7', 'q8']

            ret = grid5x6step
        elif model_name == 'grid20x20step' or model_name == 'grid20x20step_ep_eq_al' or model_name == 'grid20x20step_FITBETA':
            grid20x20step= {}

            grid20x20step['q1h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q1-spin-grid-20x20/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid20x20step['q2h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q2-spin-grid-20x20/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid20x20step['q3h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q3-spin-grid-20x20/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid20x20step['q4h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q4-spin-grid-20x20/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid20x20step['q5h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q5-spin-grid-20x20/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid20x20step['q6h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q6-spin-grid-20x20/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid20x20step['q7h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q7-spin-grid-20x20/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'
            grid20x20step['q8h5path0'] = rootdir + 'projects/phenEOB/condor_runs/q8-spin-grid-20x20/bf_model_coeffs/p1_m1_p2_m2_None_0.04_0.02_0.1.h5'

            grid20x20step['h5paths'] = [
                grid20x20step['q1h5path0'],
                grid20x20step['q2h5path0'],
                grid20x20step['q3h5path0'],
                grid20x20step['q4h5path0'],
                grid20x20step['q5h5path0'],
                grid20x20step['q6h5path0'],
                grid20x20step['q7h5path0'],
                grid20x20step['q8h5path0']
            ]

            grid20x20step['h5dict_names'] = ['q1', 'q2', 'q3', 'q4', 'q5', 'q6', 'q7', 'q8']

            ret = grid20x20step
        else:
            print("model_name = {0} not recognised".format(model_name))

        return ret

class phenEOBalpha(object):
    def __init__(self, interpolated_alpha):

        self.interpolated_alpha = interpolated_alpha

        self.p1_func = self._get_ansatz(self.interpolated_alpha, 'part1')
        self.p2_func = self._get_ansatz(self.interpolated_alpha, 'part2')

    def _get_ansatz(self, ia, part):
        """given part:string
        return the function/ansatz used to fit the data.
        This is thus the function needed to evaluate the
        intepolated coefficients"""
        if part == "part1":
            func = ia.bfp['q1'].p1_model_func
        elif part == "part2":
            func = ia.bfp['q1'].p2_model_func
        return func

    def get_coeffs(self, q, chimag, theta):

        p1_coeff_names = ['p1_0', 'p1_1', 'p1_2', 'p1_3']
        p2_coeff_names = ['p2_0', 'p2_1']

        p1_coeff_values = []
        for c_name in p1_coeff_names:
            p1_coeff_values.append( self.interpolated_alpha.bfturps[c_name](chimag, theta, q) )

        p2_coeff_values = []
        for c_name in p2_coeff_names:
            p2_coeff_values.append( self.interpolated_alpha.bfturps[c_name](chimag, theta, q) )

        return np.asarray(p1_coeff_values), np.asarray(p2_coeff_values)

    def alpha_at_any_omega_ref(self, omega, q, chi1x, chi1z, pnorder=-1):
        from phenom import CartToPolar
        chimag, theta = CartToPolar(chi1x, chi1z)

        #check boundary for model
        chimag_low = 0.01
        chimag_high = 0.99
        theta_low = 2.*np.pi/180.
        theta_high = 178*np.pi/180.


        #TODO!
        #Add propery boundary checks



#         chi1x_low, chi1z_low = PolarToCart(chimag_low, theta_low)
#         chi1x_high, chi1z_high = PolarToCart(chimag_high, theta_high)
#         print chi1x_low, chi1z_low
#         print chi1x_high, chi1z_high

        # commented out 23/8/16
        if np.sqrt(chi1x**2.+chi1z**2.) < chimag_low:
            print "np.sqrt(chi1x**2.+chi1z**2.) = {0} is invalid. lower boundary = {1}".format(np.sqrt(chi1x**2.+chi1z**2.), chimag_low)
            import sys
            sys.exit(1)
        if np.sqrt(chi1x**2.+chi1z**2.) > chimag_high:
            print "np.sqrt(chi1x**2.+chi1z**2.) = {0} is invalid. upper boundary = {1}".format(np.sqrt(chi1x**2.+chi1z**2.), chimag_high)
            import sys
            sys.exit(1)


        p1coeffs, p2coeffs = self.get_coeffs(q, chimag, theta)

        xjoin = 0.04
        #need to correct for xjoin...
        if omega <= xjoin:
            #then use part 1
            alpha = self.p1_func(omega, q, chi1x, chi1z, pnorder, *p1coeffs)
        elif omega > xjoin:
            #then use part 2
            alpha = self.p2_func(omega, *p2coeffs)

            #simply translational alignment at xjoin
            #I want this function to be callable at any frequency so that it can
            #be used in the ROQs
            #so applying the shift like this

            p1_alpha_xjoin = self.p1_func(xjoin, q, chi1x, chi1z, pnorder, *p1coeffs)
            p2_alpha_xjoin = self.p2_func(xjoin, *p2coeffs)

            shift = p1_alpha_xjoin - p2_alpha_xjoin
            alpha += shift
        else:
            print("Should never get to here")

        return alpha


    # def alpha_at_range_omgea(self, omega, q, chi1x, chi1z, pnorder=-1):




# class phenEOBepsilon(object):
#     def __init__(self):
#         pas
#     def epsilon_at_reference(self):
#         """need to be able to compute the values at any frequency, even
#         those not in list"""
#         pass
