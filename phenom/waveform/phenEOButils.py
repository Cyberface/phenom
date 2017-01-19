
# this file `phenEOButils` is the home of the `BestfitPrediction` and `InterpolateAlpha` classes
# Abbriviation: InterpolateAlpha = `ia`

#BestfitPrediction organises the coefficients
#InterpolateAlpha interpolates them


import numpy as np
import h5py
from scipy.interpolate import LinearNDInterpolator


class BestfitPrediction(object):
    """
    h5path : string. Path to h5 file containing model coefficients
    modelpath : string. Path to models.py in phenEOB git repo.
    """
    def __init__(self, h5path, modelpath):

        #need to append to path so that I can import the models.py
        self.modelpath = modelpath
        import sys
        if self.modelpath not in sys.path:
            sys.path.append(self.modelpath)
        # from models import *
        import models

        self.h5path = h5path

        self.p1_bf, self.p2_bf = self._get_bf_coeffs(self.h5path)

        self.p1_model_name, self.p2_model_name, self.p1_model_func, self.p2_model_func = self._get_funcs(self.h5path, self.modelpath)

        self.q, self.chimag, self.theta, self.chi1x, self.chi1z = self._get_physpars(self.h5path)

    def _get_bf_coeffs(self, h5path):
        f = h5py.File(h5path, 'r')
        #parts 1 and 2 best fit coefficients
        #turn best fit coefficients into nice arrays where each row corresponds
        #to each calibration waveform.
        p1_bf = []
        p2_bf = []
        for p1name in f['p1'].keys():
            p1_bf.append(f['p1'][p1name].value)
        for p2name in f['p2'].keys():
            p2_bf.append(f['p2'][p2name].value)

        #p1_bf and p2_bf have numpy array shapes
        #(num waveforms, num coeffs)
        p1_bf = np.asarray(p1_bf).T
        p2_bf = np.asarray(p2_bf).T

        f.close()

        return p1_bf, p2_bf

    def _get_funcs(self, h5path, modelpath):
        #annoyingly I had to put this import here so tht it would import models
        import sys
        if modelpath not in sys.path:
            sys.path.append(modelpath)
        import models

        f = h5py.File(h5path, 'r')

        #names of the models used to make fits
        p1_model_name = f['models']['p1'].value
        p2_model_name = f['models']['p2'].value

        #using the GLOBAL_FUNC_DICT_P1/2 get actual python functions
        p1_model_func = models.GLOBAL_FUNC_DICT_P1[p1_model_name]
        p2_model_func = models.GLOBAL_FUNC_DICT_P2[p2_model_name]

        f.close()

        return p1_model_name, p2_model_name, p1_model_func, p2_model_func

    def _get_physpars(self, h5path):
        f = h5py.File(h5path, 'r')

        physpars = f['physpars']
        q, chimag, theta, chi1x, chi1z = (physpars['physpars'].value).T

        f.close()

        return q, chimag, theta, chi1x, chi1z

    def prediction(self, x1, x2, dx, xjoin, i):
        """
        uses best-fit coefficients
        """
        q, chi1x, chi1z = self.q[i], self.chi1x[i], self.chi1z[i]
        p1x = np.arange(x1, xjoin, dx)
        p2x = np.arange(xjoin, x2, dx)

        #predition part 1
        pnorder=-1
        p1y = self.p1_model_func(p1x, q, chi1x, chi1z, pnorder, *self.p1_bf[i])

        #predition part 2
        p2y = self.p2_model_func(p2x, *self.p2_bf[i])

        #simply translational alignment
        shift = p1y[-1] - p2y[0]
        p2y += shift

        x = np.concatenate((p1x, p2x[1:])) #OmegaOrb
        y = np.concatenate((p1y, p2y[1:])) #alpha

        return x, y




class InterpolateAlpha(object):
    def __init__(self, h5paths, h5dict_names):
        # self.modelpath = '/home/spx8sk/git/phenEOB/convergence/2d/'
        # self.modelpath = '/Users/sebastian/git/phenEOB/convergence/2d/'
        self.modelpath = '/Users/sebastian/work/git/phenEOB/convergence/2d/'
        #bfp dict is the 'best fit prediction' dictionary
        #it contains the calibration points and the values of the best fit coefficients
        #each entry is for a particular mass-ratio spin-grid which corresponds to
        #a particular h5 file in the 'h5paths' input list.

        self.h5paths = h5paths
        self.h5dict_names = h5dict_names

        self.bfp = self._populate_bfp(self.h5paths, self.h5dict_names)

        #get number of best-fit coeffs. Will be useful later
        self.p1_num_bf_coeffs = len(self.bfp[self.h5dict_names[0]].p1_bf[0])
        self.p2_num_bf_coeffs = len(self.bfp[self.h5dict_names[0]].p2_bf[0])

#         print self.p1_num_bf_coeffs
#         print self.p2_num_bf_coeffs

        self.grid = self.construct_ThreeD_grid(self.h5dict_names, self.bfp)

        #bfc stands for best-fit coefficients.
        #this is a dict with keys prefixed with either  'p1_' or 'p2_' for part 1 or part 2.
        #then the coefficient number eg 'p1_0' is the first coefficient for part 1.
        #this is then a 1D numpy array for all the mass-ratio spin-grids essentially concatonated
        self.bfc = self.construct_bestfit_coeff_values(self.bfp, self.h5dict_names, self.p1_num_bf_coeffs, self.p2_num_bf_coeffs)


        #dictionary of interpolants. labeled the same as self.bfc
        #each interpolant is a function of (chimag, theta, q)
        self.bfturps = self.construct_interpolants(self.grid, self.bfc, self.p1_num_bf_coeffs, self.p2_num_bf_coeffs)

        #compute residuals of interpolant
        # self.residuals = self.compute_residuals(self.grid, self.bfc, self.p1_num_bf_coeffs, self.p2_num_bf_coeffs, self.bfturps)


    def _populate_bfp(self, h5paths, h5dict_names):
        bfp={}
        for i, h5path in enumerate(h5paths):
            bfp[h5dict_names[i]] = BestfitPrediction(h5path, self.modelpath)
        return bfp

    def construct_ThreeD_grid(self, h5dict_names, bfp):
        #next construct the 3D grid of q, chi1x, chi1z which will be the interpolation nodes

        #setup empty lists
        chimag_array = []
        theta_array = []
        q_array = []
        #append all the calibration points
        for h5dict_name in h5dict_names:
            chimag_array.append(bfp[h5dict_name].chimag)
            theta_array.append(bfp[h5dict_name].theta)
            q_array.append(bfp[h5dict_name].q)
        #flatten them to achieve concatenation
        chimag_array = np.asarray(chimag_array).flatten()
        theta_array = np.asarray(theta_array).flatten()
        q_array = np.asarray(q_array).flatten()


        # I NEED TO DO THIS TO AVOID nan returns when evaluating the interpolant
        for i, chi in enumerate(chimag_array):
            chimag_array[i] = float("{0:.5f}".format(chi))
        for i, theta in enumerate(theta_array):
            theta_array[i] = float("{0:.5f}".format(theta))
        for i, q in enumerate(q_array):
            q_array[i] = float("{0:.4f}".format(q))


        #correct 3D grid is obtained via a column_stack
        ThreeDgrid = np.column_stack(
            (
                chimag_array,
                theta_array,
                q_array
            )
        )

        return ThreeDgrid

    def plot_calibration_points(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(*(self.grid.T))


    def construct_bestfit_coeff_values(self, bfp, h5dict_names, p1_num_bf_coeffs, p2_num_bf_coeffs):
        #this dict is labeled by the coefficient number '0', '1', ..., upto either p1_num_bf_coeffs or p2_num_bf_coeffs
        #each entry in the dict is a 1D array for each best-fit coefficient
        #the same is achieved by concatonating each a particular best-fit coefficient for each spin-grid
        bfcoeffs = {}

        #get part 1 coefficients
        for i in range(p1_num_bf_coeffs):
            p1_bf_array = []
            for h5dict_name in h5dict_names:
                p1_bf_array.append(bfp[h5dict_name].p1_bf.T[i])
            p1_bf_array = np.asarray(p1_bf_array).flatten()
            bfcoeffs['p1_' + str(i)] = p1_bf_array

        #get part 2 coefficients
        for i in range(p2_num_bf_coeffs):
            p2_bf_array = []
            for h5dict_name in h5dict_names:
                p2_bf_array.append(bfp[h5dict_name].p2_bf.T[i])
            p2_bf_array = np.asarray(p2_bf_array).flatten()
            bfcoeffs['p2_' + str(i)] = p2_bf_array

        return bfcoeffs


    def construct_interpolants(self, grid, bfc, p1_num_bf_coeffs, p2_num_bf_coeffs):
        #bfturps is a dictionary labeled in the same way as bfcoeffs.
        #i.e. 'p1_0' is the first coefficient for part 1.
        #this then returns a 3D interpolant indexed by (chimag, theta, q)
        #which predicts the particular coefficient
        bfturps = {}

        #part 1 interpolants
        for i in range(p1_num_bf_coeffs):
            bfturps['p1_' + str(i)] = LinearNDInterpolator(grid, bfc['p1_' + str(i)])

        #part 2 interpolants
        for i in range(p2_num_bf_coeffs):
            bfturps['p2_' + str(i)] = LinearNDInterpolator(grid, bfc['p2_' + str(i)])

        return bfturps

    def compute_residuals(self, grid, bfc, p1_num_bf_coeffs, p2_num_bf_coeffs, bfturps):
        residuals = {}

        num_calibration_points = grid.shape[0]

        #part 1 interpolants
        for i in range(p1_num_bf_coeffs):
            res=[]
            for j in range(num_calibration_points):
                res.append( bfturps['p1_' + str(i)](*grid[j]) - bfc['p1_' + str(i)][j] )
            residuals['p1_' + str(i)] = np.asarray(res)

        #part 2 interpolants
        for i in range(p2_num_bf_coeffs):
            res=[]
            for j in range(num_calibration_points):
                res.append( bfturps['p2_' + str(i)](*grid[j]) - bfc['p2_' + str(i)][j] )
            residuals['p2_' + str(i)] = np.asarray(res)

        return residuals
