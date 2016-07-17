import unittest

class OutcomesTest(unittest.TestCase):

    def test_compare_with_reference_phenomp(self):
        import phenom
        import numpy as np

        m1 = 16.
        m2 = 12.
        chi1x = 0.6
        chi1y = 0.
        chi1z = 0.
        chi2x = 0.
        chi2y = 0.
        chi2z = 0.
        f_min = 30.
        fRef = 30.
        delta_f = 1/8.
        inclination = np.pi / 8.

        input_params = {}
        input_params.update({'m1' : m1})
        input_params.update({'m2' : m2})
        input_params.update({'chi1x' : chi1x})
        input_params.update({'chi1y' : chi1y})
        input_params.update({'chi1z' : chi1z})
        input_params.update({'chi2x' : chi2x})
        input_params.update({'chi2y' : chi2y})
        input_params.update({'chi2z' : chi2z})
        input_params.update({'f_min' : f_min})
        input_params.update({'fRef' : fRef})
        input_params.update({'inclination' : inclination})
        input_params.update({'delta_f' : delta_f})

        phenp_new = phenom.PhenomP(**input_params)
        f_new=phenp_new.flist_Hz
        a_new = np.absolute(phenp_new.hp)
        p_new = np.unwrap(np.angle(phenp_new.hp))

        f, a, p = np.loadtxt('./phenp_ref.dat').T

        try:
            self.assertTrue(len(f)==len(f_new))
        except:
            raise ValueError('new frequency series does not have the same length as the reference frequency series')

        tol = 6
        try:
            np.testing.assert_almost_equal(a, a_new, tol)
        except:
            import matplotlib
            import matplotlib.pylab as plt
            plt.figure()
            plt.plot(np.abs(a-a_new))
            plt.yscale('log')
            plt.savefig('./comp_amp.png')
            raise ValueError('new amplitude does not equal the reference amplitude')

        try:
            np.testing.assert_almost_equal(p, p_new, tol)
        except:
            import matplotlib
            import matplotlib.pylab as plt
            plt.figure()
            plt.plot(np.abs(p-p_new))
            plt.yscale('log')
            plt.savefig('./comp_phase.png')
            raise ValueError('new phase does not equal the reference phase')


# assertAlmostEqual(a, b, places=7, msg=None, delta=None)


if __name__ == '__main__':
    unittest.main()
