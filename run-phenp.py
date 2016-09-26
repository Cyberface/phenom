#to run with interactive plots us
#frameworkpython

import phenom

import matplotlib
# matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import numpy as np
from phenom.utils.utils import Constants, HztoMf

#A bad case - should be improved when new model is done
# PhenomPv2_FD(50, 0.1, -0.86, -0.86, 0.5, iota=pi/3., deltaF=0.5, f_min=10, f_max=1024)
# eta = phenom.eta_from_q(1.2)
# eta = 0.1
# print eta
# m1, m2 = phenom.m1_m2_M_eta(50., eta)
# print m1, m2
# chi1x = 0.5
# chi1z = -0.86
# chi2z = -0.86
# # chi1z = 0.
# # chi2z = 0.
# f_min = 10.
# fRef=0.



# Set I used to test with
# m1 = 12.
# m2 = 8.31
# chi1x = 0.8
# chi1y = 0.
# chi1z = 0.2
# chi2x = 0.
# chi2y = 0.
# chi2z = 0.4
# f_min = 30.
# fRef = 30.
# delta_f = 1/16.
#
# inclination = np.pi / 8.

m1 = 12.
# m2 = 8.31
m2 = 16.
chi1x = 0.6   # chi1x = 0.8
chi1y = 0.   # chi1y = 0.
chi1z = 0.   # chi1z = 0.2
chi2x = 0.   # chi2x = 0.
chi2y = 0.   # chi2y = 0.
chi2z = 0.   # chi2z = 0.4
f_min = 30.
fRef = 30.
# delta_f = 1/64.
delta_f = 1/32.


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

phenp = phenom.PhenomP(**input_params)

#get phenomD waveform
# phend = phenom.Waveform(approximant='IMRPhenomD',m1=m1, m2=m2, delta_f=1/2., f_min=10, inclination=inc)
phend = phenom.Waveform(approximant='IMRPhenomD', **input_params)

print phenp.p['omega_Ref']


fl_Hz = np.arange(phenp.p['f_min'], phenp.p['f_max'], phenp.p['delta_f'])
#NOTE: CONVERTED TO Mf
fl_Hz = HztoMf(fl_Hz, phenp.p['Mtot'])

# fl_Mf = phenom.HztoMf(fl_Hz, phenp.p['Mtot'])

#
#time correction from lalsimulation
#
# t_corr = -0.00533
# t_corr = 0.00533
# # phase_corr = np.exp(-2*Constants.LAL_PI * 1.j * fl_Hz * t_corr * Constants.MTSUN_SI*phenp.p['Mtot'])
# phase_corr = np.exp(-2*Constants.LAL_PI * 1.j * fl_Hz * t_corr)
# phenp.hpold = phenp.hp
# phenp.hcold = phenp.hc
# phenp.hp *= phase_corr
# phenp.hc *= phase_corr

plt.figure()
plt.plot(fl_Hz, phenp.alpha)
# plt.plot(fl_Mf, phenp.alpha)
# plt.xscale('log')
plt.savefig('alpha.png')
plt.close()

plt.figure()
plt.plot(fl_Hz, phenp.epsilon)
# plt.plot(fl_Mf, phenp.epsilon)
# plt.xscale('log')
plt.savefig('epsilon.png')
plt.close()


#get the correct amp scaling for 22 mode, aligned spin.
# amp0 = 2. * np.sqrt(5. / (64.*np.pi)) * phenp.p['Mtot'] * Constants.MRSUN_SI * phenp.p['Mtot']
#using the aligned spin waveform from phenomp i.e. phenp.hP doesn't have the correct scaling.
#easier to just call Waveform(approximant="IMRPhenomD")


plt.figure()
plt.plot(fl_Hz, np.abs(phenp.hp), label='hp')
plt.plot(fl_Hz, np.abs(phenp.hc), label='hc')
# plt.plot(fl_Hz, amp0*np.abs(phenp.hP), label='h_aligned')
plt.plot(HztoMf(phend.flist_Hz, phend.input_params['m1']+phend.input_params['m2']), np.absolute(phend.hptilde), label='phenomD:hp', ls='--')
plt.plot(HztoMf(phend.flist_Hz, phend.input_params['m1']+phend.input_params['m2']), np.absolute(phend.hctilde), label='phenomD:hc', ls='--')
plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
plt.savefig('phenp_1.png')
plt.close()


#plot phenomP from LAL

import lal
import lalsimulation as lalsim


_lshp, _lshc = lalsim.SimInspiralChooseFDWaveform(0, delta_f, m1*lal.MSUN_SI, m2*lal.MSUN_SI,
                chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, f_min, 0, fRef, 1e6*phenom.Constants.PC_SI, inclination,
                0, 0, None, None, -1, -1, lalsim.IMRPhenomPv2)

lsf = np.arange(_lshp.data.length) * _lshp.deltaF

#NOTE: CONVERTED TO Mf
lsf = HztoMf(lsf, m1+m2)


lshp = _lshp.data.data
lshc = _lshc.data.data


plt.figure()
plt.plot(fl_Hz, np.abs(phenp.hp), label='hp-phenom')
plt.plot(lsf, np.abs(lshp), label='hp-lalsim')
# plt.plot(fl_Hz, np.abs(phenp.hc), label='hc-phenom')
# plt.plot(lsf, np.abs(lshc), label='hc-lalsim')
# plt.plot(fl_Hz, np.abs(phenp.hP), label='h_aligned')
plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
# plt.xlim(10,3000)
# plt.show()
plt.savefig('lal_phenp_1.png')
plt.close()


#comparison
#interpolate lalsimulation data onto phenp frequencies: fl_Hz
from scipy.interpolate import interp1d
lshp_i = interp1d(lsf, lshp)
lshc_i = interp1d(lsf, lshc)

lshp_flHz = lshp_i(fl_Hz)
lshc_flHz = lshc_i(fl_Hz)

plt.figure()
# plt.plot(fl_Hz, np.abs(np.abs(phenp.hp) - np.abs(lshp_flHz))/np.abs(lshp_flHz), label='hp: phenom-lalsim')
plt.plot(fl_Hz, np.abs(np.abs(phenp.hp) - np.abs(lshp_flHz)), label='hp: phenom-lalsim')
plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
# plt.xlim(10,1000)
# plt.show()
plt.savefig('lal_phenp_1_ampdiff.png')
plt.close()


phenp_phase = np.unwrap(np.angle(phenp.hp))
lalsim_phase = np.unwrap(np.angle(lshp_flHz))

phenp_phase_hc = np.unwrap(np.angle(phenp.hc))
lalsim_phase_hc = np.unwrap(np.angle(lshc_flHz))

print "initial phase phenp = ", phenp_phase[0]
print "initial phase lalsim_phase = ", lalsim_phase[0]

# shift = phenp_phase[0] - lalsim_phase[0]
# print "shift = ", shift
# phenp_phase -= shift

#subtract off initial phase form each
phenp_phase -= phenp_phase[0]
lalsim_phase -= lalsim_phase[0]

print "subtracted** initial phase phenp = ", phenp_phase[0]
print "subtracted** initial phase lalsim_phase = ", lalsim_phase[0]

phenp_phase_hc -= phenp_phase_hc[0]
lalsim_phase_hc -= lalsim_phase_hc[0]



np.savetxt('/Users/sebastian/Desktop/phenp_phase.dat', np.column_stack([fl_Hz, phenp_phase]))
np.savetxt('/Users/sebastian/Desktop/lalsim_phase.dat', np.column_stack([fl_Hz, lalsim_phase]))

np.savetxt('/Users/sebastian/Desktop/phenp_beta.dat', np.column_stack([fl_Hz, phenp.beta]))
np.savetxt('/Users/sebastian/Desktop/phenp_cexp_i_alpha.dat', np.column_stack([fl_Hz, np.real(phenp.cexp_i_alpha), np.imag(phenp.cexp_i_alpha)]))

np.savetxt('/Users/sebastian/Desktop/phenp_alpha.dat', np.column_stack([fl_Hz, phenp.alpha]))
np.savetxt('/Users/sebastian/Desktop/phenp_epsilon.dat', np.column_stack([fl_Hz, phenp.epsilon]))

# np.savetxt('/Users/sebastian/Desktop/phenp_beta.dat', phenp.beta)


plt.figure()
plt.plot(fl_Hz, phenp_phase, label='hp: phenom')
plt.plot(fl_Hz, lalsim_phase, label='hp: lalsim')
plt.legend(loc='best')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim(10,1000)
# plt.show()
plt.savefig('lal_phenp_1_phase.png')
plt.close()

plt.figure()
plt.plot(fl_Hz, np.abs( (phenp_phase - lalsim_phase) ), label='hp: phenom-lalsim')
plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
# plt.xlim(10,1000)
# plt.show()
plt.savefig('lal_phenp_1_phasediff.png')
plt.close()

plt.figure()
plt.plot(fl_Hz, np.abs( (phenp_phase_hc - lalsim_phase_hc) ), label='hc: phenom-lalsim')
plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
# plt.xlim(10,1000)
# plt.show()
plt.savefig('lal_phenp_hc_phasediff.png')
plt.close()


plt.figure()
plt.plot(fl_Hz, phenp.dphilist, label='phase-deriv-phenp')
plt.legend(loc='best')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim(10,1000)
plt.ylim(-10*np.abs(phenp.t_corr),10*np.abs(phenp.t_corr))
plt.axvline(x=phenp.MfRD, linewidth=2, color='r')
# plt.show()
plt.savefig('lal_phenp_dphi.png')
plt.close()





#
