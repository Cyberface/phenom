import phenom

#get instances of phenomD
# ph1 = phenom.PhenomD(m1=100, m2=10, chi1z=0., f_min=10, delta_f=1./128.)
wf1 = phenom.Waveform(approximant='IMRPhenomD', m1=100, m2=10, chi1z=0., f_min=10, delta_f=1./128.)
ph2 = phenom.PhenomD(m1=100, m2=10.1, chi1z=0., f_min=10, delta_f=1./128.)

#generate htilde
# ph1.IMRPhenomDGenerateFD()
ph2.IMRPhenomDGenerateFD()

#load psd function
# psd_fun = phenom.read_and_interp_psd_from_txt(phenom.PHENOM_PACKAGE_PATH + '/psd/data/ZERO_DET_high_P.txt')
psd_fun = phenom.read_and_interp_psd_from_txt('/Users/sebastian/git/phenom/phenom' + '/psd/data/ZERO_DET_high_P.txt')

#get instance of Match
match = phenom.Match()

# print "match = ", match.match(ph1, ph2)
# print "match = ", match.match(ph1, ph2, 9, 100)
# print "match = ", match.match(ph1, ph2, fmin=10, fmax=0, psd_fun=psd_fun)
#
import time

t0 = time.time()
print "[without psd] match = ", match.match(wf1, ph2, fmin=10, fmax=0, psd_fun=None)
print "[with psd] match = ", match.match(wf1, ph2, fmin=10, fmax=0, psd_fun=psd_fun)
t1 = time.time()

total = t1-t0
print "time taken = {0} seconds".format(total)

# import numpy as np
# m2list = np.linspace(9.9,10.1,100)
# phlist = [None] * len(m2list)
# mlist = [None] * len(m2list)
# for i in range(len(phlist)):
#     phlist[i] = phenom.PhenomD(m1=100, m2=m2list[i], chi1z=0., f_min=10, delta_f=1./16.)
#     phlist[i].IMRPhenomDGenerateFD()
#     mlist[i] = match.my_match(ph1, phlist[i], fmin=10, fmax=0, psd_fun=psd_fun)
#
# import matplotlib
# matplotlib.use('MacOSX')
# import matplotlib.pyplot as plt
#
#
# plt.figure()
# plt.plot(m2list, mlist)
# plt.show()
