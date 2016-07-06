import phenom
from phenom.utils.utils import __MSUN_SI__ as MSUN
ph = phenom.PhenomD(m1=10 * MSUN, m2=50 * MSUN, chi1z=0.3, chi2z=0.2, f_min=10)

for i in range(len(ph.pn)):
    print str(ph.pn.keys()[i]) + " = " + str(ph.pn.values()[i])

print "Now comparing to LAL"
print "Note, they will differ in the v[6] (3PN) term"
print "because in PhenomD we do not use the 3PN spin-spin term"

import lal
import lalsimulation as lalsim

lalpn = lalsim.SimInspiralTaylorF2AlignedPhasing(10, 50, 0.3, 0.2, 1, 1, 7)
# print lalpn.v

print "\n\n\n"
print "lal v"

for i in range(len(lalpn.v)):
    print str(i) + " = " + str(lalpn.v[i])

print "\n\n\n"
print "lal vlogv"

for i in range(len(lalpn.vlogv)):
    print str(i) + " = " + str(lalpn.vlogv[i])
