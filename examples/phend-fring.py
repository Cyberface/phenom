
import phenom

eta, chi1z, chi2z = 0.25, 0., 0.
fin_spin = phenom.remnant.FinalSpin0815(eta, chi1z, chi2z)

fring = phenom.remnant.fring(eta, chi1z, chi2z, fin_spin)

fdamp = phenom.remnant.fdamp(eta, chi1z, chi2z, fin_spin)

print "ringdown frequency in geometric units = ", fring
print "imaginary part of the ringdown frequency = ", fdamp

Mtot = 100. #Msol

print "ringdown frequency in Hz = ", phenom.MftoHz(fring, Mtot)



