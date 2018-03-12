import phenom

# you can convert from Hz to Mf using
#phenom.MftoHz(freq_in_mf, total_mass_in_Msun)

# get an instance of the phenomD model
phend=phenom.phenomD.PhenomD(m1=50, m2=50, chi1z=0., chi2z=0., f_min=10)

# phend.ChirpTime returns the time in seconds or in dimensionless units
# of the above system with a given start frequency (in Hz) that you
# specify below


ts = phend.ChirpTime(fStart=30, units='s')
tM = phend.ChirpTime(fStart=30, units='M')


print("the chirp time for this system is = {0} seconds".format(ts))
print("the chirp time for this system is = {0} M".format(tM))

