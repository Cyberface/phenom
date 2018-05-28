import numpy as np
from phenom.testing import prec_angle_helper as pah
from phenom import HztoMf, m1_m2_M_q
import lal

f_gw_min=1.
f_gw_max=700
df_gw=0.1
f_gw_ref=f_gw_min*2

Npts = int(np.ceil((f_gw_max - f_gw_min)/df_gw))


f_gw_list = lal.CreateREAL8Sequence(Npts)
f_gw_list.data = np.arange(f_gw_min, f_gw_max, df_gw)


# m1=60*lal.MSUN_SI/10.
# m2=10*lal.MSUN_SI/10.
# s1x=0.7393
# s1y=0.1902
# s1z=-0.4953
# s2x=-0.1896
# s2y=5.103e-2
# s2z=-0.2268

# s1x=0.
# s1y=0
# s1z=0
# s2x=0
# s2y=0
# s2z=0

# s1x=1.
# s1y=0
# s1z=0
# s2x=1.
# s2y=0
# s2z=0



m1,m2 =m1_m2_M_q(20., 3)

m1=m1*lal.MSUN_SI
m2=m2*lal.MSUN_SI
s1x=1.
s1y=0
s1z=0
s2x=0.
s2y=0
s2z=0

pv3pars = {
    "flist":f_gw_list,
    "m1":m1,
    "m2":m2,
    "s1x":s1x,
    "s1y":s1y,
    "s1z":s1z,
    "s2x":s2x,
    "s2y":s2y,
    "s2z":s2z,
    "fref":f_gw_ref,
    "ExpansionOrder":5,
    "PN":"3"
}

pv3angles={}
pv3angles['alpha'],pv3angles['beta'],pv3angles['epsilon']= pah.evaluate_phenomPv3_angles(**pv3pars)


pv2pars = {
    "flist":f_gw_list,
    "m1":m1,
    "m2":m2,
    "s1x":s1x,
    "s1y":s1y,
    "s1z":s1z,
    "s2x":s2x,
    "s2y":s2y,
    "s2z":s2z,
    "fref":f_gw_ref
}

pv2angles={}
pv2angles['alpha'],pv2angles['beta'],pv2angles['epsilon']= pah.evaluate_phenomPv2_angles(**pv2pars)

# mf_gw_list = HztoMf(f_gw_list.data, (m1+m2)/lal.MSUN_SI)

import matplotlib.pyplot as plt

# plt.figure()
# plt.plot(f_gw_list.data, pv2angles['alpha'], label='pv2')
# plt.plot(f_gw_list.data, pv3angles['alpha'], label='pv3')
# plt.title("alpha")
# plt.legend()
# plt.xscale('log')
# plt.show()

# plt.figure()
# plt.plot(f_gw_list.data, pv2angles['beta'], label='pv2')
# plt.plot(f_gw_list.data, pv3angles['beta'], label='pv3')
# plt.title("beta")
# plt.legend()
# plt.yscale('log')
# plt.xscale('log')
# plt.show()


# plt.figure()
# plt.plot(f_gw_list.data, pv2angles['epsilon'], label='pv2')
# plt.plot(f_gw_list.data, pv3angles['epsilon'], label='pv3')
# plt.title("epsilon")
# plt.legend()
# plt.xscale('log')
# plt.show()

fig, axes = plt.subplots(1,3, figsize=(14,4))

axes[0].plot(f_gw_list.data, pv2angles['alpha'], label='pv2')
axes[0].plot(f_gw_list.data, pv3angles['alpha'], label='pv3')
axes[0].set_title("alpha")
axes[0].set_xscale('log')

axes[1].plot(f_gw_list.data, pv2angles['beta'], label='pv2')
axes[1].plot(f_gw_list.data, pv3angles['beta'], label='pv3')
axes[1].set_title("beta")
axes[1].set_yscale('log')
axes[1].set_xscale('log')

axes[2].plot(f_gw_list.data, pv2angles['epsilon'], label='pv2')
axes[2].plot(f_gw_list.data, pv3angles['epsilon'], label='pv3')
axes[2].set_title("epsilon")
axes[2].set_xscale('log')
axes[2].legend()
plt.show()
