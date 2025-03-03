%matplotlib inline
%config InlineBackend.figure_format = 'retina'

from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model
from scipy import interpolate

#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
pars2 = camb.CAMBparams()

#This function sets up with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=56.81393873457129, ombh2=0.0332213595801413, omch2=0.1756429346485721, mnu=0.06, omk=0, tau=0.051491469744757455, TCMB=3.0982108859530806)
pars.InitPower.set_params(As=2.0892528918873754e-09, ns=0.9675096235400634)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
pars.set_matter_power(redshifts=[0.], kmax=2.0)
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=1e-3, maxkh=10, npoints = 200)
s8 = np.array(results.get_sigma8())

print(s8)




pars2.set_cosmology(H0=67.36, ombh2=0.02237, omch2=0.1200, mnu=0.06, omk=0, tau=0.054, TCMB=2.7255)
pars2.InitPower.set_params(As=2.0989e-09, ns=0.968)
pars2.set_for_lmax(2500, lens_potential_accuracy=0);
pars2.set_matter_power(redshifts=[0.], kmax=2.0)
pars2.NonLinear = model.NonLinear_none
results2 = camb.get_results(pars2)
kh2, z2, pk2 = results2.get_matter_power_spectrum(minkh=1e-3, maxkh=10, npoints = 200)
s82 = np.array(results2.get_sigma8())
print(s82)


fig, axs = plt.subplots(2,1,figsize=(7,7), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

plt.subplots_adjust(wspace=0, hspace=0.02)
axs[0].tick_params(which = 'both',direction="in", left="off",labelleft="on")
axs[1].tick_params(which = 'both',direction="in", left="off",labelleft="on")

f = interpolate.interp1d(kh*np.exp(-0.12710492134186938)*(0.6736/0.5681393873457129)**(-1), pk[0,:]*np.exp(3*0.12710492134186938)*(0.6736/0.5681393873457129)**3)
f2 = interpolate.interp1d(kh2, pk2[0,:])
xnew = np.arange(1e-3, 1, 0.01)

ynew = f(xnew)
ynew2 = f2(xnew)


axs[0].loglog(kh*np.exp(-0.12710492134186938)*(0.6736/0.5681393873457129)**(-1), pk[0,:]*np.exp(3*0.12710492134186938)*(0.6736/0.5681393873457129)**3, color='red', label='Conformal')
axs[0].loglog(kh2, pk2[0,:], color='k', label='Standard')
axs[1].plot(xnew, ynew/ynew2, c='red', label='Conformal')
axs[1].plot(xnew, ynew/ynew, c='black', linestyle='--', linewidth=0.7, label='Conformal')


axs[0].set_ylabel('P(k)$[h^{-3}Mpc^{3}]$')
axs[1].set_xlabel('k $[h^{-1}Mpc]$')
axs[0].legend()
axs[0].set_ylim(9, 1e5)
axs[0].set_xlim(1e-3, 1)
axs[1].set_ylim(1.025, 1.055)

plt.savefig('CMB_flatlCDM_matter.pdf', dpi=500, bbox_inches='tight', pad_inches=0.1)
