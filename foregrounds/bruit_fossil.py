import numpy as np
import healpy as hp
import sympy as sym
import pysm3.units as u
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'

_, f, std = np.loadtxt("foregrounds/FOSSIL_M8_STEP2_sensitivity.txt",
                       unpack=True, skiprows=9)
f /= 10 ** 9

nfreq = len(f)
points_per_bin = 2
nbins = nfreq // points_per_bin

sigma_bin_conv = np.zeros(nbins)
f_bin = np.zeros(nbins)

for i in range(nbins):
    inv_var = np.sum(1 / std[i*points_per_bin: (i+1)*points_per_bin] ** 2)
    sigma_bin_conv[i] = np.sqrt(1 / inv_var)
    f_bin[i] = np.mean(f[i*points_per_bin:(i+1)*points_per_bin])

# conversion to uK_CMB
sigma_bin = np.zeros(len(f_bin))
sigma_bin[:] = (sigma_bin_conv[:] * u.Jy/u.sr).to(u.uK_CMB,
                                                  equivalencies=u.cmb_equivalencies(f_bin*u.GHz)).value

freq_fossil = f_bin
