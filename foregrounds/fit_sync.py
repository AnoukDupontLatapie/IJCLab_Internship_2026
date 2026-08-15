from bruit_fossil import freq_fossil, sigma_bin, sigma_bin_conv
from conversion import conv
from func_theo import PowerLawpysm, mu_dist, y_dist
from foregrounds.sync_SED import nu0_sync, data_s1_20, data_s1_80, data_s5_20, data_s5_80, data_s7_20, data_s7_80
import numpy as np
import pysm3.units as u
import matplotlib.pyplot as plt
import scipy.optimize as sp
plt.rcParams['font.family'] = 'sans-serif'

def fitPL(data):
    popt, pcov = sp.curve_fit(
    lambda freq, A, beta: PowerLawpysm(freq, A, beta, None, nu0_sync),
    freq_fossil,
    data, maxfev=10000)
    return popt, pcov

def fitPL7(data):
    popt, pcov = sp.curve_fit(
    lambda freq, A, beta, gamma: PowerLawpysm(freq, A, beta, gamma, nu0_sync),
    freq_fossil,
    data, p0 = [10**2, -1, 0.1], maxfev=10000)
    return popt, pcov

params1_20,_ = fitPL(data_s1_20)
params1_80,_ = fitPL(data_s1_80)
bestfits1_20= PowerLawpysm(freq_fossil,params1_20[0],params1_20[1],None,nu0_sync)
bestfits1_80= PowerLawpysm(freq_fossil,params1_80[0],params1_80[1],None,nu0_sync)

params5_20,_ = fitPL(data_s5_20)
params5_80,_ = fitPL(data_s5_80)
bestfits5_20= PowerLawpysm(freq_fossil,params5_20[0],params5_20[1],None,nu0_sync)
bestfits5_80= PowerLawpysm(freq_fossil,params5_80[0],params5_80[1],None,nu0_sync)

params7_20,_ = fitPL7(data_s7_20)
params7_80,_ = fitPL7(data_s7_80)
bestfits7_20= PowerLawpysm(freq_fossil,params7_20[0],params7_20[1],params7_20[2],nu0_sync)
bestfits7_80= PowerLawpysm(freq_fossil,params7_80[0],params7_80[1],params7_80[2],nu0_sync)

def resid(data, bestfit):
    """|data - bestfit| converted to MJy/sr."""
    return np.abs(data - bestfit) * conv

models = {'s1':  ('blue',   data_s1_20,  bestfits1_20,  data_s1_80,  bestfits1_80),
          's5':  ('orange', data_s5_20,  bestfits5_20,  data_s5_80,  bestfits5_80),
          's7':  ('green',  data_s7_20,  bestfits7_20,  data_s7_80,  bestfits7_80),}

for name, (color, d20, bf20, d80, bf80) in models.items():
    y20, y80 = resid(d20, bf20), resid(d80, bf80)
    plt.plot(freq_fossil, y20, c=color, label=f'{name}_20')
    plt.plot(freq_fossil, -y20, c=color, linestyle='dashed')
    plt.plot(freq_fossil, y80, c=color, linewidth=2, label=f'{name}_80')
    plt.plot(freq_fossil, -y80, c=color, linewidth=2, linestyle='dashed')
    plt.fill_between(freq_fossil, y20, y80, color=color, alpha=0.2)
    plt.fill_between(freq_fossil, -y20, -y80, color=color, alpha=0.2)

plt.xlabel(r'$\nu\,$ [GHz]')
plt.ylabel(r'$|data-model|$ [Jy/sr]')

freq_pour_dist = np.arange(57, 2001, 1) #GHz
plt.plot(freq_pour_dist, mu_dist(freq_pour_dist), color='tomato',  label='mu-distortion')
plt.plot(freq_pour_dist, -mu_dist(freq_pour_dist), color='tomato', linestyle='dashed')
plt.plot(freq_pour_dist, y_dist(freq_pour_dist), color='purple', label='y-distortion')
plt.plot(freq_pour_dist, -y_dist(freq_pour_dist), color='purple', linestyle='dashed')

plt.plot(freq_fossil, sigma_bin_conv, '.k', label='FOSSIL sensitivity')

plt.xlabel(r'$\nu\,$ [GHz]')
plt.ylabel(r'$|data-model|$ [Jy/sr]')

plt.legend(fontsize=7)
plt.loglog()