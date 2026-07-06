from bruit_fossil import freq_fossil, sigma_bin, sigma_bin_conv
from conversion import conv
from func_theo import MBBpysm, mu_dist, y_dist
from foregrounds.dust_SED import nu0_dust, data_d1_20, data_d1_80, data_d10_20, data_d10_80, data_d12_20, data_d12_80
import numpy as np
import pysm3.units as u
import matplotlib.pyplot as plt
import scipy.optimize as sp
plt.rcParams['font.family'] = 'sans-serif'

val = (1e6 * u.Jy/u.sr).to(u.uK_CMB,
                           equivalencies=u.cmb_equivalencies(nu0_dust*u.GHz)).value


def fitMBB(data):
    popt, pcov = sp.curve_fit(
        lambda freq, A, beta, T: MBBpysm(freq, A, beta, T, nu0_dust),
        freq_fossil,
        data,
        p0=[val, 1.5, 20.0], maxfev=5000, sigma=sigma_bin)  # valeur de paramètres à convertir en uK_CMB
    return popt, pcov


paramd1_20, _ = fitMBB(data_d1_20)
paramd1_80, _ = fitMBB(data_d1_80)
bestfitd1_20 = MBBpysm(
    freq_fossil, paramd1_20[0], paramd1_20[1], paramd1_20[2], nu0_dust)
bestfitd1_80 = MBBpysm(
    freq_fossil, paramd1_80[0], paramd1_80[1], paramd1_80[2], nu0_dust)

paramd10_20, _ = fitMBB(data_d10_20)
paramd10_80, _ = fitMBB(data_d10_80)
bestfitd10_20 = MBBpysm(
    freq_fossil, paramd10_20[0], paramd10_20[1], paramd10_20[2], nu0_dust)
bestfitd10_80 = MBBpysm(
    freq_fossil, paramd10_80[0], paramd10_80[1], paramd10_80[2], nu0_dust)

paramd12_20, _ = fitMBB(data_d12_20)
paramd12_80, _ = fitMBB(data_d12_80)
bestfitd12_20 = MBBpysm(
    freq_fossil, paramd12_20[0], paramd12_20[1], paramd12_20[2], nu0_dust)
bestfitd12_80 = MBBpysm(
    freq_fossil, paramd12_80[0], paramd12_80[1], paramd12_80[2], nu0_dust)


def resid(data, bestfit):
    """|data - bestfit| converted to MJy/sr."""
    return (data - bestfit) * conv


models = {
    'd1':  ('blue',   data_d1_20,  bestfitd1_20,  data_d1_80,  bestfitd1_80),
    'd10': ('orange', data_d10_20, bestfitd10_20, data_d10_80, bestfitd10_80),
    'd12': ('green',  data_d12_20, bestfitd12_20, data_d12_80, bestfitd12_80),
}

for name, (color, d20, bf20, d80, bf80) in models.items():
    y20, y80 = resid(d20, bf20), resid(d80, bf80)
    plt.plot(freq_fossil, y20, c=color, label=f'{name}_20')
    plt.plot(freq_fossil, -y20, c=color, linestyle='dashed')
    plt.plot(freq_fossil, y80, c=color, linewidth=2, label=f'{name}_80')
    plt.plot(freq_fossil, -y80, c=color, linestyle='dashed')
    plt.fill_between(freq_fossil, y20, y80, color=color, alpha=0.2)
    plt.fill_between(freq_fossil, -y20, -y80, color=color, alpha=0.2)

freq_pour_dist = np.arange(57, 2001, 1)  # GHz
plt.plot(freq_pour_dist, mu_dist(freq_pour_dist),
         color='tomato',  label='mu-distortion')
plt.plot(freq_pour_dist, -mu_dist(freq_pour_dist),
         color='tomato', linestyle='dashed')
plt.plot(freq_pour_dist, y_dist(freq_pour_dist),
         color='purple', label='y-distortion')
plt.plot(freq_pour_dist, -y_dist(freq_pour_dist),
         color='purple', linestyle='dashed')

plt.plot(freq_fossil, sigma_bin_conv, '.k', label='FOSSIL sensitivity')
plt.xlabel(r'$\nu\,$ [GHz]')
plt.ylabel(r'$data-model$ [Jy/sr]')

plt.legend(fontsize=7)
plt.loglog()
