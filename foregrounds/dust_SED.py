
from bruit_fossil import freq_fossil
from conversion import conv
from func_theo import MBBpysm
from mask_file_download import mask_file
import numpy as np
import healpy as hp
import pysm3.units as u
import pysm3
import matplotlib.pyplot as plt
from astropy import constants
import gc
plt.rcParams['font.family'] = 'sans-serif'


def map_pysm(fgtype, freq_reference):  # Generate map at instrument frequencies

    sky = pysm3.Sky(nside=nside, preset_strings=['%s' % fgtype])

    dust = sky.components[0]
    if fgtype == 'd12':
        betamap = dust.mbb_index.value  # get the beta map
        tempmap = dust.mbb_temperature.value  # get the temp map
        nu0 = dust.freq_ref.value  # get reference frequency
        AI = dust.layers.to(
            u.uK_CMB, equivalencies=u.cmb_equivalencies(nu0*u.GHz)).value
    elif fgtype == 'd0':
        betamap = np.ones(npix)*np.array([dust.mbb_index.value])
        betamap = np.reshape(betamap, shape=(1, npix))
        tempmap = np.ones(npix)*np.array([dust.mbb_temperature.value])
        tempmap = np.reshape(tempmap, shape=(1, npix))
        nu0 = dust.freq_ref_I.value
        AI = np.array([dust.get_emission(freq_reference * u.GHz).to(u.uK_CMB,
                                                                    equivalencies=u.cmb_equivalencies(freq_reference*u.GHz)).value])
    else:
        nu0 = dust.freq_ref_I.value
        betamap = np.array([dust.mbb_index.value])
        tempmap = np.array([dust.mbb_temperature.value])
        AI = np.array([dust.get_emission(freq_reference * u.GHz).to(u.uK_CMB,
                                                                    equivalencies=u.cmb_equivalencies(freq_reference*u.GHz)).value])

    del sky, dust
    gc.collect()
    return AI, betamap, tempmap, nu0


def dust_abi(nu):
    """Thermal dust emission. (nu en GHz)"""
    A_D = 1.36 * 10**6 * u.Jy / u.sr
    beta_D = 1.53
    T_D = 21 * u.K
    x = constants.h.to(u.J/u.GHz) * nu*u.GHz / (constants.k_B * T_D)
    return A_D * x ** beta_D * x**3 / (np.exp(x) - 1)


# general keywords
nside = 512  # resolution
npix = hp.nside2npix(nside)
instr = 'Planck'
Pathload = './maps/'


nu0_dust = 353.0  # GHz

# mask_file = "HFI_Mask_GalPlane-apo0_2048_R2.00.fits"
fskylist = np.array([0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.97, 0.99])

m20 = hp.read_map(mask_file, field=list(np.where(fskylist == 0.2)[0]))
mask20 = hp.ud_grade(m20, nside_out=nside)
m80 = hp.read_map(mask_file, field=list(np.where(fskylist == 0.8)[0]))
mask80 = hp.ud_grade(m80, nside_out=nside)


AI_dust_d1, betamap_dust_d1, tempmap_dust_d1, nu0_dust_d1 = map_pysm(
    'd1', nu0_dust)
I_nu_dust_d1 = MBBpysm(
    freq_fossil, AI_dust_d1[0, 0, :], betamap_dust_d1[0, :], tempmap_dust_d1[0, :], nu0_dust)

AI_dust_d10, betamap_dust_d10, tempmap_dust_d10, nu0_dust_d10 = map_pysm(
    'd10', nu0_dust)
I_nu_dust_d10 = MBBpysm(
    freq_fossil, AI_dust_d10[0, 0, :], betamap_dust_d10[0, :], tempmap_dust_d10[0, :], nu0_dust)

AI_dust_d12, betamap_dust_d12, tempmap_dust_d12, nu0_dust_d12 = map_pysm(
    'd12', nu0_dust)
I_nu_dust_d12 = 0
for i in range(6):
    I_nu_dust_d12 += MBBpysm(freq_fossil, AI_dust_d12[i, 0, :],
                             betamap_dust_d12[i, :], tempmap_dust_d12[i, :], nu0_dust)


# Mask
I_nu_d1_20 = I_nu_dust_d1*mask20
I_nu_d1_80 = I_nu_dust_d1*mask80

I_nu_d10_20 = I_nu_dust_d10*mask20
I_nu_d10_80 = I_nu_dust_d10*mask80

I_nu_d12_20 = I_nu_dust_d12*mask20
I_nu_d12_80 = I_nu_dust_d12*mask80

data_d1_20 = np.mean(I_nu_d1_20[:, mask20 != 0], axis=1)
data_d1_80 = np.mean(I_nu_d1_80[:, mask80 != 0], axis=1)
data_d10_20 = np.mean(I_nu_d10_20[:, mask20 != 0], axis=1)
data_d10_80 = np.mean(I_nu_d10_80[:, mask80 != 0], axis=1)
data_d12_20 = np.mean(I_nu_d12_20[:, mask20 != 0], axis=1)
data_d12_80 = np.mean(I_nu_d12_80[:, mask80 != 0], axis=1)

## PLOTS for detailed comparison ##

# plt.plot(freq_fossil, np.mean(I_nu_dust_d1, axis=1)*conv, label='d1')
# plt.plot(freq_fossil, np.mean(I_nu_dust_d10, axis=1)*conv, label='d10')
# plt.plot(freq_fossil, np.mean(I_nu_dust_d12, axis=1)*conv, label='d12')

# plt.plot(freq_fossil, data_d1_20*conv, label='d1, 20%', c='blue')
# plt.plot(freq_fossil, data_d10_20*conv, label='d10,20%', c='orange')
# plt.plot(freq_fossil, data_d12_20*conv, label='d12,20%', c='green')
# plt.plot(freq_fossil, data_d1_80*conv,
#          label='d1,80%', linestyle='--', c='blue')
# plt.plot(freq_fossil, data_d10_80*conv,
#          label='d10,80%', linestyle='--', c='orange')
# plt.plot(freq_fossil, data_d12_80*conv,
#          label='d12,80%', linestyle='--', c='green')

# plt.xlabel(r"$\nu$ [GHz]")
# plt.ylabel(r"$I_\nu$ [Jy/sr]")
# plt.legend()
# plt.loglog()
# plt.show()
