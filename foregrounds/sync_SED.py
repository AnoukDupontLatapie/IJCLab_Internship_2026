from bruit_fossil import freq_fossil
from conversion import conv
from func_theo import PowerLawpysm
from mask_file_download import mask_file
import numpy as np
import healpy as hp
import pysm3.units as u
import pysm3
import matplotlib.pyplot as plt
from astropy import constants
import gc
plt.rcParams['font.family'] = 'sans-serif'


def map_pysm(fgtype, freq_reference): #Generate map at instrument frequencies

    if fgtype[0]=='d':
        model_type='mbb'
    elif fgtype[0]=='s':
        model_type='powerlaw'

    sky = pysm3.Sky(nside=nside, preset_strings=['%s'%fgtype])#,'s%s'%synctype])
    
    if model_type=='powerlaw':
        if fgtype=="s1" or fgtype=="s5":
            sync = sky.components[0]
            betamap= np.array([sync.pl_index.value])
            nu0= sync.freq_ref_I.value
            AI = np.array([sync.get_emission(freq_reference * u.GHz).to(u.uK_CMB, equivalencies=u.cmb_equivalencies(freq_reference*u.GHz)).value])
            return AI, betamap, nu0
        
        elif fgtype=="s7":
            sync = sky.components[0]
            betamap= np.array([sync.pl_index.value])
            gammamap = np.array([sync.spectral_curvature.value])
            nu0= sync.freq_ref_I.value
            AI = np.array([sync.get_emission(freq_reference * u.GHz).to(u.uK_CMB, equivalencies=u.cmb_equivalencies(freq_reference*u.GHz)).value])        
            return AI, betamap, gammamap, nu0


def synchrotron_abi(nu):
    """Synchrotron emission."""
    A_S = 288.0 * u.Jy / u.sr
    alpha_S = -0.82
    w_S = 0.2
    nu_0 = 100 * u.GHz
    return A_S * (nu*u.GHz / nu_0) ** alpha_S * (1 + w_S/2 * np.log(nu*u.GHz/nu_0)**2)


# general keywords
nside = 512  # resolution
npix = hp.nside2npix(nside)
instr = 'Planck'
Pathload = './maps/'


nu0_sync = 23.0  # GHz

mask_file = "HFI_Mask_GalPlane-apo0_2048_R2.00.fits"
fskylist = np.array([0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.97, 0.99])

m20 = hp.read_map(mask_file, field=list(np.where(fskylist == 0.2)[0]))
mask20 = hp.ud_grade(m20, nside_out=nside)
m80 = hp.read_map(mask_file, field=list(np.where(fskylist == 0.8)[0]))
mask80 = hp.ud_grade(m80, nside_out=nside)

AI_sync_s1, betamap_sync_s1, nu0_sync_s1 = map_pysm('s1', nu0_sync)
I_nu_sync_s1 = PowerLawpysm(freq_fossil, AI_sync_s1[0,0,:], betamap_sync_s1[0,:], None, nu0_sync)

AI_sync_s5, betamap_sync_s5, nu0_sync_s5 = map_pysm('s5', nu0_sync)
I_nu_sync_s5 = PowerLawpysm(freq_fossil, AI_sync_s5[0,0,:], betamap_sync_s5[0,:], None, nu0_sync)

AI_sync_s7, betamap_sync_s7, gammamap_sync_s7, nu0_sync_s7 = map_pysm('s7', nu0_sync)
I_nu_sync_s7 = PowerLawpysm(freq_fossil, AI_sync_s7[0,0,:], betamap_sync_s7[0,:], gammamap_sync_s7[0,:], nu0_sync)


# Mask
I_nu_s1_20 = I_nu_sync_s1*mask20
I_nu_s1_80 = I_nu_sync_s1*mask80

I_nu_s5_20 = I_nu_sync_s5*mask20
I_nu_s5_80 = I_nu_sync_s5*mask80

I_nu_s7_20 = I_nu_sync_s7*mask20
I_nu_s7_80 = I_nu_sync_s7*mask80

data_s1_20= np.mean(I_nu_s1_20[:,mask20!=0],axis=1)
data_s1_80= np.mean(I_nu_s1_80[:,mask80!=0],axis=1)
data_s5_20= np.mean(I_nu_s5_20[:,mask20!=0],axis=1)
data_s5_80= np.mean(I_nu_s5_80[:,mask80!=0],axis=1)
data_s7_20= np.mean(I_nu_s7_20[:,mask20!=0],axis=1)
data_s7_80= np.mean(I_nu_s7_80[:,mask80!=0],axis=1)


plt.plot(freq_fossil, np.mean(I_nu_sync_s1, axis=1)*conv, label='Modèle s1')
plt.plot(freq_fossil, np.mean(I_nu_sync_s5, axis=1)*conv, label='Modèle s5')
plt.plot(freq_fossil, np.mean(I_nu_sync_s7, axis=1)*conv, label='Modèle s7')

plt.plot(freq_fossil, data_s1_20*conv, label='s1, 20%', c='blue')
plt.plot(freq_fossil, data_s5_20*conv, label='s5,20%', c='orange')
plt.plot(freq_fossil, data_s7_20*conv, label='s7,20%', c='green')
plt.plot(freq_fossil, data_s1_80*conv,
         label='s1,80%', linestyle='--', c='blue')
plt.plot(freq_fossil, data_s5_80*conv,
         label='s5,80%', linestyle='--', c='orange')
plt.plot(freq_fossil, data_s7_80*conv,
         label='s7,80%', linestyle='--', c='green')

plt.xlabel(r"$\nu$ [GHz]")
plt.ylabel(r"$I_\nu$ [Jy/sr]")
plt.legend()
plt.loglog()
plt.show()
