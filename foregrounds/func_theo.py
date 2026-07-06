import numpy as np
import pysm3.units as u
import matplotlib.pyplot as plt
from astropy import constants
plt.rcParams['font.family'] = 'sans-serif'


def B(nu, T):
    """Planck function.

    :param nu: (float) frequency in GHz at which to evaluate planck function.
    :param T: (float) temperature of black body.
    :return: float -- black body brightness.

    """
    x = constants.h.to(u.J / u.GHz)*nu*u.GHz/(constants.k_B * T*u.K)
    return 2.*constants.h * (nu*u.GHz)**3 / (constants.c)**2/np.expm1(x)


def mbb(nu, beta, T):
    return B(nu, T)*nu**beta


def MBBpysm(freq, A, beta, T, nu0):
    """
    Generate MBB to reproduce pysm maps
    A in muKCMB
    """
    factor = u.K_RJ.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(
        freq*u.GHz))/u.K_RJ.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(nu0*u.GHz))
    mapd = np.array([A*mbb(freq[f], beta-2, T)/mbb(nu0, beta-2, T)
                    * factor[f] for f in range(len(freq))])
    return mapd


def temp_shift(x):
    return x*np.exp(x)/(np.exp(x)-1)**2


def y_dist(nu):
    T0 = 2.7255 * u.K
    x = (constants.h * nu*u.GHz / (constants.k_B * T0)
         ).to(u.dimensionless_unscaled)
    I0 = (2 * constants.h * (nu*u.GHz)**3 / constants.c**2) / u.sr
    y0 = 1.77e-6
    result = I0 * y0 * temp_shift(x) * (x*(np.exp(x)+1)/(np.exp(x)-1) - 4)
    return result.to(u.Jy / u.sr)


def mu_dist(nu):
    T0 = 2.7255 * u.K
    x = (constants.h * nu*u.GHz / (constants.k_B * T0)
         ).to(u.dimensionless_unscaled)
    I0 = (2 * constants.h * (nu*u.GHz)**3 / constants.c**2) / u.sr
    mu0 = 2e-8
    result = I0 * mu0 * temp_shift(x) * (0.4561 - 1/x)
    return result.to(u.Jy / u.sr)
