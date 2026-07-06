import numpy as np
from astropy import constants
from bruit_fossil import freq_fossil


def uK_RJ_to_MJy_sr(nu):
    """
    Compute factors to convert brightness from uK_RJ to MJy/sr for all input frequencies.

    Parameters
    ----------
    nu : float or np.array
        Frequencies for which conversion factors should be computed in GHz.

    Returns
    -------
    float or np.array
        Conversion factors for all input frequencies.

    """
    k = constants.k_B.cgs.value
    c = constants.c.cgs.value

    return 2*k * (nu*1e9 / c)**2 * 1e11


def uK_RJ_to_uK_CMB(nu):
    """
    Compute factors to convert brightness from uK_RJ to uK_CMB for all input frequencies.

    Parameters
    ----------
    nu : float or np.array
        Frequencies for which conversion factors should be computed in GHz.

    Returns
    -------
    float or np.array
        Conversion factors for all input frequencies.

    """
    h = constants.h.cgs.value
    k = constants.k_B.cgs.value
    Tcmb = 2.725

    x = h*nu*1e9 / (k*Tcmb)
    return (np.exp(x) - 1)**2 / (x**2 * np.exp(x))


def unit_conversion(nu, input_unit, output_unit):
    """
    Compute factors to convert brightness from input_unit to output_unit for all input frequencies.

    Parameters
    ----------
    nu : float or np.array
        Frequencies for which conversion factors should be computed in GHz.
    input_unit : string
        Unit from which conversion factors should be computed. Can be 'MJy/sr', 'uK_RJ' or 'uK_CMB'.
    output_unit : string
        Unit to which conversion factors should be computed. Can be 'MJy/sr', 'uK_RJ' or 'uK_CMB'.

    Returns
    -------
    float or np.array
        Conversion factors for all input frequencies.

    """
    if input_unit == 'uK_CMB':
        if output_unit == 'uK_CMB':
            return np.ones_like(nu)

        elif output_unit == 'uK_RJ':
            return 1 / uK_RJ_to_uK_CMB(nu)

        elif output_unit == 'MJy/sr':
            return uK_RJ_to_MJy_sr(nu) / uK_RJ_to_uK_CMB(nu)

        else:
            raise ValueError('Incorrect output unit')

    elif input_unit == 'uK_RJ':
        if output_unit == 'uK_CMB':
            return uK_RJ_to_uK_CMB(nu)

        elif output_unit == 'uK_RJ':
            return np.ones_like(nu)

        elif output_unit == 'MJy/sr':
            return uK_RJ_to_MJy_sr(nu)

        else:
            raise ValueError('Incorrect output unit')

    elif input_unit == 'MJy/sr':
        if output_unit == 'uK_CMB':
            return uK_RJ_to_uK_CMB(nu) / uK_RJ_to_MJy_sr(nu)

        elif output_unit == 'uK_RJ':
            return 1 / uK_RJ_to_MJy_sr(nu)

        elif output_unit == 'MJy/sr':
            return np.ones_like(nu)

        else:
            raise ValueError('Incorrect output unit')

    else:
        raise ValueError('Incorrect input unit')


conv = 10**6*unit_conversion(freq_fossil, 'uK_CMB', 'MJy/sr')
