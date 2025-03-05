"""
Ideal observer of speech intelligibility used in GEDI (sEPSM)
based on listeners' vocabulary size. See also sEPSM.
Irino, T.
Created: 30 Jan 2022 IT, from mrGEDI
Modified: 30 Jan 2022 IT

INPUT:
    Metric: Metric derived in the main. SDRenv_lin when using GEDI
    Param: Vector of [k, q, m, sigma_s]

OUTPUT:
    Pcorrect: Percent correct of speech intelligibility

Note in GEDI:
    IdealObserver: Converts the overall SDRenv to percent correct.

    Usage: Pcorrect = IdealObserver(SDRenv_lin, parameters)
    Parameters: Vector with the parameters for the ideal Observer formatted as [k q m sigma_s]
    # d_prime = k * (SDRenv_lin) ** q;

    Green, D. M. and Birdsall, T. G. (1964). "The effect of vocabulary size",
    In Signal Detection and Recognition by Human Observers,
    edited by John A. Swets (John Wiley & Sons, New York)
"""

import numpy as np
from scipy.special import erfinv, erf

def IdealObserver_VocabGauss(Metric, Param):
    if len(Param) < 4:
        raise ValueError('You have to specify the k, q, m, sigma_s parameters for the IdealObserver')

    k = Param[0]
    q = Param[1]
    m = Param[2]
    sigma_s = Param[3]

    # Converting from Metric (SNRenv in GEDI) to d_prime
    d_prime = k * (Metric) ** q

    # Converting from d_prime to Percent correct, Green and Birdsall (1964)
    Un = 1 * norminv_erf(1 - (1 / m))
    mn = Un + (0.577 / Un)
    sig_n = 1.28255 / Un
    Pcorrect = normcdf_erf(d_prime, mn, np.sqrt(sigma_s ** 2 + sig_n ** 2)) * 100

    return Pcorrect

def norminv_erf(p):
    """Equivalent to norminv in Statistics Toolbox."""
    return np.sqrt(2) * erfinv(2 * p - 1)

def normcdf_erf(x, mu=0, sigma=1):
    """Equivalent to normcdf in Statistics Toolbox."""
    return (1 + erf((x - mu) / (sigma * np.sqrt(2)))) / 2
