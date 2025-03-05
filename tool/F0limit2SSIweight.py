"""

    SSI weight function from SSI weight  --- Simplified version
    Toshio IRINO
    Created:    5 Feb 2022   from CalSSIweight
    Modified:   5 Feb 2022
    Modified:   1 Sep 2022   Add NaN or negative handling

    Function: [SSIweight, SSIparam] = F0limit2SSIweight(SSIparam)
    INPUT:  SSIparam.Fr1 :  Filterbank channel frequency (== GCparam.Fr1==GCresp.Fr1)
            SSIparam.h_max = 5;
            SSIparam.F0_limit =  F0  # specified from adaptive F0 value
    OUTPUT: SSIweight, SSIparam

"""
import numpy as np


def F0limit2SSIweight(SSIparam):
    if 'Fr1' not in SSIparam:
        raise ValueError('SSIparam.Fr1 (== GCparam.Fr1) is essential.')

    if 'F0_limit' not in SSIparam:
        SSIparam['F0_limit'] = 150  # default

    if 'h_max' not in SSIparam:
        SSIparam['h_max'] = 5  # default

    if SSIparam['F0_limit'] > 0:  # when F0 is positive
        SSIparam['TI_limit'] = 1 / SSIparam['F0_limit']  # Limit of successive Glottal pulse
        SSIweight = np.minimum(SSIparam['Fr1'] * SSIparam['TI_limit'], SSIparam['h_max']) / SSIparam['h_max']
    elif SSIparam['F0_limit'] == 0:  # when F0 == 0, Uniform
        SSIparam['TI_limit'] = np.inf
        SSIweight = np.ones_like(SSIparam['Fr1'])
    elif np.isnan(SSIparam['F0_limit']):  # NaN
        raise ValueError('SSIparam.F0_limit should not be NaN.')
    else:
        raise ValueError('SSIparam.F0_limit should not be negative.')

    return SSIweight, SSIparam

