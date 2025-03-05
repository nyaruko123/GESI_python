import numpy as np
from .Freq2ERB import Freq2ERB
from .MakeAsymCmpFiltersV2 import MakeAsymCmpFiltersV2


def AsymCmpFrspV2(Frs, fs, b, c, NfrqRsl=1024, NumFilt=4):
    """
    Compute the frequency response of an Asymmetric Compensation Filter.

    INPUT:
        fs: Sampling frequency
        Frs: array of the center frequencies
        b: array or scalar of a bandwidth coefficient
        c: array or scalar of asymmetric parameters
        NfrqRsl: >64 for freq. resolution for linear frequency scale
                 0 for specify response at Frs
        NumFilt: Number of 2nd-order filters, default 4

    OUTPUT:
        ACFFrsp: abs(Frsp of ACF) (NumCh*NfrqRsl matrix)
        freq: frequency (1*NfrqRsl vector)
        AsymFunc: Original Asymmetric Function (NumCh*NfrqRsl matrix)
    """

    if NfrqRsl is None or NfrqRsl == 0:
        NfrqRsl = 1024
    if NumFilt is None:
        NumFilt = 4
    if NumFilt != 4:
        raise ValueError('NumFilter should be 4.')

    Frs = np.array(Frs, ndmin=1)
    b = np.array(b, ndmin=1)
    c = np.array(c, ndmin=1)
    NumCh = len(Frs)

    if NfrqRsl >= 64:
        freq = np.arange(NfrqRsl) / NfrqRsl * fs / 2
    elif NfrqRsl == 0:
        freq = Frs
        NfrqRsl = len(freq)
    else:
        raise ValueError('Specify NfrqRsl 0 for Frs or N>=64 for linear-freq scale')

    SwCoef = 0  # self consistency

    if SwCoef == 0:
        p0 = 2
        p1 = 1.7818 * (1 - 0.0791 * b) * (1 - 0.1655 * np.abs(c))
        p2 = 0.5689 * (1 - 0.1620 * b) * (1 - 0.0857 * np.abs(c))
        p3 = 0.2523 * (1 - 0.0244 * b) * (1 + 0.0574 * np.abs(c))
        p4 = 1.0724
    else:
        ACFcoef = MakeAsymCmpFiltersV2(fs, Frs, b, c)  # Undefined function

    _, ERBw = Freq2ERB(Frs)  # Undefined function
    ACFFrsp = np.ones((NumCh, NfrqRsl))
    freq2 = np.hstack((np.ones((NumCh, 1)) * freq, Frs[:, np.newaxis]))

    for Nfilt in range(1, NumFilt + 1):
        if SwCoef == 0:
            r = np.exp(-p1 * (p0 / p4)**(Nfilt - 1) * 2 * np.pi * b * ERBw / fs)
            delfr = (p0 * p4)**(Nfilt - 1) * p2 * c * b * ERBw
            phi = 2 * np.pi * np.maximum(Frs + delfr, 0) / fs
            psy = 2 * np.pi * np.maximum(Frs - delfr, 0) / fs
            ap = np.column_stack((np.ones(NumCh), -2 * r * np.cos(phi), r**2))
            bz = np.column_stack((np.ones(NumCh), -2 * r * np.cos(psy), r**2))
        else:
            ap = ACFcoef.ap[:, :, Nfilt - 1]
            bz = ACFcoef.bz[:, :, Nfilt - 1]

        cs1 = np.cos(2 * np.pi * freq2 / fs)
        cs2 = np.cos(4 * np.pi * freq2 / fs)
        bzz0 = (bz[:, 0]**2 + bz[:, 1]**2 + bz[:, 2]**2)[:, np.newaxis]
        bzz1 = (2 * bz[:, 1] * (bz[:, 0] + bz[:, 2]))[:, np.newaxis]
        bzz2 = (2 * bz[:, 0] * bz[:, 2])[:, np.newaxis]
        hb = bzz0 + bzz1 * cs1 + bzz2 * cs2

        app0 = (ap[:, 0]**2 + ap[:, 1]**2 + ap[:, 2]**2)[:, np.newaxis]
        app1 = (2 * ap[:, 1] * (ap[:, 0] + ap[:, 2]))[:, np.newaxis]
        app2 = (2 * ap[:, 0] * ap[:, 2])[:, np.newaxis]
        ha = app0 + app1 * cs1 + app2 * cs2

        H = np.sqrt(hb / ha)
        Hnorm = H[:, NfrqRsl]  # Normalization by fn value

        ACFFrsp *= H[:, :NfrqRsl] / Hnorm[:, np.newaxis]

    fd = np.outer(np.ones(NumCh), freq) - Frs[:, np.newaxis]
    be = np.outer(b * ERBw, np.ones(NfrqRsl))
    cc = np.outer(c, np.ones(NfrqRsl))
    AsymFunc = np.exp(cc * np.arctan2(fd, be))

    return ACFFrsp, freq, AsymFunc

