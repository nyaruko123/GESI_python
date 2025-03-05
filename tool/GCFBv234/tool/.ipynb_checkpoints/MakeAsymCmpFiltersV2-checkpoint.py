import numpy as np
from .Freq2ERB import Freq2ERB

def MakeAsymCmpFiltersV2(fs, Frs, b, c):
    if Frs.ndim != 1:
        raise ValueError('Frs should be a column vector Frs(:).')

    NumCh = len(Frs)
    ERBw = Freq2ERB(Frs)

    ACFcoef = {
        'fs': fs,
        'ap': np.zeros((NumCh, 3, 4)),
        'bz': np.zeros((NumCh, 3, 4))
    }

    NumFilt = 4
    p0 = 2
    p1 = 1.7818 * (1 - 0.0791 * b) * (1 - 0.1655 * np.abs(c))
    p2 = 0.5689 * (1 - 0.1620 * b) * (1 - 0.0857 * np.abs(c))
    p3 = 0.2523 * (1 - 0.0244 * b) * (1 + 0.0574 * np.abs(c))
    p4 = 1.0724

    for Nfilt in range(1, NumFilt + 1):
        r = np.exp(-p1 * (p0 / p4) ** (Nfilt - 1) * 2 * np.pi * b * ERBw / fs)
        delFrs = (p0 * p4) ** (Nfilt - 1) * p2 * c * b * ERBw
        phi = 2 * np.pi * np.maximum(Frs + delFrs, 0) / fs
        psy = 2 * np.pi * np.maximum(Frs - delFrs, 0) / fs
        fn = Frs

        ap = np.column_stack((np.ones_like(r), -2 * r * np.cos(phi), r ** 2))
        bz = np.column_stack((np.ones_like(r), -2 * r * np.cos(psy), r ** 2))

        vwr = np.exp(1j * 2 * np.pi * fn / fs)
        vwrs = np.column_stack((np.ones_like(vwr), vwr, vwr ** 2))
        nrm = np.abs(np.sum(vwrs * ap, axis=1) / np.sum(vwrs * bz, axis=1))
        bz *= nrm[:, np.newaxis]

        ACFcoef['ap'][:, :, Nfilt - 1] = ap
        ACFcoef['bz'][:, :, Nfilt - 1] = bz

    ACFcoefConv = {
        'ap': np.zeros((NumCh, 9)),
        'bz': np.zeros((NumCh, 9))
    }

    for nch in range(NumCh):
        ap1 = np.array([1])
        bz1 = np.array([1])
        for Nfilt in range(4):
            ap1 = np.convolve(ACFcoef['ap'][nch, :, Nfilt], ap1)
            bz1 = np.convolve(ACFcoef['bz'][nch, :, Nfilt], bz1)
        ACFcoefConv['ap'][nch, :9] = ap1
        ACFcoefConv['bz'][nch, :9] = bz1

    return ACFcoef, ACFcoefConv


