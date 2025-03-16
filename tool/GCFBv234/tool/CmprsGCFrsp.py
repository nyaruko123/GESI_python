import numpy as np
from .Fr2Fpeak import Fr2Fpeak
from .GammaChirpFrsp import GammaChirpFrsp
from .AsymCmpFrspV2 import AsymCmpFrspV2



def CmprsGCFrsp(Fr1, fs=48000, n=4, b1=1.81, c1=-2.96, frat=1, b2=2.17, c2=2.20, NfrqRsl=1024):
    """
    Compute compressive GammaChirp frequency response.

    INPUT:
        Fr1: Resonance Frequency (vector)
        fs: Sampling Frequency
        n: Order of Gamma function t^(n-1) (vector)
        b1: b1 for exp(-2*pi*b1*ERB(f))
        c1: c1 for exp(j*2*pi*Fr + c1*ln(t))
        frat: frequency ratio. Fr2 = frat*Fp1
        b2: b2 for HP-AF
        c2: c2 for HP-AF
        NfrqRsl: frequency resolution

    OUTPUT:
        cGCresp: dict for cGC response
            - pGCFrsp: passive gc frequency response (NumCh*NfrqRsl matrix)
            - cGCFrsp: compressive gc frequency response (NumCh*NfrqRsl matrix)
            - cGCNrmFrsp: Normalized cGCFrsp (NumCh*NfrqRsl matrix)
            - ACFrsp: Asym Compensation Filter frequency response
            - AsymFunc: Asym Function
            - freq: frequency (1 * NfrqRsl vector)
            - Fp2: peak frequency
            - ValFp2: peak value
    """

    Fr1 = np.array(Fr1, ndmin=1)
    NumCh = len(Fr1)

    if np.isscalar(n):
        n = np.full(NumCh, n)
    if np.isscalar(b1):
        b1 = np.full(NumCh, b1)
    if np.isscalar(c1):
        c1 = np.full(NumCh, c1)
    if np.isscalar(frat):
        frat = np.full(NumCh, frat)
    if np.isscalar(b2):
        b2 = np.full(NumCh, b2)
    if np.isscalar(c2):
        c2 = np.full(NumCh, c2)

    #pGCFrsp, freq = GammaChirpFrsp(Fr1, fs, n, b1, c1, 0, NfrqRsl)  # Undefined function\
    pGCFrsp, freq, _, _, _ = GammaChirpFrsp(Fr1, fs, n, b1, c1, 0, NfrqRsl)

    Fp1 = Fr2Fpeak(n, b1, c1, Fr1)  # Undefined function
    Fr2 = frat * Fp1
    ACFFrsp, freq, AsymFunc = AsymCmpFrspV2(Fr2, fs, b2, c2, NfrqRsl)  # Undefined function
    cGCFrsp = pGCFrsp * AsymFunc
    ValFp2 = np.max(cGCFrsp, axis=1)
    NormFactFp2 = 1.0 / ValFp2

    cGCresp = {
        'Fr1': Fr1,
        'n': n,
        'b1': b1,
        'c1': c1,
        'frat': frat,
        'b2': b2,
        'c2': c2,
        'NfrqRsl': NfrqRsl,
        'pGCFrsp': pGCFrsp,
        'cGCFrsp': cGCFrsp,
        'cGCNrmFrsp': cGCFrsp * (NormFactFp2[:, np.newaxis]),
        'ACFFrsp': ACFFrsp,
        'AsymFunc': AsymFunc,
        'Fp1': Fp1,
        'Fr2': Fr2,
        'Fp2': freq[np.argmax(cGCFrsp, axis=1)],
        'ValFp2': ValFp2,
        'NormFctFp2': NormFactFp2,
        'freq': freq
    }

    return cGCresp

