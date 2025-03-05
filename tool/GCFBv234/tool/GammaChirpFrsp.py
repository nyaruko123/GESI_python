import numpy as np
from .Freq2ERB import Freq2ERB

def GammaChirpFrsp(Frs, SR=48000, OrderG=4, CoefERBw=1.019, CoefC=0, Phase=None, NfrqRsl=1024):
    """
    Compute the frequency response of a GammaChirp filter.

    INPUT:
        Frs: Resonance Frequency (vector)
        SR: Sampling Frequency
        OrderG: Order of Gamma function t^(OrderG-1) (vector)
        CoefERBw: Coefficient -> exp(-2*pi*CoefERBw*ERB(f))
        CoefC: Coefficient -> exp(j*2*pi*Fr + CoefC*ln(t))
        Phase: Start Phase (0 ~ 2*pi)
        NfrqRsl: frequency resolution

    OUTPUT:
        AmpFrsp: abs(Response) (NumCh*NfrqRsl matrix)
        freq: frequency (1 * NfrqRsl vector)
        Fpeak: Peak frequency (NumCh * 1 vector)
        GrpDly: Group delay (NumCh*NfrqRsl matrix)
        PhsFrsp: angle(Response) (NumCh*NfrqRsl matrix)
    """

    if Frs is None:
        raise ValueError("Frs is a required input.")
    if SR is None or SR == 0:
        raise ValueError("Specify Sampling Frequency.")

    Frs = np.array(Frs, ndmin=1)
    NumCh = len(Frs)

    if np.isscalar(OrderG):
        OrderG = np.full(NumCh, OrderG)
    if np.isscalar(CoefERBw):
        CoefERBw = np.full(NumCh, CoefERBw)
    if np.isscalar(CoefC):
        CoefC = np.full(NumCh, CoefC)
    if Phase is None:
        Phase = np.zeros(NumCh)
    if NfrqRsl < 256:
        raise ValueError("NfrqRsl must be >= 256.")

    ERBrate, ERBw = Freq2ERB(Frs)  # Undefined function
    freq = np.arange(NfrqRsl) / NfrqRsl * SR / 2

    one1 = np.ones(NfrqRsl)
    bh = (CoefERBw * ERBw)[:, np.newaxis] * one1
    fd = (np.ones((NumCh, 1)) * freq) - Frs[:, np.newaxis]
    cn = (CoefC / OrderG)[:, np.newaxis] * one1
    n = OrderG[:, np.newaxis] * one1
    c = CoefC[:, np.newaxis] * one1
    Phase = Phase[:, np.newaxis] * one1

    AmpFrsp = ((1 + cn**2) / (1 + (fd / bh)**2))**(n / 2) * np.exp(c * (np.arctan(fd / bh) - np.arctan(cn)))

    Fpeak = Frs + CoefERBw * ERBw * CoefC / OrderG
    GrpDly = (1 / (2 * np.pi)) * (n * bh + c * fd) / (bh**2 + fd**2)
    PhsFrsp = -n * np.arctan(fd / bh) - c / 2 * np.log((2 * np.pi * bh)**2 + (2 * np.pi * fd)**2) + Phase

    return AmpFrsp, freq, Fpeak, GrpDly, PhsFrsp

