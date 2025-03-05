import numpy as np
from scipy.signal import freqz
from scipy.fftpack import next_fast_len
from .Freq2ERB  import Freq2ERB
from .Fr2Fpeak  import Fr2Fpeak


def GammaChirp(Frs, SR, OrderG=None, CoefERBw=None, CoefC=None, Phase=None, SwCarr=None, SwNorm=None):
    """
    Generate GammaChirp signals.

    INPUT:
        Frs: Asymptotic Frequency (vector)
        SR: Sampling Frequency
        OrderG: Order of Gamma function t^(OrderG-1) == n
        CoefERBw: Coefficient -> exp(-2*pi*CoefERBw*ERB(f)) == b
        CoefC: Coefficient -> exp(j*2*pi*Frs + CoefC*ln(t)) == c
        Phase: Start Phase (0 ~ 2*pi)
        SwCarr: Carrier ('cos', 'sin', 'complex', 'envelope': 3 letters)
        SwNorm: Normalization of peak spectrum level ('no', 'peak')

    OUTPUT:
        GC: GammaChirp (matrix)
        LenGC: Length of GC for each channel (vector)
        Fps: Peak Frequency (vector)
        InstFreq: Instantaneous Frequency (matrix)
    """

    if Frs is None or SR is None:
        raise ValueError("Frs and SR are required inputs.")

    Frs = np.array(Frs, ndmin=1)
    NumCh = len(Frs)

    if OrderG is None:
        OrderG = 4
    if np.isscalar(OrderG):
        OrderG = np.full(NumCh, OrderG)

    if CoefERBw is None:
        CoefERBw = 1.019
    if np.isscalar(CoefERBw):
        CoefERBw = np.full(NumCh, CoefERBw)

    if CoefC is None:
        CoefC = 0
    if np.isscalar(CoefC):
        CoefC = np.full(NumCh, CoefC)

    if Phase is None:
        Phase = 0
    if np.isscalar(Phase):
        Phase = np.full(NumCh, Phase)

    if SwCarr is None:
        SwCarr = 'cos'

    if SwNorm is None:
        SwNorm = 'no'

    ERBrate, ERBw = Freq2ERB(Frs)  # Undefined function
    LenGC1kHz = (40 * np.max(OrderG) / np.max(CoefERBw) + 200) * SR / 16000
    _, ERBw1kHz = Freq2ERB(1000)  # Undefined function

    if SwCarr.startswith('sin'):
        Phase -= np.pi / 2

    Phase += CoefC * np.log(Frs / 1000)
    LenGC = np.fix(LenGC1kHz * ERBw1kHz / ERBw).astype(int)

    GC = np.zeros((NumCh, np.max(LenGC)))
    Fps = None
    InstFreq = None

    if 'Fps' in locals():
        Fps = Fr2Fpeak(OrderG, CoefERBw, CoefC, Frs)  # Undefined function

    if 'InstFreq' in locals():
        InstFreq = np.zeros((NumCh, np.max(LenGC)))

    for nch in range(NumCh):
        t = np.arange(1, LenGC[nch]) / SR
        GammaEnv = t**(OrderG[nch] - 1) * np.exp(-2 * np.pi * CoefERBw[nch] * ERBw[nch] * t)
        GammaEnv = np.concatenate(([0], GammaEnv / np.max(GammaEnv)))

        if SwCarr.startswith('env'):
            Carrier = np.ones_like(GammaEnv)
        elif SwCarr.startswith('com'):
            Carrier = np.concatenate(([0], np.exp(1j * (2 * np.pi * Frs[nch] * t + CoefC[nch] * np.log(t) + Phase[nch]))))
        else:
            Carrier = np.concatenate(([0], np.cos(2 * np.pi * Frs[nch] * t + CoefC[nch] * np.log(t) + Phase[nch])))

        GC[nch, :LenGC[nch]] = GammaEnv * Carrier

        if InstFreq is not None:
            InstFreq[nch, :LenGC[nch]] = np.concatenate(([0], Frs[nch] + CoefC[nch] / (2 * np.pi * t)))

        if SwNorm == 'peak':
            frsp, freq = freqz(GC[nch, :LenGC[nch]], 1, next_fast_len(LenGC[nch]), SR)
            fp = Fr2Fpeak(OrderG[nch], CoefERBw[nch], CoefC[nch], Frs[nch])  # Undefined function
            np_idx = np.argmin(np.abs(freq - fp))
            GC[nch, :] /= np.abs(frsp[np_idx])

    return GC, LenGC, Fps, InstFreq

