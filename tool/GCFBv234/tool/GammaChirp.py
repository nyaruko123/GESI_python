import numpy as np
from scipy.signal import freqz
from scipy.fftpack import next_fast_len
from .Freq2ERB import Freq2ERB
from .Fr2Fpeak import Fr2Fpeak

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

    Frs = np.asarray(Frs, dtype=np.float64).flatten()
    NumCh = len(Frs)

    if OrderG is None:
        OrderG = np.full(NumCh, 4, dtype=np.float64)
    else:
        OrderG = np.asarray(OrderG, dtype=np.float64).flatten()
    
    if CoefERBw is None:
        CoefERBw = np.full(NumCh, 1.019, dtype=np.float64)
    else:
        CoefERBw = np.asarray(CoefERBw, dtype=np.float64).flatten()
    
    if CoefC is None:
        CoefC = np.zeros(NumCh, dtype=np.float64)
    else:
        CoefC = np.asarray(CoefC, dtype=np.float64).flatten()
    
    # if Phase is None:
    #     Phase = np.zeros(NumCh, dtype=np.float64)
    # else:
    #     Phase = np.asarray(Phase, dtype=np.float64).flatten()

    if Phase is None:
        Phase = 0  # 对应 MATLAB 代码中的 `if length(Phase) == 0, Phase = 0;`   
    if isinstance(Phase, (int, float)):
        Phase = np.full(NumCh, Phase, dtype=np.float64)  # 对应 `Phase = Phase*ones(NumCh,1);`
    else:
        Phase = np.asarray(Phase, dtype=np.float64).flatten()  # 确保是 NumPy 数组
        
    
    if SwCarr is None:
        SwCarr = 'cos'
    
    if SwNorm is None:
        SwNorm = 'no'
    
    ERBrate, ERBw = Freq2ERB(Frs)
    LenGC1kHz = (40 * np.max(OrderG) / np.max(CoefERBw) + 200) * SR / 16000
    _, ERBw1kHz = Freq2ERB(1000)
    
    if SwCarr.startswith('sin'):
        Phase -= np.pi / 2
    
    Phase = np.asarray(Phase, dtype=np.float64)  # Ensure Phase is float64
    Frs = np.asarray(Frs, dtype=np.float64)  # Ensure Frs is float64
    Phase += CoefC * np.log(Frs / 1000)
    
    LenGC = np.fix(LenGC1kHz * ERBw1kHz / ERBw).astype(int)
    GC = np.zeros((NumCh, np.max(LenGC)))
    Fps = None
    InstFreq = None
    
    if 'Fps' in locals():
        Fps = Fr2Fpeak(OrderG, CoefERBw, CoefC, Frs)
    
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
            fft_size = next_fast_len(LenGC[nch])
            frsp, freq = freqz(GC[nch, :LenGC[nch]], 1, fft_size, SR)
            fp = Fr2Fpeak(OrderG[nch], CoefERBw[nch], CoefC[nch], Frs[nch])
            np_idx = np.argmin(np.abs(freq - fp))
            np_idx = min(np_idx, len(frsp) - 1)  # 确保索引不会超界
            GC[nch, :] /= np.abs(frsp[np_idx])
    
    return GC, LenGC, Fps, InstFreq
