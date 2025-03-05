"""
Function [Frs, WFvals] = EqualFreqScale(NameScale, NumCh, RangeFreq)
    INPUT:
        NameScale: 'ERB', 'mel', 'log', 'linear'
        NumCh: Number of channels
        RangeFreq: Frequency Range
    OUTPUT:
        Frs: Frequency vector
        WFvals: Warped Frequency values
"""

import numpy as np
import librosa
from .Freq2ERB import Freq2ERB
from .ERB2Freq import ERB2Freq

def EqualFreqScale(NameScale, NumCh, RangeFreq):

    if len(RangeFreq) < 2 or RangeFreq[0] >= RangeFreq[1]:
        raise ValueError('RangeFreq(1) should be less than RangeFreq(2).')

    # **确保 RangeFreq 是 numpy 数组**
    RangeFreq = np.asarray(RangeFreq, dtype=np.float64)

    if NameScale.lower() == 'linear':
        dWF = (RangeFreq[1] - RangeFreq[0]) / (NumCh - 1)
        WFvals = np.arange(RangeFreq[0], RangeFreq[1] + np.finfo(float).eps * 1000, dWF)
        Frs = WFvals

    elif NameScale.lower() == 'mel':
        RangeWF = librosa.hz_to_mel(RangeFreq)
        dWF = (RangeWF[1] - RangeWF[0]) / (NumCh - 1)
        WFvals = np.arange(RangeWF[0], RangeWF[1] + np.finfo(float).eps * 1000, dWF)
        Frs = librosa.mel_to_hz(WFvals)

    elif NameScale.lower() == 'erb':
        ERBrate, _ = Freq2ERB(RangeFreq)  # 只取 ERBrate
        dWF = (ERBrate[1] - ERBrate[0]) / (NumCh - 1)
        WFvals = np.arange(ERBrate[0], ERBrate[1] + np.finfo(float).eps * 1000, dWF)
        Frs = ERB2Freq(WFvals)

    elif NameScale.lower() == 'log':
        if min(RangeFreq) < 50:
            print('min(RangeFreq) < 50. Replaced by 50.')
            RangeFreq[0] = 50
        RangeWF = np.log10(RangeFreq)
        dWF = (RangeWF[1] - RangeWF[0]) / (NumCh - 1)
        WFvals = np.arange(RangeWF[0], RangeWF[1] + np.finfo(float).eps * 1000, dWF)
        Frs = 10 ** WFvals

    else:
        raise ValueError('Specify NameScale correctly')

    return Frs, WFvals
