"""
Alignment of signals using cross-correlation
    --> It seems all right within the range of SpIntel SNR
        in the case of additive noise
    --> This may be replaced by a more sophisticated method
        (e.g., using GCFB representation)
"""

import numpy as np

def TimeAlignXcorr(SndTestIn, SndRefIn):
    print('Time Alignment using Xcorr. Equalizing length(SndTest) to length(SndRef).')

    if len(SndTestIn) < len(SndRefIn):
        raise ValueError('SndTest should be longer than SndRef.')

    LenSnd = len(SndRefIn)

    XcorrSnd = np.correlate(SndTestIn, SndRefIn, mode='full')
    Lag = np.arange(-len(SndRefIn) + 1, len(SndTestIn))
    IndxMaxLag = np.argmax(np.abs(XcorrSnd))
    NumTimeLag = Lag[IndxMaxLag]

    if NumTimeLag >= 0:
        SndTest1 = np.concatenate((SndTestIn, np.zeros(NumTimeLag)))
        SndTestOut = SndTest1[NumTimeLag:NumTimeLag + LenSnd]
    else:
        SndTest1 = np.concatenate((np.zeros(abs(NumTimeLag)), SndTestIn))
        SndTestOut = SndTest1[:LenSnd]

    ParamTA = {
        'XcorrSnd': XcorrSnd,
        'Lag': Lag,
        'IndxMaxLag': IndxMaxLag,
        'NumTimeLag': NumTimeLag
    }

    print(f'Time Lag in sample point = {NumTimeLag}')

    return SndTestOut, ParamTA
