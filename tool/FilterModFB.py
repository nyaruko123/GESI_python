"""
function [OutMFB, ParamMFB] = FilterModFB(Env, ParamMFB)
INPUT:
        Env:  The envelope to be filtered
        ParamMFB.fs: sampling frequency of the envelope
        ParamMFB.fc: center frequencies of the modulation filters
                              default: [1 2 4 8 16 32 64 128 256]; --> [1 2 4 8 16 32 64 128 256 512];
        ParamMFB.SwPlot:  Plot frequency response of MFB

OUTPUT:
        OutMFB:  Temporal outputs for each of the modulation filters
        ParamMFB: Parameter

See:
function x_filt = modFbank_YK_v2(Env, fsEnv, cf_mod) in mrGEDI
simply modified some variable names
"""

import numpy as np
from scipy.signal import butter, freqz, lfilter
import matplotlib.pyplot as plt

MFcoefIIR = None

def FilterModFB(Env, ParamMFB):
    global MFcoefIIR

    ParamMFB.setdefault('fc', [1, 2, 4, 8, 16, 32, 64, 128, 256, 512])  # v110 4 Aug 2022

    if 'fs' not in ParamMFB:
        raise ValueError('Specify ParamMFB.fs')
    if 'fc' not in ParamMFB:
        ParamMFB['fc'] = ParamMFB['fc_default']
    if 'SwPlot' not in ParamMFB:
        ParamMFB['SwPlot'] = 0

    if not MFcoefIIR:
        MFcoefIIR = MkCoefModFilter(ParamMFB)  # Making modulation filter

    LenFc = len(ParamMFB['fc'])
    NumEnv, LenEnv = Env.shape
    if NumEnv > 1:
        raise ValueError('Env should be a monoaural row vector.')

    OutMFB = np.zeros((LenFc, LenEnv))
    for nfc in range(LenFc):
        OutMFB[nfc, :] = lfilter(MFcoefIIR['b'][nfc], MFcoefIIR['a'][nfc], Env)

    ParamMFB['MFcoefIIR'] = MFcoefIIR
    return OutMFB, ParamMFB

def MkCoefModFilter(ParamMFB):
    print('--- Making modulation filter coefficients ---')
    LenFc = len(ParamMFB['fc'])

    IIR_b = np.zeros((LenFc, 4))
    IIR_a = np.zeros((LenFc, 4))

    for nfc in range(LenFc):
        if ParamMFB['fc'][nfc] == 1:  # when 1 Hz
            b, a = butter(3, ParamMFB['fc'][nfc] / (ParamMFB['fs'] / 2))
            b4 = b / a[0]
            a4 = a / a[0]

            IIR_b[nfc, :] = b4
            IIR_a[nfc, :] = a4
        else:  # Bandpass filter
            w0 = 2 * np.pi * ParamMFB['fc'][nfc] / ParamMFB['fs']
            W0 = np.tan(w0 / 2)

            Q = 1
            B0 = W0 / Q
            b = np.array([B0, 0, -B0])
            a = np.array([1 + B0 + W0**2, 2 * W0**2 - 2, 1 - B0 + W0**2])
            b3 = b / a[0]
            a3 = a / a[0]

            IIR_b[nfc, :3] = b3
            IIR_a[nfc, :3] = a3

    MFcoefIIR = {'a': IIR_a, 'b': IIR_b}

    if ParamMFB['SwPlot'] == 1:
        PlotFrspMF(ParamMFB, MFcoefIIR)
        input('Return to continue > ')

    return MFcoefIIR

def PlotFrspMF(ParamMFB, MFcoefIIR):
    plt.figure()
    Nrsl = 1024 * 4

    for nfc in range(len(ParamMFB['fc'])):
        frsp, freq = freqz(MFcoefIIR['b'][nfc], MFcoefIIR['a'][nfc], Nrsl, ParamMFB['fs'])
        plt.plot(freq, 20 * np.log10(np.abs(frsp)))

    plt.box(True)
    plt.axis([0.25, max(ParamMFB['fc']) * 2, -40, 5])
    plt.grid(True)
    plt.xscale('log')
    plt.xticks(ParamMFB['fc'])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Filter attenuation (dB)')
    plt.legend([str(fc) for fc in ParamMFB['fc']], loc='southwest')
    plt.title('Modulation filterbank')
    plt.show()
