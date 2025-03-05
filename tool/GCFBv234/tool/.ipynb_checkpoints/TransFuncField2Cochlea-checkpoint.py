"""
Function to calculate the transfer function from field to cochlea.
INPUT:
    ParamIn: input parameters
        TypeField2EarDrum: String for transfer function / headphones
        TypeMidEar2Cochlea: Middle ear transfer function from the ear drum to cochlea
        NfrqRsl: Number of data points, if zero, then direct out.
        fs: Sampling Frequency
        FreqCalib: Frequency at which SPL at ear drum is calibrated
        SwPlot: Switch for plot
OUTPUT:
    TransFunc: Transfer function structure
        freq: Corresponding Frequency at the data point
        Transfer functions in dB
    ParamOut: = ParamIn + DirData
"""

import numpy as np
import matplotlib.pyplot as plt
from .Freq2ERB import Freq2ERB
from .TransFuncMiddleEar_Moore16 import TransFuncMiddleEar_Moore16
from .TransFuncField2EarDrum_Set import TransFuncField2EarDrum_Set
from scipy.signal import freqz
import librosa



def TransFuncField2Cochlea(ParamIn):
    if not ParamIn:
        print("Set Parameter to default value")
        ParamIn = {}

    ParamIn.setdefault('TypeField2EarDrum', 'FreeField')
    ParamIn.setdefault('TypeMidEar2Cochlea', 'MiddleEar_Moore16')
    ParamIn.setdefault('fs', 48000)
    ParamIn.setdefault('NfrqRsl', 2048)
    ParamIn.setdefault('FreqCalib', 1000)
    ParamIn.setdefault('SwPlot', 0)
    ParamIn.setdefault('SwGetTable', 0)

    ParamOut = ParamIn
    TransFunc = {'fs': ParamIn['fs'], 'FreqCalib': ParamIn['FreqCalib']}

    freq = np.arange(ParamIn['NfrqRsl']) / ParamIn['NfrqRsl'] * ParamIn['fs'] / 2
    TransFunc['freq'] = freq

    # Data directory
    DirThisProgram = __file__
    DirData = DirThisProgram.replace('.py', '_Data/')
    ParamOut['DirData'] = DirData

    # Field to Ear Drum
    TypeField2EarDrumList = [
        'FreeField2EarDrum_Moore16',
        'DiffuseField2EarDrum_Moore16',
        'ITUField2EarDrum',
        'NoField2EarDrum',
        'HD580_L_AMLAB15', 'HD580_R_AMLAB15',
        'HD650_L_AMLAB16', 'HD650_R_AMLAB16'
    ]

    SwCrct = None
    for nc, type_name in enumerate(TypeField2EarDrumList):
        LenMatch = 3 if nc < 4 else 7
        if ParamIn['TypeField2EarDrum'][:LenMatch].upper() == type_name[:LenMatch].upper():
            SwCrct = nc
            break

    if SwCrct is None:
        print('Select "TypeField2EarDrum" from one of:')
        for type_name in TypeField2EarDrumList:
            print(f"   - {type_name}")
        raise ValueError(f'You specified "TypeField2EarDrum" as "{ParamIn["TypeField2EarDrum"]}"  <--- Check name')

    TransFunc['TypeField2EarDrum'] = TypeField2EarDrumList[SwCrct]

    # Interpolation method
    StrInterp1 = 'linear'

    if SwCrct <= 2:
        FreqTbl, FrspdBTbl = TransFuncField2EarDrum_Set(TransFunc['TypeField2EarDrum'])
        if ParamIn['fs'] / 2 > max(FreqTbl):
            FreqTbl = np.append(FreqTbl, ParamIn['fs'] / 2)
            FrspdBTbl = np.append(FrspdBTbl, FrspdBTbl[-1])
        Field2EarDrumdB = np.interp(Freq2ERB(freq), Freq2ERB(FreqTbl), FrspdBTbl, left=None, right=None)
    elif SwCrct == 3:
        Field2EarDrumdB = np.zeros_like(freq)
    else:
        NameImpRsp = f"ImpRsp_{TypeField2EarDrumList[SwCrct]}.wav"
        print(f"Read Impulse response: {NameImpRsp}")
        SndImpRsp, fsIR = librosa.load(f"{DirData}{NameImpRsp}", sr=None)
        frspIR, freqIR = freqz(SndImpRsp, 1, len(SndImpRsp), fsIR)
        Field2EarDrumdB = np.interp(Freq2ERB(freq), Freq2ERB(freqIR), 20 * np.log10(np.abs(frspIR)), left=None, right=None)

    # Compensate to 0dB at ParamIn.FreqCalib
    NumFreqCalib = np.argmin(np.abs(freq - ParamIn['FreqCalib']))
    TransFunc['freq_AtFreqCalib'] = freq[NumFreqCalib]
    TransFunc['Field2EarDrumdB'] = Field2EarDrumdB - Field2EarDrumdB[NumFreqCalib]
    TransFunc['Field2EarDrumdB_AtFreqCalib'] = TransFunc['Field2EarDrumdB'][NumFreqCalib]
    TransFunc['Field2EarDrumdB_CmpnstdB'] = Field2EarDrumdB[NumFreqCalib]

    # Ear Drum to cochlea
    if not ParamIn['TypeMidEar2Cochlea'].startswith('MiddleEar'):
        raise ValueError(f"Not Prepared yet: {ParamIn['TypeMidEar2Cochlea']}")
    else:
        FreqTbl2, FrspdBTbl2 = TransFuncMiddleEar_Moore16()
        if ParamIn['fs'] / 2 > max(FreqTbl2):
            FreqTbl2 = np.append(FreqTbl2, ParamIn['fs'] / 2)
            FrspdBTbl2 = np.append(FrspdBTbl2, FrspdBTbl2[-1])
        MidEar2CochleadB = np.interp(Freq2ERB(freq), Freq2ERB(FreqTbl2), FrspdBTbl2, left=None, right=None)

        TransFunc['MidEar2CochleadB'] = MidEar2CochleadB
        TransFunc['MidEar2CochleadB_AtFreqCalib'] = MidEar2CochleadB[NumFreqCalib]
        TransFunc['TypeMidEar2Cochlea'] = 'MiddleEar_Moore16'

    # Total: Field to cochlea
    TransFunc['Field2CochleadB'] = TransFunc['Field2EarDrumdB'] + TransFunc['MidEar2CochleadB']
    TransFunc['Field2CochleadB_AtFreqCalib'] = TransFunc['Field2CochleadB'][NumFreqCalib]
    TransFunc['TypeField2CochleadB'] = f"{TransFunc['TypeField2EarDrum']} + {TransFunc['TypeMidEar2Cochlea']}"

    print(f"TypeField2CochleadB: {TransFunc['TypeField2CochleadB']}")
    print(f"TransFunc.freq_AtFreqCalib = {TransFunc['freq_AtFreqCalib']} Hz  ( <-- {ParamIn['FreqCalib']} Hz )")
    print(f"TransFunc.Field2EarDrumdB_AtFreqCalib = {TransFunc['Field2EarDrumdB_AtFreqCalib']:.3f} dB")
    print(f"                        (Compensated for {TransFunc['Field2EarDrumdB_CmpnstdB']:.3f} dB)")
    print(f"TransFunc.MidEar2CochleadB_AtFreqCalib = {TransFunc['MidEar2CochleadB_AtFreqCalib']:.3f} dB")
    print(f"TransFunc.Field2CochleadB_AtFreqCalib = {TransFunc['Field2CochleadB_AtFreqCalib']:.3f} dB")

    # Plot data
    if ParamIn['SwPlot'] == 1:
        plt.figure()
        plt.subplot(3, 1, 1)
        plt.semilogx(freq, TransFunc['Field2EarDrumdB'], label='Field2EarDrumdB')
        plt.plot(TransFunc['freq_AtFreqCalib'], TransFunc['Field2EarDrumdB_AtFreqCalib'], 'rx')
        plt.grid(True)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Gain (dB)')
        plt.title(f"Frequency response: {TransFunc['TypeField2EarDrum']}, Gain normalized at {TransFunc['freq_AtFreqCalib']} Hz")
        plt.axis([10, 30000, -20, 20])

        plt.subplot(3, 1, 2)
        plt.semilogx(freq, TransFunc['MidEar2CochleadB'], label='MidEar2CochleadB')
        plt.grid(True)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Gain (dB)')
        plt.title(f"Frequency response: {TransFunc['TypeMidEar2Cochlea']}")
        plt.axis([10, 30000, -30, 10])

        plt.subplot(3, 1, 3)
        plt.semilogx(freq, TransFunc['Field2CochleadB'], label='Field2CochleadB')
        plt.grid(True)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Gain (dB)')
        plt.title("Frequency response: Total Transfer Function from Field to cochlea")
        plt.axis([10, 30000, -30, 10])

        plt.show()

    return TransFunc, ParamOut


