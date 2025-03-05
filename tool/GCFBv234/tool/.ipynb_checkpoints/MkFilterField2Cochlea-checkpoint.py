from scipy.fft import ifft
import numpy as np
from .TransFuncField2Cochlea import TransFuncField2Cochlea
from scipy.signal import freqz
from .OutMidCrct import OutMidCrct
from scipy.signal import remez
from .TaperWindow import TaperWindow
def rceps(signal):
    spectrum = np.fft.fft(signal)
    log_spectrum = np.log(np.abs(spectrum))
    cepstrum = np.real(ifft(log_spectrum))
    return spectrum, cepstrum

def MkFilterField2Cochlea(StrCrct, fs, SwFwdBwd=1, SwPlot=0):
    global Param_Keep, fs_Keep, FIRCoefFwd_Keep, FIRCoefBwd_Keep
    global TypeField2EarDrum_Keep, TypeMidEar2Cochlea_Keep

    # Initial setup
    if fs > 48000:
        print(f'MkFilterField2Cochlea: Sampling rate of {fs} (Hz) (> 48000 (Hz)) is not recommended.')
        print('<-- Transfer function is only defined below 16000 (Hz).')
    Param = {'fs': fs}

    if StrCrct.lower() in ['freefield', 'ff']:
        SwType = 1
        Param['TypeField2EarDrum'] = 'FreeField'
        Param['TypeMidEar2Cochlea'] = 'MiddleEar'

    elif StrCrct.lower() in ['diffusefield', 'df']:
        SwType = 2
        Param['TypeField2EarDrum'] = 'DiffuseField'
        Param['TypeMidEar2Cochlea'] = 'MiddleEar'

    elif StrCrct.lower() == 'itu':
        SwType = 3
        Param['TypeField2EarDrum'] = 'ITU'
        Param['TypeMidEar2Cochlea'] = 'MiddleEar'

    elif StrCrct.lower() in ['eardrum', 'ed']:
        SwType = 4
        Param['TypeField2EarDrum'] = 'NoField2EarDrum'
        Param['TypeMidEar2Cochlea'] = 'MiddleEar'

    elif StrCrct.lower() == 'elc':
        SwType = 10
        Param['TypeField2CochleadB'] = 'ELC'
        Param['TypeField2EarDrum'] = 'NoUse_ELC'
        Param['TypeMidEar2Cochlea'] = 'NoUse_ELC'

    else:
        raise ValueError('Specify: FreeField (FF) / DiffuseField (DF) / ITU / EarDrum (ED) / ELC')

    if SwFwdBwd == 1:
        Param['NameFilter'] = '(Forward) FIR minimum phase filter'
    elif SwFwdBwd == -1:
        Param['NameFilter'] = '(Backward) FIR minimum phase inverse filter'
    else:
        raise ValueError('Specify SwFwdBwd: (1) Forward, (-1) Backward.')

    Param['NameFilter'] = f'[{StrCrct}] {Param["NameFilter"]}'

    # No Calculation. Restoring from the kept data
    if (TypeField2EarDrum_Keep == Param['TypeField2EarDrum'] and
        TypeMidEar2Cochlea_Keep == Param['TypeMidEar2Cochlea'] and
        fs_Keep == fs and SwPlot == 0):
        if SwFwdBwd == 1 and len(FIRCoefFwd_Keep) > 20:
            FIRCoef = FIRCoefFwd_Keep
            print(f'*** MkFilterField2Cochlea: Restoring {Param["NameFilter"]} ***')
            Param = Param_Keep
            return FIRCoef, Param

        elif SwFwdBwd == -1 and len(FIRCoefBwd_Keep) > 20:
            FIRCoef = FIRCoefBwd_Keep
            print(f'*** MkFilterField2Cochlea: Restoring {Param["NameFilter"]} ***')
            Param = Param_Keep
            return FIRCoef, Param

    # Generating filter at the first time
    print(f'*** MkFilterField2Cochlea: Generating {Param["NameFilter"]} ***')

    if SwType <= 4:
        Param['TypeMidEar2Cochlea'] = 'MiddleEar'
        TransFunc = TransFuncField2Cochlea(Param)  # Undefined function
        FrspCrct = 10 ** (TransFunc['Field2CochleadB'] / 20)
        freq = TransFunc['freq']
        Param['TypeField2CochleadB'] = TransFunc['TypeField2CochleadB']
    elif SwType == 10:
        Nrslt = 2048
        crctPwr, freq = OutMidCrct(StrCrct, Nrslt, fs, 0)  # Undefined function
        FrspCrct = np.sqrt(crctPwr)

    if SwFwdBwd == -1:
        FrspCrct = 1 / np.maximum(FrspCrct, 0.1)

    try:
        LenCoef = 200
        NCoef = int(LenCoef / 16000 * fs / 2) * 2
        FIRCoef = remez(NCoef, freq / fs * 2, FrspCrct)
    except:
        print('-- For octave compatibility --')
        LenCoef = 50
        NCoef = int(LenCoef / 16000 * fs / 2) * 2
        FIRCoef = remez(NCoef, freq / fs * 2, FrspCrct)

    Win = TaperWindow(len(FIRCoef), 'han', LenCoef / 10)  # Undefined function
    FIRCoef = Win * FIRCoef

    _, x_mp = rceps(FIRCoef)
    FIRCoef = x_mp[:len(x_mp) // 2]

    # Keep records for fast processing
    if SwFwdBwd == 1:
        FIRCoefFwd_Keep = FIRCoef
    elif SwFwdBwd == -1:
        FIRCoefBwd_Keep = FIRCoef

    fs_Keep = fs
    Param_Keep = Param
    TypeField2EarDrum_Keep = Param['TypeField2EarDrum']
    TypeMidEar2Cochlea_Keep = Param['TypeMidEar2Cochlea']

    # Plot
    if SwPlot == 1:
        Nrsl = len(FrspCrct)
        frsp, freq2 = freqz(FIRCoef, 1, Nrsl, fs)
        import matplotlib.pyplot as plt
        plt.subplot(2, 1, 1)
        plt.plot(FIRCoef)
        plt.xlabel('Sample')
        plt.ylabel('Amplitude')
        plt.title(f'Type: {Param["TypeField2EarDrum"]}')

        plt.subplot(2, 1, 2)
        plt.plot(freq2, np.abs(frsp), freq, FrspCrct, '--')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude (linear term)')
        ELCError = np.mean((np.abs(frsp) - FrspCrct) ** 2) / np.mean(FrspCrct ** 2)
        ELCErrordB = 10 * np.log10(ELCError)

        print(f'Fitting Error: {ELCErrordB} (dB)')
        if ELCErrordB > -30:
            print(f'Warning: Error in ELC correction = {ELCErrordB} dB > -30 dB')

    return FIRCoef, Param
