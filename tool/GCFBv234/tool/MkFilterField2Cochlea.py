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
    # 初始化
    TypeField2EarDrum_Keep = None
    TypeMidEar2Cochlea_Keep = None
    fs_Keep = None
    FIRCoefFwd_Keep = []
    FIRCoefBwd_Keep = []
    
    if fs > 48000:
        print(f'MkFilterField2Cochlea: Sampling rate of {fs} Hz (> 48000 Hz) is not recommended.')
        print('<-- Transfer function is only defined below 16000 Hz.')
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

    print(f'*** MkFilterField2Cochlea: Generating {Param["NameFilter"]} ***')

    if SwType <= 4:
        TransFunc, ParamOut = TransFuncField2Cochlea(Param)
        FrspCrct = 10 ** (TransFunc['Field2CochleadB'] / 20)
        freq = TransFunc['freq']
        Param['TypeField2CochleadB'] = TransFunc['TypeField2CochleadB']
    elif SwType == 10:
        Nrslt = 2048
        crctPwr, freq = OutMidCrct(StrCrct, Nrslt, fs, 0)
        FrspCrct = np.sqrt(crctPwr)

    if SwFwdBwd == -1:
        FrspCrct = 1 / np.maximum(FrspCrct, 0.1)

    try:
        LenCoef = 200
        NCoef = max(int(LenCoef / 16000 * fs / 2) * 2, 2)  # 保证 NCoef 是偶数，且至少为 2

        freq_norm = np.array(freq) / (fs / 2)  # 归一化频率
        freq_norm = np.clip(freq_norm, 0, 1)  # 确保频率在 0~1 之间

        #===== 修改关键部分：修正频带与期望响应的维度关系 =====#
        bands = []
        desired = []
        # 构造连续的频带区间
        for i in range(len(freq_norm)-1):
            bands.extend([freq_norm[i], freq_norm[i+1]])  # 每个频带两个边界
            desired.append(FrspCrct[i])  # 每个频带对应一个期望值
        
        weights = np.ones(len(desired))  # 权重数组长度与desired一致
        
        # 添加维度验证
        if len(bands) != 2 * len(desired):
            raise ValueError(f"频带维度异常 bands({len(bands)}) != 2*desired({len(desired)})")

        FIRCoef = remez(NCoef, bands, desired, weight=weights, fs=fs)

    except Exception as e:
        print('-- For octave compatibility --')
        print(f'Error: {e}')
        # 异常处理时重新生成有效参数
        LenCoef = 50
        NCoef = max(int(LenCoef / 16000 * fs / 2) * 2, 2)
        # 构造降采样后的频带参数
        downsample_ratio = max(1, len(freq_norm)//100)  # 防止频点过多
        freq_norm_simple = freq_norm[::downsample_ratio]
        FrspCrct_simple = FrspCrct[::downsample_ratio]
        
        bands = []
        desired = []
        for i in range(len(freq_norm_simple)-1):
            bands.extend([freq_norm_simple[i], freq_norm_simple[i+1]])
            desired.append(FrspCrct_simple[i])
        
        FIRCoef = remez(NCoef, bands, desired, weight=weights[:len(desired)], fs=fs)

    Win = TaperWindow(len(FIRCoef), 'han', LenCoef / 10)
    FIRCoef = Win * FIRCoef

    _, x_mp = rceps(FIRCoef)
    FIRCoef = x_mp[:len(x_mp) // 2]

    if SwFwdBwd == 1:
        FIRCoefFwd_Keep = FIRCoef
    elif SwFwdBwd == -1:
        FIRCoefBwd_Keep = FIRCoef

    fs_Keep = fs
    Param_Keep = Param
    TypeField2EarDrum_Keep = Param['TypeField2EarDrum']
    TypeMidEar2Cochlea_Keep = Param['TypeMidEar2Cochlea']

    if SwPlot == 1:
        Nrsl = len(FrspCrct)
        frsp, freq2 = freqz(FIRCoef, 1, Nrsl, fs)
        import matplotlib.pyplot as plt
        plt.plot(freq2, np.abs(frsp), label="Filter Response")
        plt.plot(freq, FrspCrct, '--', label="Target Response")
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        plt.show()

    return FIRCoef, Param
