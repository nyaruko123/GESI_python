import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz
import os
from scipy.io import wavfile
from scipy.signal import resample

def OutMidCrct(StrCrct, NfrqRsl=None, fs=None, SwPlot=None, SwInterp=None):
    """
    Correction function for ELC/MAF/MAP/MidEar/ER4B.

    INPUT:
        StrCrct: String for Correction ELC/MAF/MAP/MidEar/ER4B
        NfrqRsl: Number of data points, if zero, then direct out.
        fs:      Sampling Frequency
        SwPlot:  Switch for plot
        SwInterp: Switch for interpolation method (24 Jul 2021)

    OUTPUT:
        CrctLinPwr: Correction value in LINEAR POWER
                    This is defined as: CrctLiPwr = 10^(-FreqChardB_toBeCmpnstd/10);
        freq: Corresponding Frequency at the data point
        FreqChardB_toBeCmpnstd: Frequency char of ELC/MAP... dB to be compensated for filterbank
                                (defined By Glasberg and Moore.)
    """

    if StrCrct is None:
        raise ValueError("Input string StrCrct is required.")

    if NfrqRsl is None:
        NfrqRsl = 0
    if fs is None:
        fs = 32000  # Default sampling frequency
    if SwPlot is None:
        SwPlot = 1
    if SwInterp is None:
        SwInterp = 1  # Original Spline interpolation

    if StrCrct.upper() == 'ER4B':
        return OutMidCrct_ER4B(NfrqRsl, fs, SwPlot)

    # Conventional ELC/MAF/MidEar
    f1 = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60, 70,
                   80, 90, 100, 125, 150, 177, 200, 250, 300, 350,
                   400, 450, 500, 550, 600, 700, 800, 900, 1000, 1500,
                   2000, 2500, 2828, 3000, 3500, 4000, 4500, 5000,
                   5500, 6000, 7000, 8000, 9000, 10000, 12748, 15000])

    ELC = np.array([31.8, 26.0, 21.7, 18.8, 17.2, 15.4, 14.0, 12.6, 11.6, 10.6,
                    9.2, 8.2, 7.7, 6.7, 5.3, 4.6, 3.9, 2.9, 2.7, 2.3,
                    2.2, 2.3, 2.5, 2.7, 2.9, 3.4, 3.9, 3.9, 3.9, 2.7,
                    0.9, -1.3, -2.5, -3.2, -4.4, -4.1, -2.5, -0.5, 2.0, 5.0,
                    10.2, 15.0, 17.0, 15.5, 11.0, 22.0])

    MAF = np.array([73.4, 65.2, 57.9, 52.7, 48.0, 45.0, 41.9, 39.3, 36.8, 33.0,
                    29.7, 27.1, 25.0, 22.0, 18.2, 16.0, 14.0, 11.4, 9.2, 8.0,
                    6.9, 6.2, 5.7, 5.1, 5.0, 5.0, 4.4, 4.3, 3.9, 2.7,
                    0.9, -1.3, -2.5, -3.2, -4.4, -4.1, -2.5, -0.5, 2.0, 5.0,
                    10.2, 15.0, 17.0, 15.5, 11.0, 22.0])

    f2 = np.array([125, 250, 500, 1000, 1500, 2000, 3000,
                   4000, 6000, 8000, 10000, 12000, 14000, 16000])

    MAP = np.array([30.0, 19.0, 12.0, 9.0, 11.0, 16.0, 16.0,
                    14.0, 14.0, 9.9, 24.7, 32.7, 44.1, 63.7])

    f3 = np.array([1, 20, 25, 31.5, 40, 50, 63, 80, 100, 125,
                   160, 200, 250, 315, 400, 500, 630, 750, 800, 1000,
                   1250, 1500, 1600, 2000, 2500, 3000, 3150, 4000, 5000, 6000,
                   6300, 8000, 9000, 10000, 11200, 12500, 14000, 15000, 16000, 20000])

    MID = np.array([50, 39.15, 31.4, 25.4, 20.9, 18, 16.1, 14.2, 12.5, 11.13,
                    9.71, 8.42, 7.2, 6.1, 4.7, 3.7, 3.0, 2.7, 2.6, 2.6,
                    2.7, 3.7, 4.6, 8.5, 10.8, 7.3, 6.7, 5.7, 5.7, 7.6,
                    8.4, 11.3, 10.6, 9.9, 11.9, 13.9, 16.0, 17.3, 17.8, 20.0])

    frqTbl = []
    TblFreqChar = []
    ValHalfFs = None

    if StrCrct.upper() == 'ELC':
        frqTbl = f1
        TblFreqChar = ELC
        ValHalfFs = 130
    elif StrCrct.upper() == 'MAF':
        frqTbl = f1
        TblFreqChar = MAF
        ValHalfFs = 130
    elif StrCrct.upper() == 'MAP':
        frqTbl = f2
        TblFreqChar = MAP
        ValHalfFs = 180
    elif StrCrct.upper() == 'MIDEAR':
        frqTbl = f3
        TblFreqChar = MID
        ValHalfFs = 23
    elif StrCrct.upper() == 'NO':
        pass
    else:
        raise ValueError("Specify correction: ELC / MAF / MAP / MidEar or NO correction.")

    # Additional dummy data for high sampling frequency
    if fs > 32000:
        frqTbl = np.append(frqTbl, fs / 2)
        TblFreqChar = np.append(TblFreqChar, ValHalfFs)
        frqTbl, indx = np.unique(frqTbl, return_index=True)
        TblFreqChar = TblFreqChar[indx]

    if NfrqRsl <= 0:
        freq = frqTbl
        FreqChardB_toBeCmpnstd = TblFreqChar
    else:
        freq = np.arange(NfrqRsl) / NfrqRsl * fs / 2
        if StrCrct.upper().startswith('NO'):
            FreqChardB_toBeCmpnstd = np.zeros_like(freq)
        else:
            if SwInterp == 1:  # Original version
                FreqChardB_toBeCmpnstd = np.interp(freq, frqTbl, TblFreqChar)  # Spline interpolation
            elif SwInterp == 2:  # Added 24 Jul 2021
                FreqChardB_toBeCmpnstd = np.interp(np.log10(freq), np.log10(frqTbl), TblFreqChar, left=None,
                                                   right=None)  # Linear interpolation on log scale
            else:
                raise ValueError("Specify SwInterp")

    if SwPlot == 1:
        str_msg = f"*** Frequency Characteristics ({StrCrct.upper()}): Its inverse will be corrected. ***"
        print(str_msg)
        plt.plot(frqTbl, TblFreqChar, label='Original')
        plt.plot(freq, FreqChardB_toBeCmpnstd, '--', label='Interpolated')
        plt.title(str_msg)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Level (dB)')
        plt.legend()
        plt.show()

    CrctLinPwr = 10 ** (-FreqChardB_toBeCmpnstd / 10)  # in Linear Power
    return CrctLinPwr, freq, FreqChardB_toBeCmpnstd


def OutMidCrct_ER4B(NfrqRsl, fs, SwPlot):
    """
    Correction function for ER-4B.

    INPUT:
        NfrqRsl: Number of data points.
        fs:      Sampling Frequency.
        SwPlot:  Switch for plot.

    OUTPUT:
        CrctLinPwr: Correction value in LINEAR POWER.
        freq: Corresponding Frequency at the data point.
        FreqChardB_toBeCmpnstd: Frequency char of ER-4B dB to be compensated.
    """

    DirEarPhone = f"{os.getenv('USERPROFILE')}/Data/Measurement/MsrRsp_ER-4B_9Feb09/"
    NameImpRsp = 'TSP48-24_ImpRsp1.wav'

    # 读取 WAV 文件
    fsEP, ImpRspEP = wavfile.read(os.path.join(DirEarPhone, NameImpRsp))

    # 计算比特深度
    Nbit = ImpRspEP.dtype.itemsize * 8

    ImpRspEP = ImpRspEP[1024:]  # Skip the first 1024 samples

    if fsEP != fs:
        ImpRspEP = resample(ImpRspEP, fs, fsEP)

    Frsp, freq = freqz(ImpRspEP, 1, NfrqRsl, fs)
    CrctLinPwr = abs(Frsp ** 2)
    FreqChardB_toBeCmpnstd = -10 * np.log10(CrctLinPwr)  # It is inverse

    if SwPlot == 1:
        plt.plot(freq, FreqChardB_toBeCmpnstd)
        plt.title("ER-4B Correction")
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Level (dB)')
        plt.show()

    return CrctLinPwr, freq, FreqChardB_toBeCmpnstd

# 注意：在此代码中，`wavread`, `resample`, 和 `freqz` 是未定义的函数，需要在使用前定义或导入。
