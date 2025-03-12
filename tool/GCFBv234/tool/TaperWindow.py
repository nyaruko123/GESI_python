
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import hamming, hann, blackman

def TaperWindow(LenWin, TypeTaper, LenTaper=None, RangeSigma=3, SwPlot=0):

    if LenTaper is None:
        LenTaper = LenWin // 2
    else:
        LenTaper = int(LenTaper)  # 确保 LenTaper 是整数
    
    if (LenTaper * 2 + 1) >= LenWin:
        print('Caution (TaperWindow): No flat part.')
        if LenTaper != LenWin // 2:
            print('Caution (TaperWindow): LenTaper <-- LenWin // 2')
        LenTaper = LenWin // 2

    TypeTaper = TypeTaper[:3].upper()

    if TypeTaper == 'HAM':
        Taper = hamming(LenTaper * 2 + 1)
        TypeTaper = 'Hamming'
    elif TypeTaper == 'HAN' or TypeTaper == 'COS':
        Taper = hann(int(LenTaper * 2 + 1))
        TypeTaper = 'Hanning/Cosine'
    elif TypeTaper == 'BLA':
        Taper = blackman(LenTaper * 2 + 1)
        TypeTaper = 'Blackman'
    elif TypeTaper == 'GAU':
        nn = np.arange(-LenTaper, LenTaper + 1)
        Taper = np.exp(-(RangeSigma * nn / LenTaper) ** 2 / 2)
        TypeTaper = 'Gauss'
    else:
        Taper = np.concatenate((np.arange(1, LenTaper + 1), [LenTaper + 1], np.arange(LenTaper, 0, -1))) / (LenTaper + 1)
        TypeTaper = 'Line'

    TaperWin = np.concatenate((Taper[:LenTaper], np.ones(LenWin - LenTaper * 2), Taper[LenTaper + 1:]))

    if SwPlot == 1:
        plt.plot(TaperWin)
        plt.xlabel('Points')
        plt.ylabel('Amplitude')
        plt.title(f'TypeTaper = {TypeTaper}')
        plt.show()

    return TaperWin, TypeTaper


