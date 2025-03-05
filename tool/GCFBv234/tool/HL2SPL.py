import numpy as np
from .SPLatHL0dB_Table import SPLatHL0dB_Table

def HL2SPL(freq, HLdB):
    Table1 = SPLatHL0dB_Table()  # Retrieve the table data
    FreqRef = Table1['freq']
    SPLdBatHL0dB = Table1['SPLatHL0dB']

    # **确保 freq 和 HLdB 是 numpy 数组**
    freq = np.atleast_1d(freq)  # 如果 freq 是整数，转换成数组
    HLdB = np.atleast_1d(HLdB)  # 同理转换 HLdB

    if len(freq) != len(HLdB):
        raise ValueError('Length of freq & HLdB should be the same.')

    SPLdB = []
    for nf in range(len(freq)):
        nfreq = np.where(FreqRef == freq[nf])[0]
        if len(nfreq) == 0:
            raise ValueError('Frequency should be one of 125*2^n & 750*n (Hz) <= 8000.')

        SPLdB.append(HLdB[nf] + SPLdBatHL0dB[nfreq[0]])

    return SPLdB
