"""
Test Filtering by Modulation Filterbank (MFB)
Irino, T.
Created: 15 Oct 2022
Modified: 15 Oct 2022
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz

def FilterModFB(ImpLP, ParamMFB):
    # Placeholder function for FilterModFB
    pass

ParamMFB = {
    'fs': 2000,  # sampling rate
    'fcutEnv': 150  # cutoff frequency for the envelope LPF
}

Impulse = np.concatenate(([0, 1], np.zeros(ParamMFB['fs'] - 1)))

# Envelope LPF
b, a = butter(1, ParamMFB['fcutEnv'] / (ParamMFB['fs'] / 2))
ImpLP = lfilter(b, a, Impulse)

OutMFB, ParamMFB = FilterModFB(ImpLP, ParamMFB)

NchM = OutMFB.shape[0]

AmpPeakdB = []

plt.figure()
for nch in range(NchM):
    frsp, freq = freqz(OutMFB[nch, :], worN=2**np.ceil(np.log2(ParamMFB['fs'])), fs=ParamMFB['fs'])
    plt.plot(freq, 20 * np.log10(np.abs(frsp)))
    AmpPeakdB.append(np.max(20 * np.log10(np.abs(frsp))))
    plt.hold(True)

plt.xlabel('Modulation Freq. (Hz)')
plt.ylabel('Gain (dB)')
plt.axis([0, ParamMFB['fs'] / 2, -30, 5])
plt.grid(True)

plt.plot(ParamMFB.get('fc', []), AmpPeakdB, '--')
plt.show()

AmpPeakdB = np.array(AmpPeakdB)
print(AmpPeakdB)
