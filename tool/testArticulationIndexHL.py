"""
Some tests on Articulation index, HL0dB
Irino, T.
Created: 8 Feb 22
Modified: 8 Feb 22

AI from
Daniel R. Raichel, "THE SCIENCE AND APPLICATIONS OF ACOUSTICS (2nd Ed.)",
Springer. 2006
"""

import numpy as np
import matplotlib.pyplot as plt

def Freq2ERB(freq):
    # Placeholder function for Freq2ERB
    pass

def ERB2Freq(erb):
    # Placeholder function for ERB2Freq
    pass

def F0limit2SSIweight(SSIparam):
    # Placeholder function for F0limit2SSIweight
    pass

NumOneThird = 7 * 3 + 1
FreqOneThird = np.exp(np.linspace(np.log(125), np.log(8000), NumOneThird))
NumOneSixth = 7 * 6 + 1
FreqOneSixth = np.exp(np.linspace(np.log(125), np.log(8000), NumOneSixth))
n1000 = np.argmin(np.abs(FreqOneSixth - 1000))
Bw1000 = FreqOneSixth[n1000 + 1] - FreqOneSixth[n1000 - 1]

AIfreq = np.array([200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000])
AIweight = np.array([4, 10, 10, 14, 14, 20, 20, 24, 30, 37, 37, 34, 34, 24, 20])

ERBband = np.linspace(Freq2ERB(100), Freq2ERB(8000), 100)
Fr1 = ERB2Freq(ERBband)
_, ERBw = Freq2ERB(1000)
AIweightERB = np.interp(Fr1, AIfreq, AIweight, left=None, right=None) * ERBw / Bw1000
pwr = 0.5
AIweightERB = np.maximum(AIweightERB, 0) ** pwr * 5

SSIparam = {'Fr1': Fr1, 'F0_limit': 110}
SSIweight = F0limit2SSIweight(SSIparam)

plt.subplot(2, 1, 1)
plt.semilogx(AIfreq, AIweight, 'o-', label='AI weight')
plt.semilogx(Fr1, AIweightERB, label='AI weight ERB')
plt.semilogx(Fr1, SSIweight * 10, label='SSI weight')
plt.grid(True)
plt.xlabel('Frequency (Hz)')
plt.ylabel('AI weight (ERB interp. pwr0.3)')
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(Freq2ERB(AIfreq), AIweight, 'o-', label='AI weight')
plt.plot(ERBband, AIweightERB, label='AI weight ERB')
plt.plot(ERBband, SSIweight * 10, label='SSI weight')
plt.grid(True)
plt.xlabel('ERB_N number (Cam)')
plt.ylabel('AI weight (ERB interp. pwr0.3)')
plt.legend()

plt.show()

indices_below_200 = np.where(Fr1 < 200)[0]
print(indices_below_200)
