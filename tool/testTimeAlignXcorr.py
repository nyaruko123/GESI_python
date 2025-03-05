"""
Test Time Alignment for GESI
Irino, T.
Created: 6 Jun 2022 IT from mrGEDI GEDI_TimeAlign
Modified: 6 Jun 2022 IT
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import correlate
from numpy.random import default_rng

def TimeAlignXcorr(SndTestIn, SndRefIn):
    # Placeholder function for TimeAlignXcorr
    pass

# Read audio files
fs, SndTest1 = wavfile.read('../wav_sample/sample_sp1.wav')
_, SndRef1 = wavfile.read('../wav_sample/sample_sp_clean.wav')
rng = default_rng(1234)

SndTestIn = SndTest1.flatten()
# Add additional noise (optional)
if False:
    RmsSndTest = np.sqrt(np.mean(SndTestIn**2))
    SndTestIn += 10 * RmsSndTest * rng.standard_normal(SndTestIn.shape)
    SndTestIn = RmsSndTest * SndTestIn / np.sqrt(np.mean(SndTestIn**2))

SndRefIn = SndRef1.flatten()
nz = 0.1 * rng.standard_normal(10)
SndTestIn = np.concatenate([nz, -SndRefIn, nz])
SndTestIn = np.concatenate([nz, SndRefIn, nz])

# Play sound (requires sounddevice or similar library)
# import sounddevice as sd
# sd.play(SndTestIn, fs)
# sd.wait()

SndTestOut, ParamTA = TimeAlignXcorr(SndTestIn, SndRefIn)

Error = np.sqrt(np.mean((SndRefIn - SndTestOut)**2))
print("RMS Error:", Error)

plt.figure(1)
plt.clf()
plt.plot(ParamTA['Lag'], ParamTA['XcorrSnd'])
plt.title('Cross-correlation')

plt.figure(2)
plt.clf()
bias = 0.3
plt.subplot(2, 1, 1)
plt.plot(range(len(SndRefIn)), SndRefIn, label='Reference')
plt.plot(range(len(SndTestIn)), SndTestIn + bias, label='Test + Bias')
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(range(len(SndRefIn)), SndRefIn, label='Reference')
plt.plot(range(len(SndTestOut)), SndTestOut + bias, label='Aligned + Bias')
plt.legend()

plt.show()
