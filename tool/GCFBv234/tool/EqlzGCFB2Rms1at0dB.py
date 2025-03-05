import numpy as np

def EqlzGCFB2Rms1at0dB(GCval, StrFloor=None):
    """
    Convert GC output to be relative to absolute threshold 0 dB.

    INPUT:
        GCval: The output of GCFBv231, where rms(snd) == 1 corresponds to 30 dB
        StrFloor:
            'NoiseFloor' - Add Gaussian noise (rms(randn) == 1)
            'ZeroFloor' - Set values less than 1 to 0

    OUTPUT:
        GCreAT: GC relative to AbsThreshold 0 dB (rms(snd) == 1 corresponds to 0 dB)

    Note:
        Snd --> Eqlz2MeddisHCLevel --> GCFB
        GC output level is the same as the MeddisHCLevel.
        This function converts the level from MeddisHCLevel to rms(s(t)) == 1 --> 0 dB.
        Use this when the absolute threshold is set to 0 dB as in GCFBv231.
        GCFB --> EqlzGCFB2Rms1at0dB --> GCFBeqlz
    """

    MeddisHCLeveldB_RMS1 = 30  # used in GCFB level set
    GCreAT = 10 ** (MeddisHCLeveldB_RMS1 / 20) * GCval

    if StrFloor is not None:
        if StrFloor == 'NoiseFloor':
            GCreAT = GCreAT + np.random.randn(*GCreAT.shape)  # adding Gaussian noise
        elif StrFloor == 'ZeroFloor':
            GCreAT = np.maximum(GCreAT - 1, 0)  # cutoff values less than 1
        else:
            raise ValueError("Specify StrFloor properly: 'NoiseFloor' or 'ZeroFloor'")

    return GCreAT
