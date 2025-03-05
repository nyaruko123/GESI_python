import numpy as np

def Eqlz2MeddisHCLevel(SndIn, OutLeveldB=None, InputRms1SPLdB=None):
    """
    Equalize sound to Meddis Hair Cell Level.

    INPUT:
        SndIn: Input sound
        OutLeveldB: Output level (RMS level), conventional method
        InputRms1SPLdB: SPL(dB) of input sound digital level rms(s(t))=1

    OUTPUT:
        SndEqMds: Equalized Sound (rms value of 1 is 30 dB SPL)
        AmpdB: 3 values in dB [OutputLevel_dB, CompensationValue_dB, SourceLevel_dB]

    Reference: Meddis (1986), JASA, 79(3), pp.702-711.

    rms(s(t)) == sqrt(mean(s^2)) == 1   --> 30 dB SPL
    rms(s(t)) == sqrt(mean(s^2)) == 10  --> 50 dB SPL
    rms(s(t)) == sqrt(mean(s^2)) == 100 --> 70 dB SPL

    Usage:
        1. Conventional: SndEqMds, AmpdB = Eqlz2MeddisHCLevel(SndIn, 65)
        2. InputRms1SPLdB: SndEqMds, AmpdB = Eqlz2MeddisHCLevel(SndIn, None, 65 + 26)
    """

    if OutLeveldB is not None and InputRms1SPLdB is None:
        # Conventional method
        SourceLevel = np.sqrt(np.mean(SndIn ** 2)) * 10 ** (30 / 20)
        AmpCmpnst = (10 ** (OutLeveldB / 20)) / SourceLevel
        SndEqMds = AmpCmpnst * SndIn

        SourceLeveldB = 20 * np.log10(SourceLevel)
        CmpnstdB = 20 * np.log10(AmpCmpnst)

    elif InputRms1SPLdB is not None:
        if OutLeveldB is not None:
            raise ValueError("You need to set OutLevel = None when using InputRms1SPLdB. See Usage.")

        SourceLeveldB = 20 * np.log10(np.sqrt(np.mean(SndIn ** 2))) + InputRms1SPLdB
        OutLeveldB = SourceLeveldB

        CmpnstdB = InputRms1SPLdB - 30
        AmpCmpnst = 10 ** (CmpnstdB / 20)
        SndEqMds = AmpCmpnst * SndIn

    else:
        raise ValueError("Invalid input: either OutLeveldB or InputRms1SPLdB must be provided.")

    AmpdB = [OutLeveldB, CmpnstdB, SourceLeveldB]

    return SndEqMds, AmpdB
