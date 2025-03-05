import numpy as np


def GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, CompressionHealth, PindB):
    """
    Compute asymmetric function output in dB scale.

    Parameters:
        GCparam: dict
            Parameters for the gammatone filter bank.
        GCresp: dict
            Response data of the filter bank.
        Fr1query: float
            Frequency query specified by Fr1, usually used in specifying filter bank frequency (not Fp1).
        CompressionHealth: float
            Compression health parameter.
        PindB: float
            Input level in dB.

    Returns:
        AFoutdB: float
            Output of the asymmetric function in dB.
        IOfuncdB: float
            Input-output function in dB.
        GCparam: dict
            Updated GCparam with the normalization parameter.
    """
    # Default normalization value (changing it to 200 has little impact on GCFBv231 output, only a shift in dB scale)
    GCparam['AsymFunc_NormdB'] = 100

    # Compute asymmetric function output
    AFoutLin = CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, PindB)
    AFoutLinNorm = CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, GCparam['AsymFunc_NormdB'])

    # Convert to dB scale
    AFoutdB = 20 * np.log10(AFoutLin / AFoutLinNorm)
    IOfuncdB = AFoutdB + PindB

    return AFoutdB, IOfuncdB, GCparam


def CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, PindB):
    """
    Compute the asymmetric function in linear scale.

    Parameters:
        GCparam: dict
            Parameters for the gammatone filter bank.
        GCresp: dict
            Response data of the filter bank.
        Fr1query: float
            Frequency query specified by Fr1.
        CompressionHealth: float
            Compression health parameter.
        PindB: float
            Input level in dB.

    Returns:
        AFoutLin: float
            Output of the asymmetric function in linear scale.
    """
    # Find the closest frequency index
    nch = np.argmin(np.abs(GCparam['Fr1'] - Fr1query))

    # Get parameters
    Fp1 = GCresp['Fp1'][nch]
    frat = GCresp['frat0Pc'][nch] + GCresp['frat1val'][nch] * (PindB - GCresp['PcHPAF'][nch])
    Fr2 = frat * Fp1

    # Undefined function: Freq2ERB (Needs implementation or user-provided function)
    # _, ERBw2 = Freq2ERB(Fr2)
    ERBw2 = None  # Placeholder, should be replaced with actual function call

    # Compute b2E and c2CH
    b2E = GCresp['b2val'][nch] * ERBw2  # ERBw2 needs to be computed correctly
    c2CH = CompressionHealth * GCresp['c2val'][nch]

    # Compute the asymmetric function output in linear scale
    AFoutLin = np.exp(c2CH * np.arctan2(Fp1 - Fr2, b2E))

    return AFoutLin

