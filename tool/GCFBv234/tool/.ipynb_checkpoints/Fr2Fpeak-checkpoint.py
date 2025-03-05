"""
Function [fpeak, ERBw] = Fr2Fpeak(n, b, c, fr)
    INPUT:
        n, b, c: Gammachirp parameters
        fr: Frequency
    OUTPUT:
        fpeak: Peak frequency
        ERBw: ERB width at fr
"""

import numpy as np
from .Freq2ERB import Freq2ERB
def Fr2Fpeak(n, b, c, fr):
    if n is None or b is None or c is None or fr is None:
        raise ValueError("All parameters (n, b, c, fr) must be provided.")

    n = np.array(n).flatten()
    b = np.array(b).flatten()
    c = np.array(c).flatten()
    fr = np.array(fr).flatten()

    _, ERBw = Freq2ERB(fr)
    fpeak = fr + c * ERBw * b / n

    return fpeak, ERBw


