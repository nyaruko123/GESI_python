"""
Function [ERBrate, ERBwidth] = Freq2ERB(cf)
    INPUT:
        cf: Center frequency
    OUTPUT:
        ERBrate: ERB_N rate
        ERBwidth: ERB_N Bandwidth

    Ref: Glasberg and Moore: Hearing Research, 47 (1990), 103-138
         For different formulae (years), see Freq2ERBYear.m
"""
import numpy as np

def Freq2ERB(cf):
    if cf is None:
        raise ValueError("Center frequency (cf) must be provided.")

    ERBrate = 21.4 * np.log10(4.37 * cf / 1000 + 1)
    ERBwidth = 24.7 * (4.37 * cf / 1000 + 1)

    # Warning for Frequency Range
    cfmin = 50
    cfmax = 12000
    if np.min(cf) < cfmin or np.max(cf) > cfmax:
        print("Warning: Min or max frequency exceeds the proper ERB range:")
        print(f"         {cfmin}(Hz) <= Fc <= {cfmax}(Hz).")

    return ERBrate, ERBwidth


