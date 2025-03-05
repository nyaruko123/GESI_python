"""
Function [cf, ERBwidth] = ERB2Freq(ERBrate)
    INPUT:
        ERBrate: ERB rate
    OUTPUT:
        cf: Center frequency (Hz)
        ERBwidth: ERB bandwidth (Hz)

    Ref: Glasberg and Moore: Hearing Research, 47 (1990), 103-138
         For different formulae (years), see ERB2FreqYear.m
"""


def ERB2Freq(ERBrate):
    if ERBrate is None:
        raise ValueError("ERB rate must be provided.")

    cf = (10 ** (ERBrate / 21.4) - 1) / 4.37 * 1000
    ERBwidth = 24.7 * (4.37 * cf / 1000 + 1)

    return cf, ERBwidth


