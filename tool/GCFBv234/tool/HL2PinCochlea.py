"""
Function [PinCchdB] = HL2PinCochlea(freq, HLdB)
    INPUT:
        freq: Frequency
        HLdB: Hearing Level in dB
    OUTPUT:
        PinCchdB: Cochlear Input level in dB
"""
from .HL2SPL import HL2SPL
from .TransFuncMiddleEar_Moore16 import TransFuncMiddleEar_Moore16

def HL2PinCochlea(freq, HLdB):
    SPLdB = HL2SPL(freq, HLdB)
    _, FrspMEdB = TransFuncMiddleEar_Moore16(freq)
    PinCchdB = SPLdB + FrspMEdB
    return PinCchdB


