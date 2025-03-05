import numpy as np

def SetFrame4TimeSequence(Snd, LenWin, LenShift=None):
    if LenShift is None:
        LenShift = LenWin // 2

    IntDivFrame = LenWin / LenShift

    if IntDivFrame % 1 != 0 or LenWin % 2 != 0:
        raise ValueError(f'LenWin = {LenWin}, LenShift = {LenShift}, Ratio = {IntDivFrame:.4f} <-- should be integer value. LenWin must be even number. LenShift must be LenWin/Integer value.')

    Snd1 = np.concatenate((np.zeros(LenWin // 2), Snd.flatten(), np.zeros(LenWin // 2)))
    LenSnd1 = len(Snd1)
    NumFrame1 = int(np.ceil(LenSnd1 / LenWin))
    nlim = LenWin * NumFrame1
    Snd1 = np.concatenate((Snd1[:min(nlim, LenSnd1)], np.zeros(nlim - LenSnd1)))
    LenSnd1 = len(Snd1)

    NumFrameAll = (NumFrame1 - 1) * IntDivFrame + 1
    SndFrame = np.zeros((LenWin, int(NumFrameAll)))
    NumSmplPnt = np.zeros(int(NumFrameAll), dtype=int)

    for nid in range(int(IntDivFrame)):
        NumFrame2 = NumFrame1 - (nid > 0)
        nSnd = LenShift * nid + np.arange(NumFrame2 * LenWin)
        Snd2 = Snd1[nSnd]
        Mtrx = Snd2.reshape(LenWin, NumFrame2)
        num = np.arange(nid + 1, NumFrameAll + 1, IntDivFrame, dtype=int)
        nIndx = (num - 1) * LenShift
        SndFrame[:, num - 1] = Mtrx
        NumSmplPnt[num - 1] = nIndx

    nValidNumSmplPnt = NumSmplPnt <= len(Snd)
    SndFrame = SndFrame[:, nValidNumSmplPnt]
    NumSmplPnt = NumSmplPnt[nValidNumSmplPnt]

    return SndFrame, NumSmplPnt
