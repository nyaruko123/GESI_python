import numpy as np
from time import time
from .tool.CmprsGCFrsp import CmprsGCFrsp
from .GCFBv23_AsymFuncInOut import GCFBv23_AsymFuncInOut
from .tool.SetFrame4TimeSequence import SetFrame4TimeSequence

def GCFBv23_FrameBase(pGCsmpl, scGCsmpl, GCparam, GCresp):
    ExpDecayFrame = GCparam['LvlEst']['ExpDecayVal'] ** GCparam['DynHPAF']['LenShift']
    NumCh, LenSnd = pGCsmpl.shape
    print('--- Frame base processing ---')
    Tstart = time()

    NfrqRsl = 1024 * 2
    c2val_CmprsHlth = GCparam['HLoss']['FB_CompressionHealth'] * GCparam['LvlEst']['c2']
    scGCresp = CmprsGCFrsp(GCparam['Fr1'], GCparam['fs'], GCparam['n'], GCresp['b1val'], GCresp['c1val'],
                           GCparam['LvlEst']['frat'], GCparam['LvlEst']['b2'], c2val_CmprsHlth, NfrqRsl)

    for nch in range(NumCh):
        pGCframeMtrx, nSmplPt = SetFrame4TimeSequence(
            pGCsmpl[nch, :], GCparam['DynHPAF']['LenFrame'], GCparam['DynHPAF']['LenShift'])
        LenWin, LenFrame = pGCframeMtrx.shape

        scGCframeMtrx, nSmplPt = SetFrame4TimeSequence(
            scGCsmpl[nch, :], GCparam['DynHPAF']['LenFrame'], GCparam['DynHPAF']['LenShift'])

        if nch == 0:
            pGCframe = np.zeros((NumCh, LenFrame))
            LvldBframe = np.zeros((NumCh, LenFrame))
            fratFrame = np.zeros((NumCh, LenFrame))
            AsymFuncGain = np.zeros((NumCh, LenFrame))
            dcGCframe = np.zeros((NumCh, LenFrame))

        pGCframe[nch, :LenFrame] = np.sqrt(GCparam['DynHPAF']['ValWin'] @ (pGCframeMtrx ** 2))
        scGCframe = np.sqrt(GCparam['DynHPAF']['ValWin'] @ (scGCframeMtrx ** 2))

        LvlLin1FrameMtrx, nSmplPt = SetFrame4TimeSequence(
            pGCsmpl[GCparam['LvlEst']['NchLvlEst'][nch], :], GCparam['DynHPAF']['LenFrame'], GCparam['DynHPAF']['LenShift'])
        LvlLin1Frame = np.sqrt(GCparam['DynHPAF']['ValWin'] @ (LvlLin1FrameMtrx ** 2))

        LvlLin2FrameMtrx, nSmplPt = SetFrame4TimeSequence(
            scGCsmpl[GCparam['LvlEst']['NchLvlEst'][nch], :], GCparam['DynHPAF']['LenFrame'], GCparam['DynHPAF']['LenShift'])
        LvlLin2Frame = np.sqrt(GCparam['DynHPAF']['ValWin'] @ (LvlLin2FrameMtrx ** 2))

        for nFrame in range(LenFrame - 1):
            LvlLin1Frame[nFrame + 1] = max(LvlLin1Frame[nFrame + 1], LvlLin1Frame[nFrame] * ExpDecayFrame)
            LvlLin2Frame[nFrame + 1] = max(LvlLin2Frame[nFrame + 1], LvlLin2Frame[nFrame] * ExpDecayFrame)

        LvlLinTtlFrame = (GCparam['LvlEst']['Weight'] *
                          GCparam['LvlEst']['LvlLinRef'] * (LvlLin1Frame / GCparam['LvlEst']['LvlLinRef']) ** GCparam['LvlEst']['Pwr'][0] +
                          (1 - GCparam['LvlEst']['Weight']) *
                          GCparam['LvlEst']['LvlLinRef'] * (LvlLin2Frame / GCparam['LvlEst']['LvlLinRef']) ** GCparam['LvlEst']['Pwr'][1])

        CmpnstHalfWaveRectify = -3
        LvldBframe[nch, :LenFrame] = 20 * np.log10(np.maximum(LvlLinTtlFrame, GCparam['LvlEst']['LvlLinMinLim'])) + \
                                     GCparam['LvlEst']['RMStoSPLdB'] + CmpnstHalfWaveRectify

        AFoutdB, IOfuncdB, GCparam = GCFBv23_AsymFuncInOut(GCparam, GCresp, GCparam['Fr1'][nch],
                                                           GCparam['HLoss']['FB_CompressionHealth'][nch], LvldBframe[nch, :])
        AsymFuncGain[nch, :LenFrame] = 10 ** (AFoutdB / 20)

        scGCframe1 = scGCresp['NormFctFp2'][nch] * scGCframe

        dcGCframe[nch, :LenFrame] = AsymFuncGain[nch, :] * scGCframe1

        if nch == 0 or nch % 50 == 0:
            print(f'Frame-based HP-AF: ch #{nch + 1} / #{NumCh}. elapsed time = {round(time() - Tstart, 1)} (sec)')

    GCresp['LvldBframe'] = LvldBframe
    GCresp['pGCframe'] = pGCframe
    GCresp['scGCframe'] = scGCframe
    GCresp['fratFrame'] = fratFrame
    GCresp['AsymFuncGain'] = AsymFuncGain

    return dcGCframe, GCresp

# 未定义的函数或变量: CmprsGCFrsp, SetFrame4TimeSequence, GCFBv23_AsymFuncInOut
