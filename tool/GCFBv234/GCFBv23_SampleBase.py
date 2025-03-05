import numpy as np
from time import time
from .tool.MakeAsymCmpFiltersV2 import MakeAsymCmpFiltersV2
from .tool.ACFilterBank import ACFilterBank


def GCFBv23_SampleBase(pGCsmpl, scGCsmpl, GCparam, GCresp):
    """
    Perform sample-by-sample processing for the Dynamic Compressive-Gammachirp filter bank.

    INPUT:
        pGCsmpl: passive GC sample
        scGCsmpl: static cGC sample
        GCparam: Gammachirp parameters
        GCresp: Gammachirp response

    OUTPUT:
        cGCsmpl: frame level output of dcGC-FB
        GCresp: updated Gammachirp response
    """

    fs = GCparam.fs
    NumCh, LenSnd = pGCsmpl.shape
    nDisp = LenSnd // 10  # display 10 times per sound
    cGCsmpl = np.zeros((NumCh, LenSnd))
    GCresp.Fr2 = np.zeros((NumCh, LenSnd))
    GCresp.fratVal = np.zeros((NumCh, LenSnd))
    LvldB = np.zeros((NumCh, LenSnd))
    LvlLinPrev = np.zeros((NumCh, 2))

    print('--- Sample base (sample-by-sample) processing ---')
    Tstart = time()

    for nsmpl in range(LenSnd):
        # Level estimation circuit
        LvlLin = np.zeros((NumCh, 2))
        LvlLin[:, 0] = np.maximum(
            np.maximum(pGCsmpl[GCparam.LvlEst.NchLvlEst, nsmpl], 0),
            LvlLinPrev[:, 0] * GCparam.LvlEst.ExpDecayVal
        )
        LvlLin[:, 1] = np.maximum(
            np.maximum(scGCsmpl[GCparam.LvlEst.NchLvlEst, nsmpl], 0),
            LvlLinPrev[:, 1] * GCparam.LvlEst.ExpDecayVal
        )
        LvlLinPrev = LvlLin

        LvlLinTtl = (GCparam.LvlEst.Weight *
                     GCparam.LvlEst.LvlLinRef * (LvlLin[:, 0] / GCparam.LvlEst.LvlLinRef) ** GCparam.LvlEst.Pwr[0] +
                     (1 - GCparam.LvlEst.Weight) *
                     GCparam.LvlEst.LvlLinRef * (LvlLin[:, 1] / GCparam.LvlEst.LvlLinRef) ** GCparam.LvlEst.Pwr[1])

        LvldB[:, nsmpl] = 20 * np.log10(np.maximum(LvlLinTtl, GCparam.LvlEst.LvlLinMinLim)) + GCparam.LvlEst.RMStoSPLdB

        # Signal path
        fratVal = (GCresp.frat0Pc + GCparam.HLoss.FB_CompressionHealth * GCresp.frat1val *
                   (LvldB[:, nsmpl] - GCresp.PcHPAF))
        Fr2Val = GCresp.Fp1 * fratVal

        GCparam.NumUpdateAsymCmp = 1
        if (nsmpl % GCparam.NumUpdateAsymCmp) == 0:
            ACFcoef = MakeAsymCmpFiltersV2(fs, Fr2Val, GCresp.b2val, GCresp.c2val)  # Undefined function

        if nsmpl == 0:
            _, ACFstatus = ACFilterBank(ACFcoef, None)  # Initialization, Undefined function

        SigOut, ACFstatus = ACFilterBank(ACFcoef, ACFstatus, pGCsmpl[:, nsmpl])  # Undefined function
        cGCsmpl[:, nsmpl] = SigOut
        GCresp.Fr2[:, nsmpl] = Fr2Val
        GCresp.fratVal[:, nsmpl] = fratVal

        if nsmpl == 0 or (nsmpl % nDisp) == 0:
            elapsed_time = time() - Tstart
            print(f'Dynamic Compressive-Gammachirp: Time {nsmpl / fs * 1000:.0f}(ms) / {LenSnd / fs * 1000:.0f}(ms). '
                  f'elapsed time = {elapsed_time:.1f} (sec)')

    GCresp.LvldB = LvldB

    return cGCsmpl, GCresp

# 注意：在此代码中，`MakeAsymCmpFiltersV2` 和 `ACFilterBank` 是未定义的函数，需要在使用前定义或导入。
