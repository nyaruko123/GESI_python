"""
Function [dcGCout, scGCsmpl, GCparam, GCresp] = GCFBv234(SndIn, GCparam)
    INPUT:
        SndIn: Input Sound
        GCparam: Gammachirp parameters
            GCparam.fs: Sampling rate (48000)
            GCparam.NumCh: Number of Channels (100)
            GCparam.FRange: Frequency Range of GCFB: default [100 6000]
            GCparam.HLoss: Struct of Hearing Loss setting

    OUTPUT:
        dcGCout: Dynamic Compressive GammaChirp Filter Output
        scGCsmpl: Static Compressive GammaChirp Filter Output
        Ppgc: Power at the output of passive GC
        GCparam: GCparam values
        GCresp: GC response result
        pGCframe: pGC frame -- empty when sample-by-sample
        scGCframe: scGC frame -- empty when sample-by-sample

Notes:
    1) This version is completely different from GCFB v.1.04 (obsolete).
       We introduced the "compressive gammachirp" to accommodate both the
       psychoacoustical simultaneous masking and the compressive
       characteristics (Irino and Patterson, 2001). The parameters were
       determined from large dataset (See Patterson, Unoki, and Irino, 2003.)

References:
    Irino, T. and Unoki, M.: IEEE ICASSP98, pp.3653-3656, May 1998.
    Irino, T. and Patterson, R.D.: JASA, Vol.101, pp.412-419, 1997.
    Irino, T. and Patterson, R.D.: JASA, Vol.109, pp.2008-2022, 2001.
    Patterson, R.D., Unoki, M. and Irino, T.: JASA, Vol.114, pp.1529-1542, 2003.
    Irino, T. and Patterson, R.D.: IEEE Trans.ASLP, Vol.14, Nov. 2006.

Level setting: See Eqlz2MeddisHCLevel
    rms(s(t)) == sqrt(mean(s**2)) == 1   --> 30 dB SPL
    rms(s(t)) == sqrt(mean(s**2)) == 10  --> 50 dB SPL
    rms(s(t)) == sqrt(mean(s**2)) == 100 --> 70 dB SPL
"""

from .GCFBv23_SetParam import GCFBv23_SetParam
from .tool.MkFilterField2Cochlea import MkFilterField2Cochlea
from .tool.MakeAsymCmpFiltersV2 import MakeAsymCmpFiltersV2
from .tool.GammaChirp import GammaChirp
from .tool.CmprsGCFrsp import CmprsGCFrsp
from .GCFBv23_SampleBase import GCFBv23_SampleBase
from .GCFBv23_FrameBase import GCFBv23_FrameBase

def GCFBv234(SndIn, GCparam):
    import numpy as np
    from time import time

    print('-------------------- GCFBv234 --------------------')

    if SndIn.ndim != 1:
        raise ValueError('Check SndIn. It should be 1 ch (Monaural) and a single row vector.')

    GCparam, GCresp = GCFBv23_SetParam(GCparam)
    GCparam.setdefault('GainRefdB', 50)  # 如果没有设置，则默认设为50 dB

    fs = GCparam['fs']
    NumCh = GCparam['NumCh']
    Tstart = time()

    if GCparam['OutMidCrct'].upper() != 'NO':
        print(f"*** Outer/Middle Ear correction (minimum phase) : {GCparam['OutMidCrct']} ***")
        CmpnstOutMid, ParamF2C = MkFilterField2Cochlea(GCparam['OutMidCrct'], fs, 1)
        Snd = np.convolve(SndIn, CmpnstOutMid, mode='same')
        GCparam['Field2Cochlea'] = ParamF2C
    else:
        GCparam['Field2Cochlea'] = 'No Outer/Middle Ear correction'
        print(f"*** {GCparam['Field2Cochlea']} ***")
        Snd = np.asarray(SndIn, dtype=np.float64)

    print('*** Gammmachirp Calculation ***')
    
    if 'Fr2' not in GCresp or GCresp['Fr2'] is None or len(GCresp['Fr2']) == 0:
        raise ValueError("GCresp['Fr2'] is empty or not properly initialized!")
    
    # 正确解包 `MakeAsymCmpFiltersV2` 返回的多个值
    ACFcoefFixed, ACFcoefConv = MakeAsymCmpFiltersV2(fs, GCresp['Fr2'], GCresp['b2val'], GCresp['c2val'])
    
    print(f"DEBUG: ACFcoefFixed type = {type(ACFcoefFixed)}")

    pGCsmpl = np.zeros((NumCh, len(SndIn)))
    scGCsmpl = np.zeros((NumCh, len(SndIn)))

    print('--- Channel-by-channel processing of static filter ---')
    for nch in range(NumCh):
        pgc = GammaChirp(GCresp['Fr1'][nch], fs, GCparam['n'], GCresp['b1val'][nch], GCresp['c1val'][nch], 0, '', 'peak')
        if isinstance(pgc, tuple):
            pgc = pgc[0]  # 只取第一个返回值

        pgc = np.asarray(pgc, dtype=np.float64)
        if pgc.ndim > 1:
            pgc = pgc.ravel()
        
        Snd = np.asarray(Snd, dtype=np.float64)
        if Snd.size == 0 or pgc.size == 0:
            raise ValueError(f"Error: Empty array detected! Snd.size={Snd.size}, pgc.size={pgc.size}")
        
        pGCsmpl[nch, :] = np.convolve(Snd, pgc, mode='same')
        GCsmpl1 = pGCsmpl[nch, :]
        for Nfilt in range(4):
            GCsmpl1 = np.convolve(GCsmpl1, ACFcoefFixed['bz'][nch, :, Nfilt], mode='same')
        scGCsmpl[nch, :] = GCsmpl1
        
        if nch == 0 or (nch + 1) % 50 == 0:
            elapsed_time = time() - Tstart
            print(f"Static (Fixed) Compressive-Gammachirp: ch #{nch + 1} / #{NumCh}. elapsed time = {elapsed_time:.1f} sec")
    
    if GCparam['Ctrl'].startswith('dyn'):
        if GCparam['DynHPAF']['StrPrc'].startswith('sample'):
            dcGCsmpl, GCresp = GCFBv23_SampleBase(pGCsmpl, scGCsmpl, GCparam, GCresp)
            cGCout = dcGCsmpl
        elif GCparam['DynHPAF']['StrPrc'].startswith('frame'):
            dcGCframe, GCresp = GCFBv23_FrameBase(pGCsmpl, scGCsmpl, GCparam, GCresp)
            cGCout = dcGCframe
        else:
            raise ValueError('Specify "GCparam.DynHPAF.StrPrc" properly: "sample" or "frame"')
    elif GCparam['Ctrl'].startswith('sta'):
        cGCout = scGCsmpl
    else:
        raise ValueError('Specify GCparam.Ctrl properly')
    
    LenOut = cGCout.shape[1]
    if isinstance(GCparam['GainRefdB'], (int, float)):
        cGCRef = CmprsGCFrsp(GCresp['Fr1'], fs, GCparam['n'], GCresp['b1val'], GCresp['c1val'], GCresp['frat'], GCresp['b2val'], GCresp['c2val'])
        GCresp['GainFactor'] = 10**(GCparam['GainCmpnstdB'] / 20) * cGCRef['NormFctFp2']
        GCresp['cGCRef'] = cGCRef
        # 修改：将 GainFactor 转换为列向量，使乘法广播正确
        dcGCout = (GCresp['GainFactor'][:, np.newaxis]) * cGCout
    elif GCparam['GainRefdB'] == 'NormIOfunc':
        GainFactor = 10**(-(GCparam['HLoss']['FB_AFgainCmpnstdB']) / 20)
        dcGCout = (GainFactor * np.ones((1, LenOut))) * cGCout
    else:
        raise ValueError('Set GCparam.GainRefdB properly')
    
    return dcGCout, scGCsmpl, GCparam, GCresp

# Note: The following functions are placeholders and need to be implemented:
# - StartupGCFB
# - GCFBv23_SetParam
# - MkFilterField2Cochlea
# - MakeAsymCmpFiltersV2
# - GammaChirp
# - filter
# - GCFBv23_SampleBase
# - GCFBv23_FrameBase
# - CmprsGCFrsp
