"""
Function GCparam = GCFBv2xx_SetParam(GCparam)
    INPUT:
        GCparam: Your preset gammachirp parameters
            GCparam.fs: Sampling rate (48000)
            GCparam.NumCh: Number of Channels (75)
            GCparam.FRange: Frequency Range of GCFB [100 6000]
            GCparam.Ctrl: 'sta[tic]'=='fix[ed]' / 'dyn[amic]'=='tim[e-varying]'

    OUTPUT:
        GCparam: GCparam values

Reference:
    Patterson, R.D., Unoki, M. and Irino, T.: JASA, Vol.114, pp.1529-1542, 2003.
"""
import numpy as np
from .tool.EqualFreqScale import EqualFreqScale
from tool.GCFBv234.tool.Freq2ERB import Freq2ERB
from .tool.Fr2Fpeak import Fr2Fpeak
from .GCFBv23_HearingLoss import GCFBv23_HearingLoss

def GCFBv23_SetParam(GCparam):

    # Handling Input Parameters
    GCparam.setdefault('fs', 48000)
    GCparam.setdefault('OutMidCrct', 'ELC')
    GCparam.setdefault('NumCh', 100)
    GCparam.setdefault('FRange', [100, 6000])

    if GCparam['FRange'][1] * 3 > GCparam['fs']:
        print(GCparam)
        print('GCFB may not work properly when max(FreqRange)*3 > fs.')
        print('---> Set fs properly. OR If you wish to continue as is, press RETURN > ')
        input()

    # Gammachirp parameters
    GCparam.setdefault('n', 4)
    GCparam.setdefault('b1', [1.81, 0])
    GCparam.setdefault('c1', [-2.96, 0])
    GCparam.setdefault('frat', [[0.466, 0], [0.0109, 0]])
    GCparam.setdefault('b2', [[2.17, 0], [0, 0]])
    GCparam.setdefault('c2', [[2.20, 0], [0, 0]])
    GCparam.setdefault('Ctrl', 'dynamic')
    GCparam.setdefault('GainCmpnstdB', -1)

    if GCparam['Ctrl'].startswith('fix'):
        GCparam['Ctrl'] = 'static'
    if GCparam['Ctrl'].startswith('tim'):
        GCparam['Ctrl'] = 'dynamic'

    if not (GCparam['Ctrl'].startswith('sta') or GCparam['Ctrl'].startswith('dyn') or GCparam['Ctrl'].startswith(
            'lev')):
        raise ValueError('Specify GCparam.Ctrl: "static", "dynamic", or "level(-estimation)"')

    # Parameters for level estimation
    if 'PpgcRef' in GCparam or 'LvlRefdB' in GCparam:
        print('The parameter "GCparam.PpgcRef" is obsolete.')
        print('The parameter "GCparam.LvlRefdB" is obsolete.')
        raise ValueError('Please change it to GCparam.GainRefdB.')

    GCparam.setdefault('GainRefdB', 'NormIOfunc')
    GCparam.setdefault('LeveldBscGCFB', 50)
    GCparam.setdefault('LvlEst', {})
    GCparam['LvlEst'].setdefault('LctERB', 1.5)
    GCparam['LvlEst'].setdefault('DecayHL', 0.5)
    GCparam['LvlEst'].setdefault('b2', GCparam['b2'][0][0])
    GCparam['LvlEst'].setdefault('c2', GCparam['c2'][0][0])
    GCparam['LvlEst'].setdefault('frat', 1.08)
    GCparam['LvlEst'].setdefault('RMStoSPLdB', 30)
    GCparam['LvlEst'].setdefault('Weight', 0.5)
    GCparam['LvlEst'].setdefault('RefdB', 50)
    GCparam['LvlEst'].setdefault('Pwr', [1.5, 0.5])
    GCparam.setdefault('NumUpdateAsymCmp', 1)

    # Sample-by-sample or Frame-base processing
    GCparam.setdefault('DynHPAF', {})
    GCparam['DynHPAF'].setdefault('StrPrc', 'sample-by-sample')

    if GCparam['DynHPAF']['StrPrc'].startswith('frame'):
        GCparam['DynHPAF']['Tframe'] = 0.001
        GCparam['DynHPAF']['Tshift'] = 0.0005
        GCparam['DynHPAF']['LenFrame'] = int(GCparam['DynHPAF']['Tframe'] * GCparam['fs'])
        GCparam['DynHPAF']['LenShift'] = int(GCparam['DynHPAF']['Tshift'] * GCparam['fs'])
        GCparam['DynHPAF']['Tframe'] = GCparam['DynHPAF']['LenFrame'] / GCparam['fs']
        GCparam['DynHPAF']['Tshift'] = GCparam['DynHPAF']['LenShift'] / GCparam['fs']
        GCparam['DynHPAF']['fs'] = 1 / GCparam['DynHPAF']['Tshift']
        GCparam['DynHPAF']['NameWin'] = 'hanning'
        GCparam['DynHPAF']['ValWin'] = np.hanning(GCparam['DynHPAF']['LenFrame'])
        GCparam['DynHPAF']['ValWin'] /= np.sum(GCparam['DynHPAF']['ValWin'])

    # GCresp
    # Placeholder for EqualFreqScale function
    Fr1, ERBrate1 = EqualFreqScale('ERB', GCparam['NumCh'], GCparam['FRange'])
    GCparam['Fr1'] = Fr1
    GCresp = {}
    GCresp['Fr1'] = Fr1
    GCresp['ERBspace1'] = np.mean(np.diff(ERBrate1))
    # Placeholder for Freq2ERB function
    ERBrate, ERBw = Freq2ERB(GCresp['Fr1'])
    ERBrate1kHz, ERBw1kHz = Freq2ERB(1000)
    GCresp['Ef'] = ERBrate / ERBrate1kHz - 1

    OneVec = np.ones(GCparam['NumCh'])
    GCresp['b1val'] = GCparam['b1'][0] * OneVec + GCparam['b1'][1] * GCresp['Ef']
    GCresp['c1val'] = GCparam['c1'][0] * OneVec + GCparam['c1'][1] * GCresp['Ef']
    # Placeholder for Fr2Fpeak function
    GCresp['Fp1'] = Fr2Fpeak(GCparam['n'], GCresp['b1val'], GCresp['c1val'], GCresp['Fr1'])
    GCresp['b2val'] = GCparam['b2'][0][0] * OneVec + GCparam['b2'][0][1] * GCresp['Ef']
    GCresp['c2val'] = GCparam['c2'][0][0] * OneVec + GCparam['c2'][0][1] * GCresp['Ef']

    # New parameters for HPAF
    GCresp['frat0val'] = GCparam['frat'][0][0] * OneVec + GCparam['frat'][0][1] * GCresp['Ef']
    GCresp['frat1val'] = GCparam['frat'][1][0] * OneVec + GCparam['frat'][1][1] * GCresp['Ef']
    GCresp['PcHPAF'] = (1 - GCresp['frat0val']) / GCresp['frat1val']
    GCresp['frat0Pc'] = GCresp['frat0val'] + GCresp['frat1val'] * GCresp['PcHPAF']

    # GC Hearing Loss
    # Placeholder for GCFBv23_HearingLoss function
    GCparam = GCFBv23_HearingLoss(GCparam, GCresp)

    # Set Params Estimation circuit
    ExpDecayVal = np.exp(-1 / (GCparam['LvlEst']['DecayHL'] * GCparam['fs'] / 1000) * np.log(2))
    NchShift = round(GCparam['LvlEst']['LctERB'] / GCresp['ERBspace1'])
    NchLvlEst = np.clip(np.arange(1, GCparam['NumCh'] + 1) + NchShift, 1, GCparam['NumCh'])
    LvlLinMinLim = 10 ** (-GCparam['LvlEst']['RMStoSPLdB'] / 20)
    LvlLinRef = 10 ** ((GCparam['LvlEst']['RefdB'] - GCparam['LvlEst']['RMStoSPLdB']) / 20)

    GCparam['LvlEst']['ExpDecayVal'] = ExpDecayVal
    GCparam['LvlEst']['ERBspace1'] = GCresp['ERBspace1']
    GCparam['LvlEst']['NchShift'] = NchShift
    GCparam['LvlEst']['NchLvlEst'] = NchLvlEst
    GCparam['LvlEst']['LvlLinMinLim'] = LvlLinMinLim
    GCparam['LvlEst']['LvlLinRef'] = LvlLinRef

    return GCparam, GCresp


