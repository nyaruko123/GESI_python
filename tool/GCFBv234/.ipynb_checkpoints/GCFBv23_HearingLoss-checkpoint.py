"""
Function [GCparam] = GCFBv23_HearingLoss(GCparam, GCresp)
    INPUT:
        GCparam: Parameters for Gammachirp filterbank
        GCresp.Fr1: Frequency response
    OUTPUT:
        GCparam.HLoss: Hearing loss parameters including PinLossdB_PAS, PinLossdB_ACT, FB_PinLossdB_PAS, etc.

    For setting hearing loss, see 'function GCparam = SetHearingLoss(GCparam)'
    Basically the same as WHIS setting
    GCparam.HLoss.FaudgramList = [125, 250, 500, 1000, 2000, 4000, 8000]
    GCparam.HLoss.HearingLeveldB: Hearing level in dB
    GCparam.HLoss.SwType: Selection of preset hearing loss
        NH: Normal hearing (= Hearing level 0 dB)
        HL0: Manual setting of GCparam.HLoss.HearingLeveldB -- necessary
        HL1: Example 1
        HL2: Tsuiki 2002 80yr (average of 80 years old listener)
        HL3: ISO7029 70yr male
        HL4: ISO7029 70yr female
        HL5: ISO7029 60yr male
        HL6: ISO7029 60yr female
        HL7: Example of otosclerosis (textbook: Understanding audiogram p.47)
        HL8: Example of noise induced HL (textbook: Understanding audiogram p.63)

Note: 21 May 2020
    Assumption: pGC is always the same for both NH and HI listeners. Only HP-AF differs.
"""

import numpy as np
from .tool.Freq2ERB import Freq2ERB
from .tool.HL2PinCochlea import HL2PinCochlea
from .GCFBv23_AsymFuncInOut import GCFBv23_AsymFuncInOut
from .GCFBv23_AsymFuncInOut_InvIOfunc import GCFBv23_AsymFuncInOut_InvIOfunc
def GCFBv23_HearingLoss(GCparam, GCresp):
    GCparam = SetHearingLoss(GCparam)  # Set hearing loss parameters

    if GCresp is None:
        print('--- GCFBv23_HearingLoss: Setting default hearing loss parameter and return. ---')
        return GCparam

    # Setting parameters of hearing loss
    GCparam['HLoss']['CompressionHealth_InitVal'] = GCparam['HLoss']['CompressionHealth']  # Keep initial value

    LenFag = len(GCparam['HLoss']['FaudgramList'])
    for nFag in range(LenFag):
        Fr1query = GCparam['HLoss']['FaudgramList'][nFag]
        HL0_PinCochleadB = HL2PinCochlea(Fr1query, 0)  # Convert to cochlear input level
        CompressionHealth = GCparam['HLoss']['CompressionHealth'][nFag]
        _, HL0_IOfuncdB_CH1 = GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, 1, HL0_PinCochleadB)
        PindB_ACTreduction = GCFBv23_AsymFuncInOut_InvIOfunc(GCparam, GCresp, Fr1query, CompressionHealth, HL0_IOfuncdB_CH1)

        PinLossdB_ACT = PindB_ACTreduction - HL0_PinCochleadB
        PinLossdB_ACT_Init = PinLossdB_ACT
        PinLossdB_PAS = max(GCparam['HLoss']['HearingLeveldB'][nFag] - PinLossdB_ACT, 0)

        if PinLossdB_PAS < np.finfo(float).eps * 10**4:
            PinLossdB_ACT = GCparam['HLoss']['HearingLeveldB'][nFag] - PinLossdB_PAS
            CmprsHlthList = np.arange(1, -0.1, -0.1)
            PinLossdB_ACT4Cmpnst = []

            for CmprsHlth in CmprsHlthList:
                PindB_CmprsHlthVal_Inv = GCFBv23_AsymFuncInOut_InvIOfunc(GCparam, GCresp, Fr1query, CmprsHlth, HL0_IOfuncdB_CH1)
                PinLossdB_ACT4Cmpnst.append(PindB_CmprsHlthVal_Inv - HL0_PinCochleadB)

            CompressionHealth = np.interp(PinLossdB_ACT, PinLossdB_ACT4Cmpnst, CmprsHlthList, left=np.nan, right=np.nan)

            if np.isnan(CompressionHealth):
                raise ValueError('Error in CompressionHealth recalculation')

            PindB_ACTreduction = GCFBv23_AsymFuncInOut_InvIOfunc(GCparam, GCresp, Fr1query, CompressionHealth, HL0_IOfuncdB_CH1)
            PinLossdB_ACT = PindB_ACTreduction - HL0_PinCochleadB
            PinLossdB_PAS = GCparam['HLoss']['HearingLeveldB'][nFag] - PinLossdB_ACT

            if abs(GCparam['HLoss']['CompressionHealth_InitVal'][nFag] - CompressionHealth) > np.finfo(float).eps:
                print(f'Compensated GCparam.HLoss.CompressionHealth ({Fr1query} Hz): {GCparam["HLoss"]["CompressionHealth_InitVal"][nFag]} --> {CompressionHealth}')

        ErrorACTPAS = GCparam['HLoss']['HearingLeveldB'][nFag] - (PinLossdB_PAS + PinLossdB_ACT)
        if abs(ErrorACTPAS) > np.finfo(float).eps * 100:
            print(ErrorACTPAS, GCparam['HLoss']['HearingLeveldB'][nFag], PinLossdB_ACT, PinLossdB_PAS)
            if not GCparam['HLoss']['Type'].startswith('NH'):
                raise ValueError('Error in HL_total = HL_ACT + HL_PAS')

        GCparam['HLoss']['CompressionHealth'][nFag] = CompressionHealth
        HLval_PinCochleadB = HL2PinCochlea(Fr1query, 0) + GCparam['HLoss']['HearingLeveldB'][nFag]
        _, HLval_IOfuncdB_CHval = GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, CompressionHealth, HLval_PinCochleadB)
        GCparam['HLoss']['AFgainCmpnstdB'][nFag] = HLval_IOfuncdB_CHval

    NHgainCmpnstBiasdB = [0, 0, 0, 0, 0, 0, 0]
    GCparam['HLoss']['AFgainCmpnstdB'] = np.add(GCparam['HLoss']['AFgainCmpnstdB'], NHgainCmpnstBiasdB)
    GCparam['HLoss']['HLval_PinCochleadB'] = HLval_PinCochleadB
    GCparam['HLoss']['PinLossdB_ACT'] = PinLossdB_ACT
    GCparam['HLoss']['PinLossdB_PAS'] = PinLossdB_PAS
    GCparam['HLoss']['PinLossdB_ACT_Init'] = PinLossdB_ACT_Init

    ERBrateFag = Freq2ERB(GCparam['HLoss']['FaudgramList'])
    ERBrateFr1 = Freq2ERB(GCresp['Fr1'])
    GCparam['HLoss']['FB_Fr1'] = GCresp['Fr1']
    GCparam['HLoss']['FB_HearingLeveldB'] = np.interp(ERBrateFr1, ERBrateFag, GCparam['HLoss']['HearingLeveldB'], left='linear', right='extrap')
    GCparam['HLoss']['FB_HLval_PinCochleadB'] = np.interp(ERBrateFr1, ERBrateFag, GCparam['HLoss']['HLval_PinCochleadB'], left='linear', right='extrap')
    GCparam['HLoss']['FB_PinLossdB_PAS'] = np.interp(ERBrateFr1, ERBrateFag, GCparam['HLoss']['PinLossdB_PAS'], left='linear', right='extrap')
    GCparam['HLoss']['FB_PinLossdB_ACT'] = np.interp(ERBrateFr1, ERBrateFag, GCparam['HLoss']['PinLossdB_ACT'], left='linear', right='extrap')
    GCparam['HLoss']['FB_CompressionHealth'] = np.clip(np.interp(ERBrateFr1, ERBrateFag, GCparam['HLoss']['CompressionHealth'], left='linear', right='extrap'), 0, 1)
    GCparam['HLoss']['FB_AFgainCmpnstdB'] = np.interp(ERBrateFr1, ERBrateFag, GCparam['HLoss']['AFgainCmpnstdB'], left='linear', right='extrap')

    return GCparam

def SetHearingLoss(GCparam):
    GCparam['HLoss']['FaudgramList'] = [125, 250, 500, 1000, 2000, 4000, 8000]
    LenFag = len(GCparam['HLoss']['FaudgramList'])

    if 'Type' not in GCparam['HLoss'] or not GCparam['HLoss']['Type'] or GCparam['HLoss']['Type'].startswith('NH'):
        GCparam['HLoss']['Type'] = 'NH_NormalHearing'
        GCparam['HLoss']['HearingLeveldB'] = np.zeros(LenFag)
        GCparam['HLoss']['PinLossdB_ACT'] = np.zeros(LenFag)
        GCparam['HLoss']['PinLossdB_PAS'] = np.zeros(LenFag)
        GCparam['HLoss']['IOfuncLossdB_PAS'] = np.zeros(LenFag)
        if 'CompressionHealth' not in GCparam['HLoss']:
            GCparam['HLoss']['CompressionHealth'] = np.ones(LenFag)

    elif GCparam['HLoss']['Type'].startswith('HL'):
        if 'CompressionHealth' not in GCparam['HLoss']:
            GCparam['HLoss']['CompressionHealth'] = 0.5 * np.ones(LenFag)

        NumHL = int(GCparam['HLoss']['Type'][2:4].strip()) if len(GCparam['HLoss']['Type']) > 3 else int(GCparam['HLoss']['Type'][2])
        GCparam['HLoss']['SwType'] = NumHL

        if GCparam['HLoss']['SwType'] == 0:
            GCparam['HLoss']['Type'] = 'HLval_ManualSet'
            if len(GCparam['HLoss']['HearingLeveldB']) < LenFag:
                raise ValueError('Set GCparam.HLoss.HearingLeveldB at FaudgramList in advance.')
            if any(h < 0 for h in GCparam['HLoss']['HearingLeveldB']):
                raise ValueError('GCparam.HLoss.HearingLeveldB must not be negative.')

        elif GCparam['HLoss']['SwType'] == 1:
            GCparam['HLoss']['Type'] = 'HL1_Example'
            GCparam['HLoss']['HearingLeveldB'] = [10, 4, 10, 13, 48, 58, 79]

        elif GCparam['HLoss']['SwType'] == 2:
            GCparam['HLoss']['Type'] = 'HL2_Tsuiki2002_80yr'
            GCparam['HLoss']['HearingLeveldB'] = [23.5, 24.3, 26.8, 27.9, 32.9, 48.3, 68.5]

        elif GCparam['HLoss']['SwType'] == 3:
            GCparam['HLoss']['Type'] = 'HL3_ISO7029_70yr_male'
            GCparam['HLoss']['HearingLeveldB'] = [8, 8, 9, 10, 19, 43, 59]

        elif GCparam['HLoss']['SwType'] == 4:
            GCparam['HLoss']['Type'] = 'HL4_ISO7029_70yr_female'
            GCparam['HLoss']['HearingLeveldB'] = [8, 8, 9, 10, 16, 24, 41]

        elif GCparam['HLoss']['SwType'] == 5:
            GCparam['HLoss']['Type'] = 'HL5_ISO7029_60yr_male'
            GCparam['HLoss']['HearingLeveldB'] = [5, 5, 6, 7, 12, 28, 39]

        elif GCparam['HLoss']['SwType'] == 6:
            GCparam['HLoss']['Type'] = 'HL6_ISO7029_60yr_female'
            GCparam['HLoss']['HearingLeveldB'] = [5, 5, 6, 7, 11, 16, 26]

        elif GCparam['HLoss']['SwType'] == 7:
            GCparam['HLoss']['Type'] = 'HL7_Example_Otosclerosis'
            GCparam['HLoss']['HearingLeveldB'] = [50, 55, 50, 50, 40, 25, 20]

        elif GCparam['HLoss']['SwType'] == 8:
            GCparam['HLoss']['Type'] = 'HL8_Example_NoiseInduced'
            GCparam['HLoss']['HearingLeveldB'] = [15, 10, 15, 10, 10, 40, 20]

        else:
            raise ValueError('Specify GCparam.HLoss.Type (HL0, HL1, HL2, ....) properly.')

    else:
        raise ValueError('Specify GCparam.HLoss.Type (NH, HL0, HL1, HL2, ....) properly.')

    if len(GCparam['HLoss']['CompressionHealth']) == 1:
        GCparam['HLoss']['CompressionHealth'] = GCparam['HLoss']['CompressionHealth'] * np.ones(LenFag)

    return GCparam


