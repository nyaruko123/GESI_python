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
    GCparam = SetHearingLoss(GCparam)  # 设置听力损失参数

    if GCresp is None:
        print('--- GCFBv23_HearingLoss: 设置默认听力损失参数并返回 ---')
        return GCparam

    GCparam['HLoss']['CompressionHealth_InitVal'] = GCparam['HLoss']['CompressionHealth']  # 备份初始值

    LenFag = len(GCparam['HLoss']['FaudgramList'])
    for nFag in range(LenFag):
        Fr1query = GCparam['HLoss']['FaudgramList'][nFag]

        # **确保 HL0_PinCochleadB 是标量**
        HL0_PinCochleadB = HL2PinCochlea(Fr1query, 0)
        HL0_PinCochleadB = float(np.asarray(HL0_PinCochleadB).flatten()[0])
        
        print(f"[DEBUG] HL0_PinCochleadB = {HL0_PinCochleadB}, type = {type(HL0_PinCochleadB)}")

        CompressionHealth = GCparam['HLoss']['CompressionHealth'][nFag]
        _, HL0_IOfuncdB_CH1, _ = GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, 1, HL0_PinCochleadB)
        PindB_ACTreduction = GCFBv23_AsymFuncInOut_InvIOfunc(GCparam, GCresp, Fr1query, CompressionHealth, HL0_IOfuncdB_CH1)

        PinLossdB_ACT = PindB_ACTreduction - HL0_PinCochleadB
        PinLossdB_ACT_Init = PinLossdB_ACT
        PinLossdB_PAS = max(GCparam['HLoss']['HearingLeveldB'][nFag] - PinLossdB_ACT, 0)

        GCparam['HLoss']['CompressionHealth'][nFag] = CompressionHealth

        # **确保 HLval_PinCochleadB 也是标量**
        HLval_PinCochleadB = HL2PinCochlea(Fr1query, 0) + GCparam['HLoss']['HearingLeveldB'][nFag]
        HLval_PinCochleadB = float(np.asarray(HLval_PinCochleadB).flatten()[0])

        print(f"[DEBUG] HLval_PinCochleadB = {HLval_PinCochleadB}, type = {type(HLval_PinCochleadB)}")

        _, HLval_IOfuncdB_CHval, _ = GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, CompressionHealth, HLval_PinCochleadB)
        GCparam['HLoss']['AFgainCmpnstdB'][nFag] = HLval_IOfuncdB_CHval

    return GCparam

def SetHearingLoss(GCparam):
    GCparam['HLoss']['FaudgramList'] = [125, 250, 500, 1000, 2000, 4000, 8000]
    LenFag = len(GCparam['HLoss']['FaudgramList'])
    
    # 听力损失类型逻辑扩展
    if 'Type' not in GCparam['HLoss'] or not GCparam['HLoss']['Type']:
        GCparam['HLoss']['Type'] = 'NH_NormalHearing'

    hl_type = GCparam['HLoss']['Type']
    
    if hl_type.startswith('NH'):
        # 正常听力配置（原有逻辑）
        GCparam['HLoss']['HearingLeveldB'] = np.zeros(LenFag)
    elif hl_type == 'HL0':
        # 平坦型轻度听力损失（全频段20dB）
        GCparam['HLoss']['HearingLeveldB'] = [20,20,20,20,20,20,20]
    elif hl_type == 'HL1':
        # 低频陡降型（125-1kHz逐步下降）
        GCparam['HLoss']['HearingLeveldB'] = [40,35,30,25,20,20,20]
    elif hl_type == 'HL2':
        # 高频缓降型（2k-8kHz逐步上升）
        GCparam['HLoss']['HearingLeveldB'] = [20,20,20,20,30,40,50]
    elif hl_type == 'HL3':
        # 噪声性聋典型曲线（4k凹陷）
        GCparam['HLoss']['HearingLeveldB'] = [25,25,30,35,45,55,50]
    elif hl_type == 'HL4':
        # 老年性聋对称型（高频陡降）
        GCparam['HLoss']['HearingLeveldB'] = [30,30,35,40,55,65,70]
    elif hl_type == 'HL5':
        # 传导性聋平坦型（骨导正常）
        GCparam['HLoss']['HearingLeveldB'] = [45,45,45,45,45,45,45]
    elif hl_type == 'HL6':
        # 混合性聋（低频传导+高频感音）
        GCparam['HLoss']['HearingLeveldB'] = [35,35,40,45,55,65,60]
    elif hl_type == 'HL7':
        # 全频重度损失（平坦型）
        GCparam['HLoss']['HearingLeveldB'] = [70,70,70,70,70,70,70]
    elif hl_type == 'HL8':
        # 极重度高频损失（陡降型）
        GCparam['HLoss']['HearingLeveldB'] = [20,20,25,30,80,90,100]
    else:
        # 无效类型回退默认
        GCparam['HLoss']['Type'] = 'NH_NormalHearing'
        GCparam['HLoss']['HearingLeveldB'] = np.zeros(LenFag)

    # 共用参数初始化（保留原有逻辑）
    GCparam['HLoss']['PinLossdB_ACT'] = np.zeros(LenFag)
    GCparam['HLoss']['PinLossdB_PAS'] = np.zeros(LenFag)
    GCparam['HLoss']['IOfuncLossdB_PAS'] = np.zeros(LenFag)
    
    if 'AFgainCmpnstdB' not in GCparam['HLoss']:
        GCparam['HLoss']['AFgainCmpnstdB'] = np.zeros(LenFag)
    
    if 'CompressionHealth' not in GCparam['HLoss']:
        GCparam['HLoss']['CompressionHealth'] = np.ones(LenFag)

    return GCparam

