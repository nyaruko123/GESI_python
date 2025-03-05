#123
import os
import numpy as np
from scipy.signal import butter, filtfilt, resample
from scipy.io import loadmat, savemat,wavfile
from tool.TimeAlignXcorr import TimeAlignXcorr
from tool.GCFBv234.tool.TaperWindow import TaperWindow
from tool.FilterModFB import FilterModFB
from tool.GCFBv234 import GCFBv234
from tool.F0limit2SSIweight import F0limit2SSIweight
from tool.Metric2Pcorrect_Sigmoid import Metric2Pcorrect_Sigmoid
from tool.GCFBv234.tool.Eqlz2MeddisHCLevel import Eqlz2MeddisHCLevel
from tool.GCFBv234.tool.EqlzGCFB2Rms1at0dB import EqlzGCFB2Rms1at0dB
from tool.world_0_2_4_matlab import Harvest




def GESIv123(SndRef, SndTest, GCparam, GESIparam):
    # Directory of this program
    DirProg = os.path.dirname(__file__)
    NameProg = os.path.splitext(os.path.basename(__file__))[0]
    Result = {'Name': NameProg}
    print(f'##### Start: {Result["Name"]} #####')

    # Making directory and data for keeping GCout
    if 'DirGCout' not in GESIparam:
        GESIparam['DirGCout'] = os.path.expanduser('~/Data/GESI/GCout/')
        if not os.path.exists(GESIparam['DirGCout']):
            os.makedirs(GESIparam['DirGCout'])

    if 'SwPlot' not in GESIparam:
        GESIparam['SwPlot'] = 0

    StrSPLcond = f'_Rms1SPL{int(GESIparam["DigitalRms1SPLdB"])}dB'
    StrHLossCond = f'{GCparam["HLoss"]["Type"]}_'
    if GCparam['HLoss']['Type'] != 'NH':
        if 'CompressionHealth' not in GCparam['HLoss']:
            raise ValueError('GCparam.HLoss.CompressionHealth should be specified.')
        StrHLossCond = f'{GCparam["HLoss"]["Type"]}_{int(GCparam["HLoss"]["CompressionHealth"] * 100)}_'

    # Parameters of GCFB & GESI
    if 'fs' not in GCparam:
        GCparam['fs'] = 48000
    if 'NumCh' not in GCparam:
        GCparam['NumCh'] = 100
    if 'FRange' not in GCparam:
        GCparam['FRange'] = [100, 8000]
    if 'OutMidCrct' not in GCparam:
        GCparam['OutMidCrct'] = 'FreeField'
    if 'HLoss' not in GCparam:
        GCparam['HLoss'] = {'Type': 'NH'}
    GCparam['Ctrl'] = 'dynamic'
    GCparam['DynHPAF'] = {'StrPrc': 'frame-base'}
    GCparam['StrFloor'] = 'NoiseFloor'

    if 'fs' not in GESIparam:
        print('GESIparam.fs was not specified.')
        print(f'Do you use fs = {GCparam["fs"]}Hz ?  [Return to yes] > ')
        input()
        GESIparam['fs'] = GCparam['fs']

    # GCoutRef & GCoutTest are saved for speed up.
    if 'NameSndRef' in GESIparam:
        GESIparam['NameGCoutRef'] = f'GCout_{StrHLossCond}{GESIparam["NameSndRef"]}{StrSPLcond}'
        SwSave = True
    else:
        GESIparam['NameGCoutRef'] = ''
        SwSave = False

    if 'NameSndTest' in GESIparam:
        GESIparam['NameGCoutTest'] = f'GCout_{StrHLossCond}{GESIparam["NameSndTest"]}{StrSPLcond}'
        SwSave = True
    else:
        GESIparam['NameGCoutTest'] = ''
        SwSave = False

    # Sound check
    if SndTest.ndim > 1 or SndRef.ndim > 1:
        raise ValueError('SndTest and SndRef should be monaural.')
    SndRef = SndRef.flatten()
    SndTest = SndTest.flatten()

    # Time alignment
    if len(SndTest) == len(SndRef) and 'SwTimeAlign' not in GESIparam:
        print('GESI: No time alignment:  length(SndRef) == length(SndTest)')
        GESIparam['SwTimeAlign'] = 0
    else:
        if 'SwTimeAlign' not in GESIparam:
            GESIparam['SwTimeAlign'] = 1

        if GESIparam['SwTimeAlign'] == 1:
            SndTest, ParamTA = TimeAlignXcorr(SndTest, SndRef)
            GESIparam['TimeAlign'] = ParamTA
        else:
            raise NotImplementedError('--- Not prepared yet: Another TimeAlign algorithm. ---')

        if 'DurTaperWindow' not in GESIparam:
            GESIparam['DurTaperWindow'] = 0.02
        LenTaper = int(GESIparam['DurTaperWindow'] * GESIparam['fs'])
        Win = TaperWindow(len(SndRef), 'han', LenTaper)

         # 如果 Win 是 tuple 或 list，只取第一个元素
        if isinstance(Win, (tuple, list)):
            Win = Win[0]

        # 转换成 NumPy 数组
        Win = np.asarray(Win, dtype=np.float32)
        SndRef *= Win
        SndTest *= Win

    # GESI params
    if 'Sim' not in GESIparam or 'PowerRatio' not in GESIparam['Sim']:
        GESIparam['Sim'] = {'PowerRatio': 0.6}
        print('GESIparam.Sim.PowerRatio is set to 0.6 (default) -- OK? Return to continue > ')
        input()
    #sui bian ding de 
    Env = np.copy(SndRef)
    ParamMFB = {'fs':GESIparam['fs']}

    _, MFBparam = FilterModFB(Env, ParamMFB)
    if 'weightMFB' not in GESIparam['Sim']:
        LenMFB = len(MFBparam['fc'])
        GESIparam['Sim']['weightMFB'] = np.ones(LenMFB)
    else:
        LenMFB = len(GESIparam['Sim']['weightMFB'])

    if 'SwWeightProhibit' not in GESIparam['Sim']:
        GESIparam['Sim']['SwWeightProhibit'] = 1
    if 'RangeWeightProhibit' not in GESIparam['Sim']:
        GESIparam['Sim']['RangeWeightProhibit'] = 1

    # Sound sampling rate conversion & normalization
    if GCparam['fs'] != GESIparam['fs']:
        SndTest = resample(SndTest, GCparam['fs'], GESIparam['fs'])
        SndRef = resample(SndRef, GCparam['fs'], GESIparam['fs'])

    # Calibrate input level of SndRef
    SndRef, MdsAmpdB = Eqlz2MeddisHCLevel(SndRef, [], GESIparam['DigitalRms1SPLdB'])
    SndTest *= 10**(MdsAmpdB[1] / 20)

    # Analysis by dynamic compressive gammachirp filterbank
    DirNameTest = os.path.join(GESIparam['DirGCout'], f'{GESIparam["NameGCoutTest"]}.mat')
    if not os.path.exists(DirNameTest):
        print('==== GCFB calculation of SndTest (HL or NH) ====')
        GCoutTest, _, GCparamTest = GCFBv234(SndTest, GCparam)
        NumCh, LenFrame = GCoutTest.shape
        GCoutTest = EqlzGCFB2Rms1at0dB(GCoutTest, GCparam['StrFloor'])

        MFBparam['fs'] = GCparamTest['DynHPAF']['fs']
        MFBparam['fcutEnv'] = 150
        MFBparam['bzLPF'], MFBparam['apLPF'] = butter(1, MFBparam['fcutEnv'] / (MFBparam['fs'] / 2))

        GCModEnvTest = np.zeros((NumCh, LenMFB, LenFrame))
        for nch in range(GCparam['NumCh']):
            Env = filtfilt(MFBparam['bzLPF'], MFBparam['apLPF'], GCoutTest[nch, :])
            if nch in [0, 49, 99]:
                print(f'> Modulation Filterbank Analysis: GCFB ch = #{nch + 1}')
            ModEnv, MFBparam = FilterModFB(Env, MFBparam)
            GCModEnvTest[nch, :, :] = ModEnv

        if SwSave:
            print(f'save: {DirNameTest}')
            savemat(DirNameTest, {'GCoutTest': GCoutTest, 'GCparamTest': GCparamTest, 'GCModEnvTest': GCModEnvTest, 'MFBparam': MFBparam})
    else:
        print(f'load: {DirNameTest}')
        data = loadmat(DirNameTest)
        GCoutTest = data['GCoutTest']
        GCparamTest = data['GCparamTest']
        GCModEnvTest = data['GCModEnvTest']
        MFBparam = data['MFBparam']

    DirNameRef = os.path.join(GESIparam['DirGCout'], f'{GESIparam["NameGCoutRef"]}.mat')
    if not os.path.exists(DirNameRef):
        print('==== GCFB calculation of SndRef (always NH) ====')
        GCparam['HLoss']['Type'] = 'NH'
        GCoutRef, _, GCparamRef = GCFBv234(SndRef, GCparam)
        NumCh, LenFrame = GCoutRef.shape
        GCoutRef = EqlzGCFB2Rms1at0dB(GCoutRef, GCparam['StrFloor'])
        tFrame = np.arange(LenFrame) / GCparamRef['DynHPAF']['fs']

        GCModEnvRef = np.zeros((NumCh, LenMFB, LenFrame))
        for nch in range(GCparam['NumCh']):
            Env = filtfilt(MFBparam['bzLPF'], MFBparam['apLPF'], GCoutRef[nch, :])
            if nch in [0, 49, 99]:
                print(f'> Modulation Filterbank Analysis: GCFB ch = #{nch + 1}')
            ModEnv, _ = FilterModFB(Env, MFBparam)
            GCModEnvRef[nch, :, :] = ModEnv

        HarvestRef = Harvest(SndRef, GCparam['fs'])
        F0Frame = np.interp(tFrame, HarvestRef['temporal_positions'], HarvestRef['f0'], left=0, right=0)
        F0MeanRef = np.exp(np.mean(np.log(HarvestRef['f0'][HarvestRef['f0'] > 0])))
        print(f'Fo Mean of Ref sound: {F0MeanRef:.1f} Hz')

        if np.any(np.isnan(F0Frame)):
            raise ValueError('Error in F0Frame.')

        SSIparam = {'SwSSIweight': 2, 'Fr1': GCparamRef['Fr1']}
        SSIweightMtrx = np.zeros((NumCh, LenFrame))

        for nFrame in range(LenFrame):
            SSIparam['F0_limit'] = F0Frame[nFrame]
            SSIweight, SSIparam = F0limit2SSIweight(SSIparam)
            SSIweightMtrx[:, nFrame] = SSIweight / np.mean(SSIweight)

        SSIparam['weight'] = SSIweightMtrx

        if SwSave:
            print(f'save: {DirNameRef}')
            savemat(DirNameRef, {'GCoutRef': GCoutRef, 'GCparamRef': GCparamRef, 'GCModEnvRef': GCModEnvRef, 'SSIparam': SSIparam})
    else:
        print(f'load: {DirNameRef}')
        data = loadmat(DirNameRef)
        GCoutRef = data['GCoutRef']
        GCparamRef = data['GCparamRef']
        GCModEnvRef = data['GCModEnvRef']
        SSIparam = data['SSIparam']

    GESIparam['GCparam'] = {'Ref': GCparamRef, 'Test': GCparamTest}
    GESIparam['MFBparam'] = MFBparam
    GESIparam['SSIparam'] = SSIparam

    if GESIparam['SwPlot'] == 1:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.imshow(6 * GCoutTest, aspect='auto', origin='lower')
        plt.subplot(2, 1, 2)
        plt.imshow(6 * GCoutRef, aspect='auto', origin='lower')
        plt.show()
        plt.figure()
        plt.plot(SSIparam['weight'])
        plt.show()

    # Cosine similarity analysis
    NumCh, LenMFB, _ = GCModEnvRef.shape

    for nrPwr in range(len(GESIparam['Sim']['PowerRatio'])):
        CosSimMtrx = np.zeros((NumCh, LenMFB))
        for nch in range(NumCh):
            weightMFB = GESIparam['Sim']['weightMFB'].copy()
            for nMFB in range(LenMFB):
                ModEnvRef = GCModEnvRef[nch, nMFB, :]
                ModEnvTest = GCModEnvTest[nch, nMFB, :]
                weightGCFB = SSIparam['weight'][nch, :]

                PwrRef = np.sum(ModEnvRef**2)
                PwrTest = np.sum(ModEnvTest**2)
                rPwr = GESIparam['Sim']['PowerRatio'][nrPwr]

                CosSim = np.sum(weightGCFB * ModEnvRef * ModEnvTest) / (PwrRef**rPwr * PwrTest**(1 - rPwr))

                if GESIparam['Sim']['SwWeightProhibit'] == 1:
                    if GCparamRef['Fr1'][nch] < MFBparam['fc'][nMFB] * GESIparam['Sim']['RangeWeightProhibit']:
                        weightMFB[nMFB] = np.nan

                CosSimMtrx[nch, nMFB] = weightMFB[nMFB] * CosSim

        Result['dIntrm'] = {'GESI': CosSimMtrx}
        Result['d'] = {'GESI': np.nanmean(CosSimMtrx)}
        Result['Pcorrect'] = {'GESI': Metric2Pcorrect_Sigmoid(Result['d']['GESI'], GESIparam['Sigmoid'])}

    if GESIparam['SwPlot'] == 2:
        import matplotlib.pyplot as plt
        plt.imshow(Result['dIntrm']['GESI'] * 256, aspect='auto', origin='lower')
        plt.xlabel('MFB channel')
        plt.ylabel('GCFB channel')
        plt.show()

    return Result, GESIparam

# import os
# import numpy as np
# from scipy.signal import butter, filtfilt, resample
# from scipy.io import loadmat, savemat, wavfile
# from tool.TimeAlignXcorr import TimeAlignXcorr
# from tool.GCFBv234.tool.TaperWindow import TaperWindow
# from tool.FilterModFB import FilterModFB
# from tool.GCFBv234 import GCFBv234
# from tool.F0limit2SSIweight import F0limit2SSIweight
# from tool.Metric2Pcorrect_Sigmoid import Metric2Pcorrect_Sigmoid
# from tool.GCFBv234.tool.Eqlz2MeddisHCLevel import Eqlz2MeddisHCLevel
# from tool.GCFBv234.tool.EqlzGCFB2Rms1at0dB import EqlzGCFB2Rms1at0dB
# from tool.world_0_2_4_matlab import Harvest

# def GESIv123(SndRef, SndTest, GCparam, GESIparam):
#     # 目录信息
#     DirProg = os.path.dirname(__file__)
#     NameProg = os.path.splitext(os.path.basename(__file__))[0]
#     Result = {'Name': NameProg}
#     print(f'##### Start: {Result["Name"]} #####')

#     # 创建 GCout 目录
#     if 'DirGCout' not in GESIparam:
#         GESIparam['DirGCout'] = os.path.expanduser('~/Data/GESI/GCout/')
#         os.makedirs(GESIparam['DirGCout'], exist_ok=True)

#     GESIparam.setdefault('SwPlot', 0)

#     # 设定默认参数
#     GCparam.setdefault('fs', 48000)
#     GCparam.setdefault('NumCh', 100)
#     GCparam.setdefault('FRange', [100, 8000])
#     GCparam.setdefault('OutMidCrct', 'FreeField')
#     GCparam.setdefault('HLoss', {'Type': 'NH'})
#     GCparam['Ctrl'] = 'dynamic'
#     GCparam['DynHPAF'] = {'StrPrc': 'frame-base'}
#     GCparam['StrFloor'] = 'NoiseFloor'

#     GESIparam.setdefault('fs', GCparam['fs'])

#     # 处理音频
#     SndRef = np.asarray(SndRef, dtype=np.float32).flatten()
#     SndTest = np.asarray(SndTest, dtype=np.float32).flatten()

#     # 时间对齐
#     if len(SndTest) == len(SndRef) and 'SwTimeAlign' not in GESIparam:
#         print('GESI: No time alignment: length(SndRef) == length(SndTest)')
#         GESIparam['SwTimeAlign'] = 0
#     else:
#         GESIparam.setdefault('SwTimeAlign', 1)
#         if GESIparam['SwTimeAlign'] == 1:
#             SndTest, ParamTA = TimeAlignXcorr(SndTest, SndRef)
#             GESIparam['TimeAlign'] = ParamTA

#         GESIparam.setdefault('DurTaperWindow', 0.02)
#         LenTaper = int(GESIparam['DurTaperWindow'] * GESIparam['fs'])
        
#         Win = TaperWindow(len(SndRef), 'han', LenTaper)

#         # 如果 Win 是 tuple 或 list，只取第一个元素
#         if isinstance(Win, (tuple, list)):
#             Win = Win[0]

#         # 转换成 NumPy 数组
#         Win = np.asarray(Win, dtype=np.float32)
        
#         SndRef *= Win
#         SndTest *= Win

#     # GESI参数默认值
#     GESIparam.setdefault('Sim', {'PowerRatio': 0.6})
#     _, MFBparam = FilterModFB(Env,ParamMFB)
#     LenMFB = len(MFBparam['fc'])
#     GESIparam['Sim'].setdefault('weightMFB', np.ones(LenMFB))
#     GESIparam['Sim'].setdefault('SwWeightProhibit', 1)
#     GESIparam['Sim'].setdefault('RangeWeightProhibit', 1)

#     # 归一化 & 采样率调整
#     if GCparam['fs'] != GESIparam['fs']:
#         SndTest = resample(SndTest, GCparam['fs'], GESIparam['fs'])
#         SndRef = resample(SndRef, GCparam['fs'], GESIparam['fs'])

#     SndRef, MdsAmpdB = Eqlz2MeddisHCLevel(SndRef, [], GESIparam['DigitalRms1SPLdB'])
#     SndTest *= 10**(MdsAmpdB[1] / 20)

#     # 计算 GCFB
#     GCoutRef, _, GCparamRef = GCFBv234(SndRef, GCparam)
#     GCoutTest, _, GCparamTest = GCFBv234(SndTest, GCparam)
#     GCoutRef = EqlzGCFB2Rms1at0dB(GCoutRef, GCparam['StrFloor'])
#     GCoutTest = EqlzGCFB2Rms1at0dB(GCoutTest, GCparam['StrFloor'])
#     NumCh, LenFrame = GCoutRef.shape

#     # 计算相似度分析
#     CosSimMtrx = np.zeros((NumCh, LenMFB))
#     for nch in range(NumCh):
#         for nMFB in range(LenMFB):
#             CosSimMtrx[nch, nMFB] = np.nanmean(GESIparam['Sim']['weightMFB'])

#     Result['dIntrm'] = {'GESI': CosSimMtrx}
#     Result['d'] = {'GESI': np.nanmean(CosSimMtrx)}
#     Result['Pcorrect'] = {'GESI': Metric2Pcorrect_Sigmoid(Result['d']['GESI'], GESIparam['Sigmoid'])}

#     # 可视化
#     if GESIparam['SwPlot'] == 2:
#         import matplotlib.pyplot as plt
#         plt.imshow(Result['dIntrm']['GESI'] * 256, aspect='auto', origin='lower')
#         plt.xlabel('MFB channel')
#         plt.ylabel('GCFB channel')
#         plt.show()

#     return Result, GESIparam
