import os
import numpy as np
from scipy.signal import butter, filtfilt, resample
from scipy.io import loadmat, savemat
from tool.TimeAlignXcorr import TimeAlignXcorr
from tool.gcfb_v234.utils import taper_window
from tool.FilterModFB import FilterModFB
from tool.gcfb_v234.gcfb_v234 import gcfb_v234
from tool.F0limit2SSIweight import F0limit2SSIweight
from tool.Metric2Pcorrect_Sigmoid import Metric2Pcorrect_Sigmoid
from tool.gcfb_v234.utils import eqlz2meddis_hc_level
from tool.gcfb_v234.utils import eqlz_gcfb2rms1_at_0db
# from tool.world_0_2_4_matlab.Harvest import Harvest
import pyworld as pw
import numpy as np



def resample_signal(signal, old_fs, new_fs):
    """ 计算目标采样点数并进行重采样 """
    num_samples = int(round(len(signal) * float(new_fs) / old_fs))  # 计算目标点数
    return resample(signal, num_samples)  # 进行重采样

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

    # 设置默认参数
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

    # 准备保存文件名
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

    # 声音维度检查
    if SndTest.ndim > 1 or SndRef.ndim > 1:
        raise ValueError('SndTest and SndRef should be monaural.')
    SndRef = SndRef.flatten()
    SndTest = SndTest.flatten()

    # ---------------------
    #  修改“自动对齐”部分
    # ---------------------
    # 若 SwTimeAlign 未设置，默认关闭对齐
    if 'SwTimeAlign' not in GESIparam:
        GESIparam['SwTimeAlign'] = 0

    if GESIparam['SwTimeAlign'] == 1:
        # 启用自动对齐
        print("TimeAlignXcorr is ENABLED. Doing cross-correlation alignment...")
        SndTest, ParamTA = TimeAlignXcorr(SndTest, SndRef)
        GESIparam['TimeAlign'] = ParamTA
    else:
        # 禁用自动对齐
        print("No time alignment performed. (SwTimeAlign=0)")
        ParamTA = {'NumTimeLag': 0, 'TimeLagInms': 0.0}
        GESIparam['TimeAlign'] = ParamTA

    # 如需加窗，则统一执行
    if 'DurTaperWindow' not in GESIparam:
        GESIparam['DurTaperWindow'] = 0.02
    LenTaper = int(GESIparam['DurTaperWindow'] * GESIparam['fs'])
    Win = taper_window(len(SndRef), 'han', LenTaper)

    # 如果 Win 是 tuple 或 list，只取第一个元素
    if isinstance(Win, (tuple, list)):
        Win = Win[0]
    Win = np.asarray(Win, dtype=np.float32)

    SndRef *= Win
    SndTest *= Win
    # ---------------------

    # GESI params
    if 'Sim' not in GESIparam or 'PowerRatio' not in GESIparam['Sim']:
        GESIparam['Sim'] = {'PowerRatio': 0.6}
        print('GESIparam.Sim.PowerRatio is set to 0.6 (default) -- OK? Return to continue > ')
        input()

    # 获取调制滤波器组参数
    _, MFBparam = FilterModFB(None, None)
    if 'weightMFB' not in GESIparam['Sim']:
        LenMFB = len(MFBparam['fc'])
        GESIparam['Sim']['weightMFB'] = np.ones(LenMFB)
    else:
        LenMFB = len(GESIparam['Sim']['weightMFB'])

    if 'SwWeightProhibit' not in GESIparam['Sim']:
        GESIparam['Sim']['SwWeightProhibit'] = 1
    if 'RangeWeightProhibit' not in GESIparam['Sim']:
        GESIparam['Sim']['RangeWeightProhibit'] = 1

    # 若 GESIparam.fs 与 GCparam.fs 不一致，则做重采样
    if GCparam['fs'] != GESIparam['fs']:
        SndTest = resample_signal(SndTest, GESIparam['fs'], GCparam['fs'])
        SndRef = resample_signal(SndRef, GESIparam['fs'], GCparam['fs'])

    # 进行听觉级别校准
    SndRef, MdsAmpdB = eqlz2meddis_hc_level(SndRef, None, GESIparam['DigitalRms1SPLdB'])
    print("SndRef (calibrated): mean =", np.mean(SndRef), ", std =", np.std(SndRef))#debug
    # 同样缩放 SndTest
    SndTest *= 10**(MdsAmpdB[1] / 20)
    print("SndTest (scaled): mean =", np.mean(SndTest), ", std =", np.std(SndTest))#debug

    # 去除 DC 分量
    SndRef = SndRef - np.mean(SndRef)
    SndTest = SndTest - np.mean(SndTest)
    print("After DC removal:")#debug
    print("SndRef: mean =", np.mean(SndRef), ", std =", np.std(SndRef))
    print("SndTest: mean =", np.mean(SndTest), ", std =", np.std(SndTest))


    # --------------------
    # 计算 GCFB 输出 (Test)
    # --------------------
    DirNameTest = os.path.join(GESIparam['DirGCout'], f'{GESIparam["NameGCoutTest"]}.mat')
    if not os.path.exists(DirNameTest):
        print('==== GCFB calculation of SndTest (HL or NH) ====')
        GCoutTest, _, GCparamTest, GCrespTest = gcfb_v234(SndTest, GCparam)
        NumCh, LenFrame = GCoutTest.shape
        GCoutTest = eqlz_gcfb2rms1_at_0db(GCoutTest, GCparam['StrFloor'])

        #print(type(GCparamTest), type(GCparamTest.get('DynHPAF'))) #debug
        #print(type(GCparamTest), type(GCparamTest.DynHPAF)) #debug


        #MFBparam['fs'] = GCparamTest['DynHPAF']['fs']

        #print(">>> DEBUG: DynHPAF keys:", GCparamTest.DynHPAF.keys()) #debug

        # print(">>> DEBUG: GCparamTest 内容:")
        # for key in vars(GCparamTest):
        #     print(f"{key}: {getattr(GCparamTest, key)}")   #debug


        #MFBparam['fs'] = GCparamTest.DynHPAF['fs']
        MFBparam['fs'] = GCparamTest.dyn_hpaf.fs


       

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
            savemat(DirNameTest, {
                'GCoutTest': GCoutTest,
                'GCparamTest': GCparamTest,
                'GCModEnvTest': GCModEnvTest,
                'MFBparam': MFBparam
            })
    else:
        print(f'load: {DirNameTest}')
        data = loadmat(DirNameTest)
        GCoutTest = data['GCoutTest']
        GCparamTest = data['GCparamTest']
        GCModEnvTest = data['GCModEnvTest']
        MFBparam = data['MFBparam']

    # --------------------
    # 计算 GCFB 输出 (Ref)
    # --------------------
    DirNameRef = os.path.join(GESIparam['DirGCout'], f'{GESIparam["NameGCoutRef"]}.mat')
    if not os.path.exists(DirNameRef):
        print('==== GCFB calculation of SndRef (always NH) ====')
        GCparam['HLoss']['Type'] = 'NH'
        GCoutRef, _, GCparamRef, GCrespRef = gcfb_v234(SndRef, GCparam)
        NumCh, LenFrame = GCoutRef.shape
        GCoutRef = eqlz_gcfb2rms1_at_0db(GCoutRef, GCparam['StrFloor'])
        #tFrame = np.arange(LenFrame) / GCparamRef['DynHPAF']['fs']
        tFrame = np.arange(LenFrame) / GCparamRef.dyn_hpaf.fs


        GCModEnvRef = np.zeros((NumCh, LenMFB, LenFrame))
        for nch in range(GCparam['NumCh']):
            Env = filtfilt(MFBparam['bzLPF'], MFBparam['apLPF'], GCoutRef[nch, :])
            if nch in [0, 49, 99]:
                print(f'> Modulation Filterbank Analysis: GCFB ch = #{nch + 1}')
            ModEnv, _ = FilterModFB(Env, MFBparam)
            GCModEnvRef[nch, :, :] = ModEnv

        # HarvestRef = Harvest(SndRef, GCparam['fs'])
        # # 计算帧对应的 F0
        # F0Frame = np.interp(tFrame, HarvestRef['temporal_positions'], HarvestRef['f0'], left=0, right=0)
        # F0Frame = np.maximum(F0Frame, 0)

        # # 计算 Ref 的平均 F0（对数平均）
        # valid_f0 = HarvestRef['f0'][HarvestRef['f0'] > 0]
        # if len(valid_f0) > 0:
        #     F0MeanRef = np.exp(np.mean(np.log(valid_f0)))
        # else:
        #     F0MeanRef = np.nan
        # print(f'Fo Mean of Ref sound: {F0MeanRef} Hz')

        

        # 使用 pyworld 提取 F0（Harvest + Stonemask 精修）
        _f0, _time_axis = pw.harvest(SndRef.astype(np.float64), GCparam['fs'])  # 初步估计
        f0_refined = pw.stonemask(SndRef.astype(np.float64), _f0, _time_axis, GCparam['fs'])  # 精细调整

        # 构造 HarvestRef，保持和原来一致的格式
        HarvestRef = {
        'f0': f0_refined,
        'temporal_positions': _time_axis
        }

        # 计算帧对应的 F0
        F0Frame = np.interp(tFrame, HarvestRef['temporal_positions'], HarvestRef['f0'], left=0, right=0)
        F0Frame = np.maximum(F0Frame, 0)

        # 计算 Ref 的平均 F0（对数平均）
        valid_f0 = HarvestRef['f0'][HarvestRef['f0'] > 0]
        if len(valid_f0) > 0:
            F0MeanRef = np.exp(np.mean(np.log(valid_f0)))
        else:
            F0MeanRef = np.nan
        print(f'Fo Mean of Ref sound: {F0MeanRef} Hz')



        if np.any(np.isnan(F0Frame)):
            raise ValueError('Error in F0Frame: F0 contains NaN.')

        #SSIparam = {'SwSSIweight': 2, 'Fr1': GCparamRef['Fr1']}
        SSIparam = {'SwSSIweight': 2, 'Fr1': GCparamRef.fr1}

        SSIweightMtrx = np.zeros((NumCh, LenFrame))

        for nFrame in range(LenFrame):
            SSIparam['F0_limit'] = F0Frame[nFrame]
            SSIweight, SSIparam = F0limit2SSIweight(SSIparam)
            # 归一化
            
            #SSIweightMtrx[:, nFrame] = SSIweight / np.mean(SSIweight)
            SSIweightMtrx[:, nFrame] = (SSIweight / np.mean(SSIweight)).flatten()


        SSIparam['weight'] = SSIweightMtrx

        if SwSave:
            print(f'save: {DirNameRef}')
            savemat(DirNameRef, {
                'GCoutRef': GCoutRef,
                'GCparamRef': GCparamRef,
                'GCModEnvRef': GCModEnvRef,
                'SSIparam': SSIparam
            })
    else:
        print(f'load: {DirNameRef}')
        data = loadmat(DirNameRef)
        GCoutRef = data['GCoutRef']
        GCparamRef = data['GCparamRef']
        GCModEnvRef = data['GCModEnvRef']
        SSIparam = data['SSIparam']

    # 将结果记录到 GESIparam
    GESIparam['GCparam'] = {'Ref': GCparamRef, 'Test': GCparamTest}
    GESIparam['MFBparam'] = MFBparam
    GESIparam['SSIparam'] = SSIparam

    # 若 SwPlot == 1, 画一些中间结果
    if GESIparam['SwPlot'] == 1:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.imshow(6 * GCoutTest, aspect='auto', origin='lower')
        plt.title("GCoutTest")
        plt.subplot(2, 1, 2)
        plt.imshow(6 * GCoutRef, aspect='auto', origin='lower')
        plt.title("GCoutRef")
        plt.show()

        plt.figure()
        plt.plot(SSIparam['weight'])
        plt.title("SSIparam['weight']")
        plt.show()

    # --------------------
    # 计算 Cosine 相似度
    # --------------------
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

                # 余弦相似度带权
                numerator = np.sum(weightGCFB * ModEnvRef * ModEnvTest)
                denominator = (PwrRef**rPwr) * (PwrTest**(1 - rPwr))
                if denominator == 0:
                    CosSim = 0.0
                else:
                    CosSim = numerator / denominator

                # 如果频率过低或过高要屏蔽
                if GESIparam['Sim']['SwWeightProhibit'] == 1:
                    #if GCparamRef['Fr1'][nch] < MFBparam['fc'][nMFB] * GESIparam['Sim']['RangeWeightProhibit']:
                    # if GCparamRef.Fr1[nch] < MFBparam['fc'][nMFB] * GESIparam['Sim']['RangeWeightProhibit']:
                    if GCparamRef.fr1[nch] < MFBparam['fc'][nMFB] * GESIparam['Sim']['RangeWeightProhibit']:


                        weightMFB[nMFB] = np.nan

                CosSimMtrx[nch, nMFB] = weightMFB[nMFB] * CosSim

        # 中间结果
        Result['dIntrm'] = {'GESI': CosSimMtrx}
        # 主指标
        Result['d'] = {'GESI': np.nanmean(CosSimMtrx)}
        # 百分比正确率（Sigmoid）
        Result['Pcorrect'] = {'GESI': Metric2Pcorrect_Sigmoid(Result['d']['GESI'], GESIparam['Sigmoid'])}

    # 若 SwPlot == 2, 显示最终 CosSim 矩阵
    if GESIparam['SwPlot'] == 2:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.imshow(Result['dIntrm']['GESI'] * 256, aspect='auto', origin='lower')
        plt.xlabel('MFB channel')
        plt.ylabel('GCFB channel')
        plt.title("CosSim Matrix (dIntrm['GESI'])")
        plt.colorbar()
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
