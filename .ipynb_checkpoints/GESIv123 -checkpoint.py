import os
import numpy as np
from scipy.signal import butter, filtfilt, resample
from scipy.io import loadmat, savemat, wavfile
from tool.TimeAlignXcorr import TimeAlignXcorr
from tool.GCFBv234.tool.TaperWindow import TaperWindow
from tool.FilterModFB import FilterModFB
from tool.GCFBv234 import GCFBv234
from tool.F0limit2SSIweight import F0limit2SSIweight
from tool.Metric2Pcorrect_Sigmoid import Metric2Pcorrect_Sigmoid
from tool.GCFBv234.tool.Eqlz2MeddisHCLevel import Eqlz2MeddisHCLevel
from tool.GCFBv234.tool.EqlzGCFB2Rms1at0dB import EqlzGCFB2Rms1at0dB
from tool.world_0_2_4_matlab import Harvest


def GESIv123(SndRef, SndTest, GCparam=None, GESIparam=None):
    """ GESI (GammaChirp Envelope Similarity Index) 计算函数 """
    
    # 确保 GCparam 和 GESIparam 不为空
    if GCparam is None:
        GCparam = {}
    if GESIparam is None:
        GESIparam = {}

    # 设定默认参数
    GCparam.setdefault('fs', 48000)
    GCparam.setdefault('NumCh', 100)
    GCparam.setdefault('FRange', [100, 8000])
    GCparam.setdefault('OutMidCrct', 'FreeField')
    GCparam.setdefault('HLoss', {'Type': 'NH'})
    GCparam['Ctrl'] = 'dynamic'
    GCparam['DynHPAF'] = {'StrPrc': 'frame-base'}
    GCparam['StrFloor'] = 'NoiseFloor'

    GESIparam.setdefault('fs', GCparam['fs'])
    GESIparam.setdefault('DigitalRms1SPLdB', 65)
    GESIparam.setdefault('SwPlot', 0)
    GESIparam.setdefault('Sim', {'PowerRatio': 0.6})

    # 读取音频文件（如果传入的是路径）
    if isinstance(SndRef, str):
        fs_ref, SndRef = wavfile.read(SndRef)
    if isinstance(SndTest, str):
        fs_test, SndTest = wavfile.read(SndTest)

    # 确保音频是 float32 类型，并且归一化
    SndRef = SndRef.astype(np.float32) / np.max(np.abs(SndRef))
    SndTest = SndTest.astype(np.float32) / np.max(np.abs(SndTest))

    # 确保音频是单声道
    if SndRef.ndim > 1:
        SndRef = np.mean(SndRef, axis=1)
    if SndTest.ndim > 1:
        SndTest = np.mean(SndTest, axis=1)

    # 进行时间对齐
    if 'SwTimeAlign' not in GESIparam:
        GESIparam['SwTimeAlign'] = 1
    if GESIparam['SwTimeAlign'] == 1:
        SndTest, ParamTA = TimeAlignXcorr(SndTest, SndRef)
        GESIparam['TimeAlign'] = ParamTA

    # 计算 GCFB
    print('==== 计算 GCFB（GammaChirp 滤波器组） ====')
    GCoutRef, _, GCparamRef = GCFBv234(SndRef, GCparam)
    GCoutTest, _, GCparamTest = GCFBv234(SndTest, GCparam)

    # 计算相似度
    CosSim = np.sum(GCoutRef * GCoutTest) / (np.sqrt(np.sum(GCoutRef**2)) * np.sqrt(np.sum(GCoutTest**2)))
    Result = {'GESI_Similarity': CosSim}

    # 结果输出
    print(f'计算完成，相似度: {Result["GESI_Similarity"]:.4f}')
    return Result


# 让脚本可以独立运行
if __name__ == "__main__":
    # 设定音频路径
    ref_audio_path = "audio_samples/ref_sound.wav"
    test_audio_path = "audio_samples/test_sound.wav"

    # 设定参数
    GCparam = {'fs': 48000, 'NumCh': 100, 'FRange': [100, 8000], 'OutMidCrct': 'FreeField', 'HLoss': {'Type': 'NH'}}
    GESIparam = {'fs': 48000, 'DigitalRms1SPLdB': 65, 'Sim': {'PowerRatio': 0.6}}

    # 运行 GESIv123
    result = GESIv123(ref_audio_path, test_audio_path, GCparam, GESIparam)
    print("最终结果:", result)
