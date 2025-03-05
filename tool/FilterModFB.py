import numpy as np
from scipy.signal import butter, freqz, lfilter

def FilterModFB(Env, ParamMFB):
    # 通过调制滤波器组进行滤波
    # 默认中心频率
    ParamMFB.setdefault('fc_default', [1, 2, 4, 8, 16, 32, 64, 128, 256, 512])

    if Env is None:
        # 如果没有输入，返回默认设置
        OutMFB = []
        ParamMFB['fc'] = ParamMFB['fc_default']
        return OutMFB, ParamMFB

    if 'fs' not in ParamMFB:
        raise ValueError('Specify ParamMFB.fs')  # 检查采样频率
    if 'fc' not in ParamMFB:
        ParamMFB['fc'] = ParamMFB['fc_default']  # 使用默认中心频率
    if 'SwPlot' not in ParamMFB:
        ParamMFB['SwPlot'] = 0  # 默认不绘图

    # **使用 ParamMFB 作为存储滤波器系数的地方**
    if 'MFcoefIIR' not in ParamMFB:
        ParamMFB['MFcoefIIR'] = MkCoefModFilter(ParamMFB)  # 生成调制滤波器系数

    MFcoefIIR = ParamMFB['MFcoefIIR']

    LenFc = len(ParamMFB['fc'])
    NumEnv, LenEnv = Env.shape
    if NumEnv > 1:
        raise ValueError('Env should be a monoaural row vector.')  # 检查输入格式

    OutMFB = np.zeros((LenFc, LenEnv))
    for nfc in range(LenFc):
        OutMFB[nfc, :] = lfilter(MFcoefIIR['b'][nfc, :], MFcoefIIR['a'][nfc, :], Env)

    return OutMFB, ParamMFB  # 返回输出和参数


def MkCoefModFilter(ParamMFB):
    # 生成调制滤波器系数
    print('--- Making modulation filter coefficients ---')
    LenFc = len(ParamMFB['fc'])

    IIR_b = np.zeros((LenFc, 4))
    IIR_a = np.zeros((LenFc, 4))

    for nfc in range(LenFc):
        if ParamMFB['fc'][nfc] == 1:  # 当频率为 1 Hz
            # 三阶低通滤波器
            b, a = butter(3, ParamMFB['fc'][nfc] / (ParamMFB['fs'] / 2))
            b4 = b / a[0]
            a4 = a / a[0]

            IIR_b[nfc, :] = b4
            IIR_a[nfc, :] = a4
        else:  # 带通滤波器
            # 预扭曲
            w0 = 2 * np.pi * ParamMFB['fc'][nfc] / ParamMFB['fs']
            W0 = np.tan(w0 / 2)

            # 二阶带通滤波器
            Q = 1
            B0 = W0 / Q
            b = np.array([B0, 0, -B0])
            a = np.array([1 + B0 + W0**2, 2 * W0**2 - 2, 1 - B0 + W0**2])
            b3 = b / a[0]
            a3 = a / a[0]

            IIR_b[nfc, :3] = b3
            IIR_a[nfc, :3] = a3

    MFcoefIIR = {'a': IIR_a, 'b': IIR_b}  # 存储滤波器系数

    if ParamMFB['SwPlot'] == 1:  # 如果需要绘图
        PlotFrspMF(ParamMFB, MFcoefIIR)
        input('Return to continue > ')  # 等待用户确认

    return MFcoefIIR  # 返回滤波器系数


def PlotFrspMF(ParamMFB, MFcoefIIR):
    # 绘制数字滤波器的频率响应
    import matplotlib.pyplot as plt

    plt.figure()
    Nrsl = 1024 * 4  # 频率响应的点数

    for nfc in range(len(ParamMFB['fc'])):
        w, h = freqz(MFcoefIIR['b'][nfc, :], MFcoefIIR['a'][nfc, :], worN=Nrsl, fs=ParamMFB['fs'])
        plt.plot(w, 20 * np.log10(abs(h)))

    plt.box(True)
    plt.axis([0.25, max(ParamMFB['fc']) * 2, -40, 5])
    plt.grid()
    plt.xscale('log')
    plt.xticks(ParamMFB['fc'])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Filter attenuation (dB)')
    Str_FcMFB = [str(fc) for fc in ParamMFB['fc']]
    plt.legend(Str_FcMFB, loc='southwest')
    plt.title('Modulation filterbank')
    plt.show()  # 显示图形
