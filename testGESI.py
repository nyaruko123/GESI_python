# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#啊哈哈哈

# %
# %       Test code for GESI
# %       Irino, T.
# %       Created :  21 Mar 2022  IT from test_mrGEDIhl
# %       Modified:  21 Mar 2022  IT
# %       Modified:  24 Mar 2022  IT Modified Location
# %       Modified:  12 May 2022  IT Modified GESIparam.Sigmoid  --> See  GESI.m
# %       Modified:  26 May 2022  IT Modification on 12 May 2022 was cancelled.
# %       Modified:  29 Jul  2022  IT introduced: GESIparam.SwWeightProhibit
# %       Modified:   4 Aug 2022   IT  v110  introduction of version number +  normalization of SSI weight
# %       Modified:  22 Aug 2022   IT  v120  The order of input arguments was replaced
# %       Modified:  31 Aug 2022   IT  v121  Introduction of time-varying SSIweight
# %       Modified:  18 Oct  2022   IT  v122  adding rng()
# %       Modified:  19 Oct  2022   IT  v122  using GCFBv234
# %       Modified:  12 Nov 2022   IT  v123  version up. Tidy up. Renamed  from GESIv122_rPwrMulti.m (YA)
# %       Modified:  18 May 2023   IT  v123  adding some comments
# %
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# clear all
# close all
#
# %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Environment settings
# % Root
# DirProg = fileparts(which(mfilename)); % Directory of this m-file
# DirRoot = [DirProg '/'];
# %
# % Essential package: dynamic compressive gammachirp filterbank
# % (GCFBv233 or the later version)
# % Please download and put it at the same level of this directory.
# % https://github.com/AMLAB-Wakayama/gammachirp-filterbank/GCFBv234
# %
# %
# %DirGCFB = [DirRoot '../GCFBv233/'];  % normal install
# DirGCFB = [DirRoot '../../../GitHub_Public/gammachirp-filterbank/GCFBv234/'];  % local use only
# %exist(DirGCFB)   % for check directory
# addpath(DirGCFB)
# StartupGCFB;   % startup GCFB
#
# % Sounds
# DirSnd = [DirRoot 'wav_sample/'];
#
#
# %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% GEDI and materials
# % Parameters of dcGC filterbank
# GCparam.Ctrl       = 'dynamic';
# GCparam.NumCh  = 100;
# GCparam.fs          = 48000;
# GCparam.DynHPAF.StrPrc = 'frame';    % frame-base introduced by IT
# GCparam.HLoss.Type = 'NH';  % NH
# % GCparam.HLoss.Type = 'HL2';  % 80yr
# % GCparam.HLoss.CompressionHealth = 0.5;  % Compression health alpha
#
#
# CalibToneSPLdB = 65;
# CalibToneRMSDigitalLeveldB = -26;
# DigitalRms1SPLdB = CalibToneSPLdB - CalibToneRMSDigitalLeveldB;
#
# %% Parameter settings for GESI
# GESIparam.DigitalRms1SPLdB = DigitalRms1SPLdB;
# GESIparam.Sigmoid = [-20, 6]; % temporal value which should be modified --- 26 May 22
# GESIparam.Sim.PowerRatio = 0.6;  % power asymmetry valid for both NH + HI listeners
#         %GESIparam.Sim.PowerRatio = 0.5;  % do not use: valid only for NH listeners
#         %GESIparam.Sim.PowerRatio = 0.55;  % intermediate : adjustment for individual listeners
#         %GESIparam.Sim.PowerRatio = 0.57;  %
# %GESIparam.Sim.PowerRatio = 0.5;  % do not use: valid only for NH listeners
#
# % controlling weight matrix introduced 29 Jul 2022
# % GESIparam.Sim.SwWeightProhibit = 0; % No prohibit region : conventional weight
# % GESIparam.Sim.SwWeightProhibit = 1; % set prohibit region (default)  -- The result changes only slightly
#
# GESIparam.SwPlot = 2; %  image(Result.dIntrm.GCMFBtmc*256)
#
# % Parameter settings for materials
# SNRList = [-6, -3, 0, 3]; %SNR between clean speech and noise
#
# rng(12345);  % simulationでの再現性確保。GCFB出力でrandnを使っているので、同じにするため。
#
# %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Start simulation
#
# for nSnd = 1:length(SNRList)
#
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     %% Test signal (enhanced/Unprocessed speech); the speech intelligiblity is calculated
#     % Name of wav-file (example: '*/GEDI_Software/wav_sample/sample_sp1')
#     NameSndTest = ['sample_sp' num2str(nSnd)];
#     disp(['SndTest: ' NameSndTest]);
#     % Read wav-file of test speech
#     [SndTest, fs] = audioread([DirSnd NameSndTest '.wav']);
#
#
#     %% Reference signal (Clean speech)
#     % Name of wav-file
#     NameSndRef = 'sample_sp_clean';
#     disp(['SndRef : ' NameSndRef]);
#     % Read wav-file of clean speech
#     [SndRef, fs2] = audioread([DirSnd NameSndRef '.wav']);
#
#     if fs ~= fs2  %  IT
#         error('Inconsistency of sampling rate.');
#     end
#     GESIparam.fs = fs;   % Samping rate of sounds.   % NG:  GESIparam.fsSnd = fs;
#
#     % GCout mat file will be kept when the name is specified.
#     % These files will be used when GESI is executed again for fast processing.
#     % GESIparam.NameSndRef  = NameSndRef;
#     % GESIparam.NameSndTest = NameSndTest;
#
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     % Preparation of sound
#     %%%%%%%%
#     SwTimeAlign = 0; % Time alignment by apriori information
#     SwTimeAlign = 1; % test for TimeAlignXcorr in GESI
#     if SwTimeAlign == 0 % preparation here
#         disp('-- Extraction of a speech segment in SndTest from apriori information.')
#         TimeSndBefore   = 0.35;
#         % SndTest = SndTest(fs*TimeSndBefore+(0:length(SndRef)-1));  not very different
#         SndTest = SndTest(fs*TimeSndBefore+(1:length(SndRef)));
#
#         % Taper window
#         GESIparam.DurTaperWindow = 0.02; % 20ms taper window
#         LenTaper = GESIparam.DurTaperWindow *GESIparam.fs;
#         Win = TaperWindow(length(SndRef),'han',LenTaper);
#         SndRef   = SndRef(:).* Win(:);  % column vector
#         SndTest  = SndTest(:).* Win(:);
#     else
#         disp('-- Extraction of a speech segment in SndTest using TimeAlignXcorr in GESI.')
#     end
#
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     %% Speech intelligibility prediction by GESI
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     % [Result, GESIparam] = GESIv120(SndRef, SndTest, GCparam, GESIparam);  % v120: SndRef, SndTest
#     % [Result, GESIparam] = GESIv121(SndRef, SndTest, GCparam, GESIparam);  % v121
#     % [Result, GESIparam] = GESIv122(SndRef, SndTest, GCparam, GESIparam);  % v122
#     [Result, GESIparam] = GESIv123(SndRef, SndTest, GCparam, GESIparam);  % v123
#     Metric(nSnd)     = Result.d.GESI;
#     Pcorrects(nSnd) = Result.Pcorrect.GESI; % temporal value. It should be changed by the sigmoid parameters.
#
#     disp('==========================================');
#     disp(['Percent correct (temporally):' num2str(Pcorrects(nSnd)) '(%)']);
#     if GESIparam.SwTimeAlign> 0
#         disp(sprintf('TimeAlign : %d',GESIparam.TimeAlign.NumTimeLag))
#     end
#     disp('==========================================');
#
#     %%%%%%%%%%%%%%%%%%%
#     % plot figures
#     %%%%%%%%%%%%%%%%%%%
#     figure(nSnd)
#     subplot(1,2,1)
#     image(Result.dIntrm.GESI*256);
#     set(gca,'YDir','normal');
#     xlabel('MFB channel')
#     ylabel('GCFB channel')
#     title(['Metric: ' num2str(Metric(nSnd),'%5.3f') ...
#         ',  Pcorrect(tmp) : ' num2str(Pcorrects(nSnd),'%4.1f')])
#     drawnow
#
#     subplot(1,2,2)
#     image(GESIparam.SSIparam.weight*256*0.8);
#     set(gca,'YDir','normal');
#     xlabel('Frame')
#     ylabel('GCFB channel')
#     title(['SSIweight'])
#     drawnow
#
#
# end
#
# disp(['Pcorrect : ' num2str(Pcorrects)])
# disp(['Metric    : ' num2str(Metric)])
#
# disp('==========================================');
#
#
# %% Plot results
# figure(nSnd+1)
# plot(SNRList,Pcorrects,'o-');
# xlim([-0.5+min(SNRList) max(SNRList)+0.5]);
# ylim([-0.5+0 100+0.5]);
# xlabel('SNR (dB)');
# ylabel('Percent correct (%)')
# legend('Unprocessed')
# title('Results of GESI for example sounds');
# grid on;
#
#
# % Keep results for comparison   10 Jan 22
# if 0
#     NameRslt = 'Rslt_GESI';
#     save([NameRslt '_Val'])
#     print([NameRslt '_Fig'],'-depsc')
# end
# %
#
#
# %%%%%%
# %% Trash / memo
# %%%%%%%%
# % TaperWindowが、SndTestにかかるかどうかでも、%は異なってくる。
# %

import os
import numpy as np
import scipy.io.wavfile as wav
import matplotlib.pyplot as plt
from tool.GCFBv234.tool.TaperWindow import TaperWindow  # 仅引用，不实现
from tool.StartupGCFB import StartupGCFB  # 仅引用，不实现
from GESIv123 import GESIv123  # 仅引用，不实现

# 目录设置
DirProg = os.path.dirname(os.path.abspath(__file__))  # 获取当前文件所在目录
DirRoot = DirProg + '/'  # 根目录

# 必要的软件包：动态压缩 gammachirp 滤波器库
DirGCFB = '/home/sun/GESI_python/tool/GCFBv234'  # 本地使用
# 添加路径并启动 GCFB
os.sys.path.append(DirGCFB)
StartupGCFB()

# 声音文件目录
DirSnd = DirRoot + 'wav_sample/'

# GEDI 和材料参数
GCparam = {
    "Ctrl": "dynamic",
    "NumCh": 100,
    "fs": 48000,
    "DynHPAF": {"StrPrc": "frame"},
    "HLoss": {"Type": "NH"}  # 正常听觉
}

CalibToneSPLdB = 65
CalibToneRMSDigitalLeveldB = -26
DigitalRms1SPLdB = CalibToneSPLdB - CalibToneRMSDigitalLeveldB

GESIparam = {
    "DigitalRms1SPLdB": DigitalRms1SPLdB,
    "Sigmoid": [-20, 6],  # 暂定值
    "Sim": {"PowerRatio": 0.6},  # 功率不对称
    "SwPlot": 2  # 画图开关
}

# 材料参数
SNRList = [-6, -3, 0, 3]  # 语音与噪声的信噪比
np.random.seed(12345)  # 设置随机种子，保证模拟结果可复现

# 开始模拟
Metric = []
Pcorrects = []

for nSnd, snr in enumerate(SNRList):
    # 测试信号（增强/未处理语音）
    NameSndTest = f'sample_sp{nSnd}'
    print(f'SndTest: {NameSndTest}')
    fs, SndTest = wav.read(os.path.join(DirSnd, NameSndTest + '.wav'))

    # 参考信号（干净语音）
    NameSndRef = 'sample_sp_clean'
    print(f'SndRef : {NameSndRef}')
    fs2, SndRef = wav.read(os.path.join(DirSnd, NameSndRef + '.wav'))

    if fs != fs2:
        raise ValueError('采样率不一致')
    GESIparam["fs"] = fs  # 设定采样率

    # 语音对齐
    print('-- 通过 TimeAlignXcorr 进行语音对齐')

    # 计算语音可懂度
    Result, GESIparam = GESIv123(SndRef.astype(np.float32), SndTest.astype(np.float32), GCparam, GESIparam)
    Metric.append(Result["d"]["GESI"])
    Pcorrects.append(Result["Pcorrect"]["GESI"])

    print('==========================================')
    print(f'百分比正确率（临时）: {Pcorrects[nSnd]} (%)')
    if "TimeAlign" in GESIparam and "NumTimeLag" in GESIparam["TimeAlign"]:
        print(f'TimeAlign: {GESIparam["TimeAlign"]["NumTimeLag"]}')
    print('==========================================')

    # 画图
    plt.figure(nSnd)
    plt.subplot(1, 2, 1)
    plt.imshow(Result["dIntrm"]["GESI"] * 256, aspect='auto', origin='lower')
    plt.xlabel('MFB 通道')
    plt.ylabel('GCFB 通道')
    plt.title(f'Metric: {Metric[nSnd]:.3f}, Pcorrect: {Pcorrects[nSnd]:.1f}')
    plt.colorbar()

    plt.subplot(1, 2, 2)
    plt.imshow(GESIparam.get("SSIparam", {}).get("weight", np.zeros((10,10))) * 256 * 0.8, aspect='auto', origin='lower')
    plt.xlabel('帧')
    plt.ylabel('GCFB 通道')
    plt.title('SSI 权重')
    plt.colorbar()

    plt.show()

print(f'Pcorrect: {Pcorrects}')
print(f'Metric: {Metric}')
print('==========================================')

# 画最终结果图
plt.figure(len(SNRList) + 1)
plt.plot(SNRList, Pcorrects, 'o-')
plt.xlim([min(SNRList) - 0.5, max(SNRList) + 0.5])
plt.ylim([0 - 0.5, 100 + 0.5])
plt.xlabel('SNR (dB)')
plt.ylabel('百分比正确率 (%)')
plt.legend(['未处理'])
plt.title('GESI 示例声音结果')
plt.grid(True)
plt.show()


