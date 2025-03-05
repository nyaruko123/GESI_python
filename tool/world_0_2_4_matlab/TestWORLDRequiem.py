# // %%  Test script for WORLD analysis/synthesis with new waveform generator
# // % 2018/04/04: First version
#
# // [x, fs] = audioread('vaiueo2d.wav');
#
# // f0_parameter = Harvest(x, fs);
# // spectrum_parameter = CheapTrick(x, fs, f0_parameter);
# // source_parameter = D4CRequiem(x, fs, f0_parameter);
#
# // seeds_signals = GetSeedsSignals(fs);
# // y = SynthesisRequiem(source_parameter, spectrum_parameter, seeds_signals);

import librosa
import D4CRequiem, GetSeedsSignals, Harvest, CheapTrick, SynthesisRequiem

# 读取音频文件
x, fs = librosa.load('vaiueo2d.wav', sr=None)

# 进行 F0 提取
f0_parameter = Harvest(x, fs)

# 计算频谱参数
spectrum_parameter = CheapTrick(x, fs, f0_parameter)

# 计算声源参数
source_parameter = D4CRequiem(x, fs, f0_parameter)

# 生成激励信号种子
seeds_signals = GetSeedsSignals(fs)

# 进行合成
y = SynthesisRequiem(source_parameter, spectrum_parameter, seeds_signals)
