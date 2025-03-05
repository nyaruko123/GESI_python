# // %%  Test script for WORLD analysis/synthesis
# // % 2014/04/29: First version (0.1.4) for release.
# // % 2015/12/02: StoneMask.m was added to the project (0.2.0_5).
# // % 2016/12/28: Refactoring.
# // % 2017/01/02: Another example was added.
#
# // [x, fs] = audioread('vaiueo2d.wav');
#
# // if 0 % You can use Dio
# //   f0_parameter = Dio(x, fs);
# //   % StoneMask is an option for improving the F0 estimation performance.
# //   % You can skip this processing.
# //   f0_parameter.f0 = StoneMask(x, fs,...
# //     f0_parameter.temporal_positions, f0_parameter.f0);
# // end;
# // f0_parameter = Harvest(x, fs);
#
# // spectrum_parameter = CheapTrick(x, fs, f0_parameter);
# // source_parameter = D4C(x, fs, f0_parameter);
#
# // y = Synthesis(source_parameter, spectrum_parameter);
#
# // return;
#
# // %% Another example (we want to modify the parameters)
# // [x, fs] = audioread('vaiueo2d.wav');
# // option_harvest.f0_floor = 40;
# // f0_parameter = Harvest(x, fs, option_harvest);
#
# // % If you modified the fft_size, you must also modify the option in D4C.
# // % The lowest F0 that WORLD can work as expected is determined by the following:
# // % 3.0 * fs / fft_size
# // option_cheaptrick.fft_size = 4096;
# // option_d4c.fft_size = option_cheaptrick.fft_size;
# // spectrum_parameter = CheapTrick(x, fs, f0_parameter, option_cheaptrick);
# // source_parameter = D4C(x, fs, f0_parameter, option_d4c);
#
# // y = Synthesis(source_parameter, spectrum_parameter);

import soundfile as sf
import numpy as np
import Dio, StoneMask, Harvest, CheapTrick, D4C, Synthesis

# 读取音频文件
x, fs = sf.read('vaiueo2d.wav')

# 选择 F0 估计算法（使用 Dio 或 Harvest）
use_dio = False  # 设置为 True 可使用 Dio

if use_dio:
    f0_parameter = Dio(x, fs)  # 使用 Dio 进行 F0 估计
    
    # StoneMask 可选用于改进 F0 估计精度，可以跳过此步
    f0_parameter['f0'] = StoneMask(x, fs, f0_parameter['temporal_positions'], f0_parameter['f0'])
else:
    f0_parameter = Harvest(x, fs)  # 使用 Harvest 进行 F0 估计

# 计算频谱包络
spectrum_parameter = CheapTrick(x, fs, f0_parameter)

# 计算声源特性
source_parameter = D4C(x, fs, f0_parameter)

# 进行合成
y = Synthesis(source_parameter, spectrum_parameter)

# 另一个示例（调整参数）
option_harvest = {'f0_floor': 40}  # 修改 Harvest 参数
f0_parameter = Harvest(x, fs, option_harvest)

# 如果修改了 FFT 大小，则必须同步修改 D4C 选项
# WORLD 期望的最低 F0 由以下公式确定： 3.0 * fs / fft_size
option_cheaptrick = {'fft_size': 4096}
option_d4c = {'fft_size': option_cheaptrick['fft_size']}

spectrum_parameter = CheapTrick(x, fs, f0_parameter, option_cheaptrick)
source_parameter = D4C(x, fs, f0_parameter, option_d4c)

# 进行合成
y = Synthesis(source_parameter, spectrum_parameter)
