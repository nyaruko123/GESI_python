# // function WorldSynthesizer(input_filename, output_filename, f0_param, spec_param, time_param)
# // % F0���V�t�g�C�X�y�N�g���̐L�k�C���b���Ԃ̐L�k�����{����֐�
# // % �����FWorldSynthesizer(input_filename, output_filename, f0_param, spec_param, time_param)
# // % ���FWorldSynthesizer('vaiueo2d.wav', 'output.wav', 1, 1, 1);
# // [x, fs] = audioread(input_filename);
#
# // f0_parameter = Harvest(x, fs);
#
# // spectrum_parameter = CheapTrick(x, fs, f0_parameter);
# // source_parameter = D4C(x, fs, f0_parameter);
#
# // % F0�̕ύX
# // source_parameter.f0 = source_parameter.f0 * f0_param;
#
# // % �X�y�N�g���̐L�k
# // fft_size = (size(spectrum_parameter.spectrogram, 1) - 1) * 2;
# // w = (0 : fft_size - 1) * fs / fft_size;
# // w2 = (0 : fft_size / 2) * fs / fft_size / spec_param;
# // for i = 1 : size(spectrum_parameter.spectrogram, 2)
# //   tmp = [spectrum_parameter.spectrogram(:, i); spectrum_parameter.spectrogram(end - 1 : -1 : 2, i)];
# //   spectrum_parameter.spectrogram(:, i) = interp1(w, tmp, w2, 'linear', 'extrap');
# // end;
#
# // % ���b���Ԃ̐L�k
# // source_parameter.temporal_positions = source_parameter.temporal_positions * time_param;
#
# // y = Synthesis(source_parameter, spectrum_parameter);
#
# // audiowrite(output_filename, y, fs);

import numpy as np
import librosa
import scipy.interpolate
import soundfile as sf
import Harvest,CheapTrick,D4C,Synthesis

def WorldSynthesizer(input_filename, output_filename, f0_param, spec_param, time_param):
    """
    F0 变换、频谱变换、时间变换的音频合成函数
    调用方式：WorldSynthesizer('vaiueo2d.wav', 'output.wav', 1, 1, 1)
    """
    # 读取音频文件
    x, fs = librosa.load(input_filename, sr=None)
    
    # 计算 F0 参数
    f0_parameter = Harvest(x, fs)
    
    # 计算频谱参数
    spectrum_parameter = CheapTrick(x, fs, f0_parameter)
    
    # 计算声源参数
    source_parameter = D4C(x, fs, f0_parameter)
    
    # F0 的变换
    source_parameter['f0'] *= f0_param
    
    # 频谱的变换
    fft_size = (spectrum_parameter['spectrogram'].shape[0] - 1) * 2
    w = np.arange(fft_size) * fs / fft_size
    w2 = np.arange(fft_size // 2 + 1) * fs / fft_size / spec_param
    
    for i in range(spectrum_parameter['spectrogram'].shape[1]):
        tmp = np.concatenate([
            spectrum_parameter['spectrogram'][:, i],
            spectrum_parameter['spectrogram'][-2:0:-1, i]
        ])
        interp_func = scipy.interpolate.interp1d(w, tmp, kind='linear', fill_value='extrapolate')
        spectrum_parameter['spectrogram'][:, i] = interp_func(w2)
    
    # 时间轴的变换
    source_parameter['temporal_positions'] *= time_param
    
    # 进行合成
    y = Synthesis(source_parameter, spectrum_parameter)
    
    # 保存合成的音频
    sf.write(output_filename, y, fs)

# 这里假设 Harvest, CheapTrick, D4C, Synthesis 是已定义的函数
# 如果需要实现它们，可以参考 WORLD 的 Python 版本，或者提供更多细节
