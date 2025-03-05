# // function y = Synthesis(source_object, filter_object)
# // % Waveform synthesis from the estimated parameters
# // % y = Synthesis(source_object, filter_object)
# // %
# // % Input
# // %   source_object : F0 and aperiodicity
# // %   filter_object : spectral envelope
# // %
# // % Output
# // %   y : synthesized waveform
# // %
# // % 2016/12/28: Refactoring
# // % 2017/05/20: DC removal was fixed.
# // % 2017/06/12: Fix phase continuation
#
# // % rng(1);  % comment out by IT 24 Jul 19
# // default_f0 = 500;
# // spectrogram = filter_object.spectrogram;
# // temporal_positions = source_object.temporal_positions;
# // time_axis =...
# //   temporal_positions(1) : 1 / filter_object.fs : temporal_positions(end);
# // y = 0 * time_axis';
#
# // [pulse_locations, pulse_locations_index, pulse_locations_time_shift,...
# //   interpolated_vuv] = ...
# //   TimeBaseGeneration(temporal_positions, source_object.f0, filter_object.fs,...
# //   source_object.vuv, time_axis, default_f0);
#
# // fft_size = (size(spectrogram, 1) - 1) * 2;
# // base_index = -fft_size / 2 + 1 : fft_size / 2;
# // latter_index = fft_size / 2 + 1 : fft_size;
#
# // temporal_position_index = interp1(temporal_positions, ...
# //   1 : length(temporal_positions), pulse_locations, 'linear', 'extrap');
# // temporal_position_index = max(1, min(length(temporal_positions),...
# //   temporal_position_index));
#
# // amplitude_aperiodic = source_object.aperiodicity .^ 2;
# // amplitude_periodic = max(0.001, (1 - amplitude_aperiodic));
#
# // dc_remover_base = hanning(fft_size);
# // dc_remover_base = dc_remover_base / sum(dc_remover_base);
# // coefficient = 2.0 * pi * filter_object.fs / fft_size;
#
# // for i = 1 : length(pulse_locations_index)
# //   [spectrum_slice, periodic_slice, aperiodic_slice] = ...
# //     GetSpectralParameters(temporal_positions, temporal_position_index(i),...
# //     spectrogram, amplitude_periodic, amplitude_aperiodic, pulse_locations(i));
#
# //   noise_size = ...
# //     pulse_locations_index(min(length(pulse_locations_index), i + 1)) -...
# //     pulse_locations_index(i);
# //   output_buffer_index = ...
# //     max(1, min(length(y), pulse_locations_index(i) + base_index));
#
# //   if interpolated_vuv(pulse_locations_index(i)) > 0.5 &&...
# //       aperiodic_slice(1) <= 0.999
# //     response = GetPeriodicResponse(spectrum_slice, periodic_slice,...
# //       fft_size, latter_index, pulse_locations_time_shift(i), coefficient);
# //     dc_remover = dc_remover_base * -sum(response);
# //     periodic_response = response + dc_remover;
# //     y(output_buffer_index) =...
# //       y(output_buffer_index) + periodic_response * sqrt(max(1, noise_size));
# //     tmp_aperiodic_spectrum = spectrum_slice .* aperiodic_slice;
# //   else
# //     tmp_aperiodic_spectrum = spectrum_slice;
# //   end;
#
# //   tmp_aperiodic_spectrum(tmp_aperiodic_spectrum == 0) = eps;
# //   aperiodic_response = GetAperiodicResponse(tmp_aperiodic_spectrum,...
# //     fft_size, latter_index, noise_size);
# //   y(output_buffer_index) = y(output_buffer_index) + aperiodic_response;
# // end;
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function response = GetAperiodicResponse(tmp_aperiodic_spectrum,...
# //   fft_size, latter_index, noise_size)
# // aperiodic_spectrum =...
# //   [tmp_aperiodic_spectrum; tmp_aperiodic_spectrum(end - 1 : -1 : 2)];
# // tmp_cepstrum = real(fft(log(abs(aperiodic_spectrum)') / 2));
# // tmp_complex_cepstrum = zeros(fft_size, 1);
# // tmp_complex_cepstrum(latter_index) = tmp_cepstrum(latter_index) * 2;
# // tmp_complex_cepstrum(1) = tmp_cepstrum(1);
# // response = fftshift(real(ifft(exp(ifft(tmp_complex_cepstrum)))));
# // noise_input = randn(max(3, noise_size), 1);
# // response = fftfilt(noise_input - mean(noise_input), response);
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function response = GetPeriodicResponse(spectrum_slice, periodic_slice,...
# //   fft_size, latter_index, fractionsl_time_shift, coefficient)
# // tmp_periodic_spectrum = spectrum_slice .* periodic_slice;
# // tmp_periodic_spectrum(tmp_periodic_spectrum == 0) = eps;
# // periodic_spectrum =...
# //   [tmp_periodic_spectrum; tmp_periodic_spectrum(end - 1 : -1 : 2)];
#
# // tmp_cepstrum = real(fft(log(abs(periodic_spectrum)) / 2));
# // tmp_complex_cepstrum = zeros(fft_size, 1);
# // tmp_complex_cepstrum(latter_index) = tmp_cepstrum(latter_index) * 2;
# // tmp_complex_cepstrum(1) = tmp_cepstrum(1);
#
# // spectrum = exp(ifft(tmp_complex_cepstrum));
# // spectrum = spectrum(1 : fft_size / 2 + 1);
# // spectrum = spectrum .*...
# //   exp(-1i * coefficient * fractionsl_time_shift * (0 : fft_size / 2)');
# // spectrum = [spectrum; conj(spectrum(end - 1 : -1 : 2))];
# // response = fftshift(real(ifft(spectrum)));
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function [pulse_locations, pulse_locations_index,...
# //   pulse_locations_time_shift,vuv_interpolated] =...
# //   TimeBaseGeneration(temporal_positions, f0, fs, vuv, time_axis, default_f0)
# // f0_interpolated_raw = ...
# //   interp1(temporal_positions, f0, time_axis, 'linear', 'extrap');
# // vuv_interpolated = ...
# //   interp1(temporal_positions, vuv, time_axis, 'linear', 'extrap');
# // vuv_interpolated = vuv_interpolated > 0.5;
#
# // f0_interpolated = f0_interpolated_raw .* vuv_interpolated;
# // f0_interpolated(f0_interpolated == 0) = ...
# //   f0_interpolated(f0_interpolated == 0) + default_f0;
#
# // total_phase = cumsum(2 * pi * f0_interpolated / fs);
# // wrap_phase = rem(total_phase, 2 * pi);
# // pulse_locations = time_axis(abs(diff(wrap_phase)) > pi);
# // pulse_locations_index = round(pulse_locations * fs) + 1;
#
# // y1 = wrap_phase(pulse_locations_index) - 2.0 * pi;
# // y2 = wrap_phase(pulse_locations_index + 1);
# // x = -y1 ./ (y2 - y1);
# // pulse_locations_time_shift = x / fs;
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function [spectrum_slice, periodic_slice, aperiodic_slice] = ...
# //   GetSpectralParameters(temporal_positions, temporal_position_index,...
# //   spectrogram, amplitude_periodic, amplitude_random, pulse_locations)
# // floor_index = floor(temporal_position_index);
# // ceil_index = ceil(temporal_position_index);
# // t1 = temporal_positions(floor_index);
# // t2 = temporal_positions(ceil_index);
#
# // if t1 == t2
# //   spectrum_slice = spectrogram(:, floor_index);
# //   periodic_slice = amplitude_periodic(:, floor_index);
# //   aperiodic_slice = amplitude_random(:, floor_index);
# // else
# //   spectrum_slice = ...
# //     interp1q([t1 t2], [spectrogram(:, floor_index) ...
# //     spectrogram(:, ceil_index)]', max(t1, min(t2, pulse_locations)))';
# //   periodic_slice = ...
# //     interp1q([t1 t2], [amplitude_periodic(:, floor_index) ...
# //     amplitude_periodic(:, ceil_index)]', max(t1, min(t2, pulse_locations)))';
# //   aperiodic_slice = ...
# //     interp1q([t1 t2], [amplitude_random(:, floor_index) ...
# //     amplitude_random(:, ceil_index)]', max(t1, min(t2, pulse_locations)))';
# // end;

import numpy as np
from scipy.fftpack import fft, ifft, fftshift
from scipy.signal import fftconvolve
from scipy.signal.windows import hann


def Synthesis(source_object, filter_object):
    """
    从估计参数合成波形
    """
    default_f0 = 500
    spectrogram = filter_object.spectrogram
    temporal_positions = source_object.temporal_positions
    time_axis = np.arange(temporal_positions[0], temporal_positions[-1], 1 / filter_object.fs)
    y = np.zeros_like(time_axis)

    pulse_locations, pulse_locations_index, pulse_locations_time_shift, interpolated_vuv = \
        TimeBaseGeneration(temporal_positions, source_object.f0, filter_object.fs,
                           source_object.vuv, time_axis, default_f0)

    fft_size = 2 ** int(np.ceil(np.log2((spectrogram.shape[0] - 1) * 2)))
    base_index = np.arange(-fft_size // 2 + 1, fft_size // 2 + 1)
    latter_index = np.arange(fft_size // 2 + 1, fft_size)

    temporal_position_index = np.interp(pulse_locations, temporal_positions, np.arange(len(temporal_positions)))
    temporal_position_index = np.clip(temporal_position_index, 1, len(temporal_positions) - 1)

    amplitude_aperiodic = source_object.aperiodicity ** 2
    amplitude_periodic = np.maximum(0.001, (1 - amplitude_aperiodic))

    dc_remover_base = hann(fft_size)
    dc_remover_base /= np.sum(dc_remover_base)
    coefficient = 2.0 * np.pi * filter_object.fs / fft_size

    for i in range(len(pulse_locations_index)):
        spectrum_slice, periodic_slice, aperiodic_slice = \
            GetSpectralParameters(temporal_positions, temporal_position_index[i],
                                  spectrogram, amplitude_periodic, amplitude_aperiodic, pulse_locations[i])

        noise_size = pulse_locations_index[min(len(pulse_locations_index) - 1, i + 1)] - pulse_locations_index[i]
        output_buffer_index = np.clip(pulse_locations_index[i] + base_index, 0, len(y) - 1)

        if interpolated_vuv[pulse_locations_index[i]] > 0.5 and aperiodic_slice[0] <= 0.999:
            response = GetPeriodicResponse(spectrum_slice, periodic_slice,
                                           fft_size, latter_index, pulse_locations_time_shift[i], coefficient)
            dc_remover = dc_remover_base * -np.sum(response)
            periodic_response = response + dc_remover
            y[output_buffer_index] += periodic_response * np.sqrt(max(1, noise_size))
            tmp_aperiodic_spectrum = spectrum_slice * aperiodic_slice
        else:
            tmp_aperiodic_spectrum = spectrum_slice

        tmp_aperiodic_spectrum[tmp_aperiodic_spectrum == 0] = np.finfo(float).eps
        aperiodic_response = GetAperiodicResponse(tmp_aperiodic_spectrum, fft_size, latter_index, noise_size)
        y[output_buffer_index] += aperiodic_response

    return y

def GetAperiodicResponse(tmp_aperiodic_spectrum, fft_size, latter_index, noise_size):
    """
    获取非周期性响应
    """
    aperiodic_spectrum = np.concatenate((tmp_aperiodic_spectrum, tmp_aperiodic_spectrum[-2:0:-1]))
    tmp_cepstrum = np.real(fft(np.log(np.abs(aperiodic_spectrum)) / 2))
    tmp_complex_cepstrum = np.zeros(fft_size)
    tmp_complex_cepstrum[latter_index] = tmp_cepstrum[latter_index] * 2
    tmp_complex_cepstrum[0] = tmp_cepstrum[0]
    response = fftshift(np.real(ifft(np.exp(ifft(tmp_complex_cepstrum)))))
    noise_input = np.random.randn(max(3, noise_size))
    response = fftconvolve(noise_input - np.mean(noise_input), response, mode='same')
    return response


def GetPeriodicResponse(spectrum_slice, periodic_slice, fft_size, latter_index, fractional_time_shift, coefficient):
    """
    获取周期性响应
    """
    tmp_periodic_spectrum = spectrum_slice * periodic_slice
    tmp_periodic_spectrum[tmp_periodic_spectrum == 0] = np.finfo(float).eps
    periodic_spectrum = np.concatenate((tmp_periodic_spectrum, tmp_periodic_spectrum[-2:0:-1]))

    tmp_cepstrum = np.real(fft(np.log(np.abs(periodic_spectrum)) / 2))
    tmp_complex_cepstrum = np.zeros(fft_size)
    tmp_complex_cepstrum[latter_index] = tmp_cepstrum[latter_index] * 2
    tmp_complex_cepstrum[0] = tmp_cepstrum[0]

    spectrum = np.exp(ifft(tmp_complex_cepstrum))
    spectrum = spectrum[:fft_size // 2 + 1]
    spectrum *= np.exp(-1j * coefficient * fractional_time_shift * np.arange(fft_size // 2 + 1))
    spectrum = np.concatenate((spectrum, np.conj(spectrum[-2:0:-1])))

    response = fftshift(np.real(ifft(spectrum)))
    return response


def TimeBaseGeneration(temporal_positions, f0, fs, vuv, time_axis, default_f0):
    """
    生成时间基准
    """
    f0_interpolated_raw = np.interp(time_axis, temporal_positions, f0, left=default_f0, right=default_f0)
    vuv_interpolated = np.interp(time_axis, temporal_positions, vuv) > 0.5
    f0_interpolated = f0_interpolated_raw * vuv_interpolated
    f0_interpolated[f0_interpolated == 0] = default_f0

    total_phase = np.cumsum(2 * np.pi * f0_interpolated / fs)
    wrap_phase = np.remainder(total_phase, 2 * np.pi)
    pulse_locations = time_axis[np.abs(np.diff(wrap_phase)) > np.pi]
    pulse_locations_index = np.round(pulse_locations * fs).astype(int)
    pulse_locations_index = np.clip(pulse_locations_index, 0, len(time_axis) - 2)

    y1 = wrap_phase[pulse_locations_index] - 2.0 * np.pi
    y2 = wrap_phase[pulse_locations_index + 1]
    x = -y1 / (y2 - y1)
    pulse_locations_time_shift = x / fs

    return pulse_locations, pulse_locations_index, pulse_locations_time_shift, vuv_interpolated


import numpy as np


def GetSpectralParameters(temporal_positions, temporal_position_index, spectrogram,
                          amplitude_periodic, amplitude_aperiodic, pulse_location):
    """
    获取频谱参数（spectrum_slice, periodic_slice, aperiodic_slice）。

    :param temporal_positions: 时间位置数组
    :param temporal_position_index: 当前脉冲的时间索引
    :param spectrogram: 频谱数据
    :param amplitude_periodic: 周期性振幅
    :param amplitude_aperiodic: 非周期性振幅
    :param pulse_location: 当前脉冲的位置
    :return: spectrum_slice, periodic_slice, aperiodic_slice
    """

    # 计算 floor 和 ceil 索引
    floor_index = int(np.floor(temporal_position_index))
    ceil_index = int(np.ceil(temporal_position_index))

    # 处理索引越界情况
    floor_index = max(0, min(len(temporal_positions) - 1, floor_index))
    ceil_index = max(0, min(len(temporal_positions) - 1, ceil_index))

    t1 = temporal_positions[floor_index]
    t2 = temporal_positions[ceil_index]

    if t1 == t2:
        spectrum_slice = spectrogram[:, floor_index]
        periodic_slice = amplitude_periodic[:, floor_index]
        aperiodic_slice = amplitude_aperiodic[:, floor_index]
    else:
        # 进行线性插值
        spectrum_slice = np.interp(pulse_location, [t1, t2],
                                   np.vstack([spectrogram[:, floor_index], spectrogram[:, ceil_index]]).T)
        periodic_slice = np.interp(pulse_location, [t1, t2],
                                   np.vstack([amplitude_periodic[:, floor_index], amplitude_periodic[:, ceil_index]]).T)
        aperiodic_slice = np.interp(pulse_location, [t1, t2],
                                    np.vstack(
                                        [amplitude_aperiodic[:, floor_index], amplitude_aperiodic[:, ceil_index]]).T)

    return spectrum_slice, periodic_slice, aperiodic_slice

