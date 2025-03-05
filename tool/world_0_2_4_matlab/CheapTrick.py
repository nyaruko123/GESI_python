# // function spectrum_paramter = CheapTrick(x, fs, source_object, option)
# // % Spectral envelope extraction based on an algorithm, CheapTrick.
# // % spectrum_paramter = CheapTrick(x, fs, source_object, option);
# // % spectrum_paramter = CheapTrick(x, fs, source_object);
# // %
# // % Input
# // %   x  : input signal
# // %   fs : sampling frequency
# // %   source_object : source information object
# // %   option    : It has two parameters (q1 and fft_size)
# // %               Parameter q1 is used for the spectral recovery.
# // %               The lowest F0 that WORLD can work as expected is determined
# // %               by the following: 3.0 * fs / fft_size.
# // %
# // % Output
# // %   spectrum_paramter : spectum infromation
# // %
# // % Caution: WORLD is not corresponding with TANDEM-STRAIGHT.
# // %          However, the difference is only the name of several parameters.
# // %
# // % 2014/04/29: First version was released.
# // % 2015/09/22: A parameter (q1) is controllable.
# // % 2016/12/28: Refactoring (default value of q1 was modified. -0.09 -> -0.15)
# // % 2017/01/02: A parameter fft_size is controllable.
# // % 2017/01/29: A bug was fixed.
# // % 2017/05/20: A safeguard was added.
#
# // % set default parameters
# // f0_low_limit = 71;
# // default_f0 = 500;
# // fft_size = 2 ^ ceil(log2(3 * fs / f0_low_limit + 1));
# // q1 = -0.15;
# // if nargin == 4
# //   if isfield(option, 'q1') == 1
# //     q1 = option.q1;
# //   end;
# //   if isfield(option, 'fft_size') == 1
# //     fft_size = option.fft_size;
# //   end;
# // end;
# // f0_low_limit = fs * 3.0 / (fft_size - 3.0);
#
# // temporal_positions = source_object.temporal_positions;
# // f0_sequence = source_object.f0;
# // if isfield(source_object, 'vuv')
# //   f0_sequence(source_object.vuv == 0) = default_f0;
# // end;
#
# // spectrogram = zeros(fft_size / 2 + 1, length(f0_sequence));
# // for i = 1:length(f0_sequence)
# //   if f0_sequence(i) < f0_low_limit; f0_sequence(i) = default_f0; end;
# //   spectrogram(:,i) = EstimateOneSlice(x, fs, f0_sequence(i),...
# //     temporal_positions(i), fft_size, q1);
# // end;
#
# // % output parameters
# // spectrum_paramter.temporal_positions = temporal_positions;
# // spectrum_paramter.spectrogram = spectrogram;
# // spectrum_paramter.fs = fs;
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function spectral_envelope  = EstimateOneSlice(x, fs, current_f0,...
# //   current_position, fft_size, q1)
# // waveform = GetWindowedWaveform(x, fs, current_f0, current_position);
# // power_spectrum = GetPowerSpectrum(waveform, fs, fft_size, current_f0);
# // smoothed_spectrum = LinearSmoothing(power_spectrum, current_f0, fs, fft_size);
# // spectral_envelope = SmoothingWithRecovery(...
# //   [smoothed_spectrum; smoothed_spectrum(end - 1 : -1 : 2)], current_f0, fs,...
# //   fft_size, q1);
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function power_spectrum = GetPowerSpectrum(waveform, fs, fft_size, f0)
# // power_spectrum = abs(fft(waveform(:), fft_size)) .^ 2;
# // % DC correction
# // frequency_axis = (0 : fft_size - 1)' / fft_size * fs;
# // low_frequency_axis = frequency_axis(frequency_axis <  f0 + fs / fft_size);
# // low_frequency_replica = interp1(f0 - low_frequency_axis,...
# //   power_spectrum(frequency_axis < f0 + fs / fft_size),...
# //   low_frequency_axis(:), 'linear', 'extrap');
# // power_spectrum(frequency_axis < f0) =...
# //   low_frequency_replica(frequency_axis < f0) +...
# //   power_spectrum(frequency_axis < f0);
# // power_spectrum(end : -1 : fft_size / 2 + 2) = power_spectrum(2 : fft_size / 2);
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function waveform = GetWindowedWaveform(x, fs, current_f0, current_position)
# // %  prepare internal variables
# // half_window_length = round(1.5 * fs / current_f0);
# // base_index = (-half_window_length : half_window_length)';
# // index = round(current_position * fs + 0.001) + 1 + base_index;
# // safe_index = min(length(x), max(1, round(index)));
#
# // %  wave segments and set of windows preparation
# // segment = x(safe_index);
# // time_axis = base_index / fs / 1.5;
# // window = 0.5 * cos(pi * time_axis * current_f0) + 0.5;
# // window = window / sqrt(sum(window .^ 2));
# // waveform = segment .* window - window * mean(segment .* window) / mean(window);
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function smoothed_spectrum = LinearSmoothing(power_spectrum, f0, fs, fft_size)
# // double_frequency_axis = (0 : 2 * fft_size - 1)' / fft_size * fs - fs;
# // double_spectrum = [power_spectrum; power_spectrum];
#
# // double_segment = cumsum(double_spectrum * (fs / fft_size));
# // center_frequency = (0 : fft_size / 2)' / fft_size * fs;
# // low_levels = interp1H(double_frequency_axis + fs / fft_size / 2,...
# //   double_segment, center_frequency - f0 / 3);
# // high_levels = interp1H(double_frequency_axis + fs / fft_size / 2,...
# //   double_segment, center_frequency + f0 / 3);
#
# // smoothed_spectrum = (high_levels - low_levels) * 1.5 / f0;
# // smoothed_spectrum =...
# //   smoothed_spectrum + abs(randn(length(smoothed_spectrum), 1)) * eps;
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // % This is the imprementation of a matlab function
# // function yi = interp1H(x, y, xi)
# // delta_x = x(2) - x(1);
# // xi = max(x(1), min(x(end), xi));
# // xi_base = floor((xi - x(1)) / delta_x);
# // xi_fraction = (xi - x(1)) / delta_x - xi_base;
# // delta_y = [diff(y); 0];
# // yi = y(xi_base + 1) + delta_y(xi_base + 1) .* xi_fraction;
#
# // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# // function spectral_envelope =...
# //   SmoothingWithRecovery(smoothed_spectrum, f0, fs, fft_size, q1)
# // quefrency_axis = (0 : fft_size - 1)' / fs;
# // smoothing_lifter = sin(pi * f0 * quefrency_axis) ./ (pi * f0 * quefrency_axis);
# // smoothing_lifter(fft_size / 2 + 2 : end) =...
# //   smoothing_lifter(fft_size / 2 : -1 : 2);
# // smoothing_lifter(1) = 1;
#
# // compensation_lifter =...
# //   (1 - 2 * q1) + 2 * q1 * cos(2 * pi * quefrency_axis * f0);
# // compensation_lifter(fft_size / 2 + 2 : end) =...
# //   compensation_lifter(fft_size / 2 : -1 : 2);
# // tandem_cepstrum = fft(log(smoothed_spectrum));
# // tmp_spectral_envelope =...
# //   exp(real(ifft(tandem_cepstrum .* smoothing_lifter .* compensation_lifter)));
# // spectral_envelope = tmp_spectral_envelope(1 : fft_size / 2 + 1);


import numpy as np
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d

def cheap_trick(x, fs, source_object, option=None):
    f0_low_limit = 71
    default_f0 = 500
    fft_size = 2 ** int(np.ceil(np.log2(3 * fs / f0_low_limit + 1)))
    q1 = -0.15
    
    if option:
        q1 = option.get('q1', q1)
        fft_size = option.get('fft_size', fft_size)
    
    f0_low_limit = fs * 3.0 / (fft_size - 3.0)
    temporal_positions = source_object['temporal_positions']
    f0_sequence = np.array(source_object['f0'])
    
    if 'vuv' in source_object:
        f0_sequence[source_object['vuv'] == 0] = default_f0
    
    spectrogram = np.zeros((fft_size // 2 + 1, len(f0_sequence)))
    
    for i in range(len(f0_sequence)):
        if f0_sequence[i] < f0_low_limit:
            f0_sequence[i] = default_f0
        spectrogram[:, i] = estimate_one_slice(x, fs, f0_sequence[i],
                                               temporal_positions[i], fft_size, q1)
    
    return {'temporal_positions': temporal_positions, 'spectrogram': spectrogram, 'fs': fs}

def estimate_one_slice(x, fs, current_f0, current_position, fft_size, q1):
    waveform = get_windowed_waveform(x, fs, current_f0, current_position)
    power_spectrum = get_power_spectrum(waveform, fs, fft_size, current_f0)
    smoothed_spectrum = linear_smoothing(power_spectrum, current_f0, fs, fft_size)
    return smoothing_with_recovery(smoothed_spectrum, current_f0, fs, fft_size, q1)

def get_power_spectrum(waveform, fs, fft_size, f0):
    power_spectrum = np.abs(fft(waveform, fft_size)) ** 2
    frequency_axis = np.arange(fft_size) / fft_size * fs
    low_freq_mask = frequency_axis < f0 + fs / fft_size
    
    interp_func = interp1d(f0 - frequency_axis[low_freq_mask], power_spectrum[low_freq_mask],
                           kind='linear', fill_value='extrapolate')
    low_frequency_replica = interp_func(frequency_axis[low_freq_mask])
    power_spectrum[low_freq_mask] += low_frequency_replica
    power_spectrum[-1:fft_size // 2:-1] = power_spectrum[1:fft_size // 2]
    return power_spectrum

def get_windowed_waveform(x, fs, current_f0, current_position):
    half_window_length = round(1.5 * fs / current_f0)
    base_index = np.arange(-half_window_length, half_window_length + 1)
    index = int(round(current_position * fs)) + base_index
    safe_index = np.clip(index, 0, len(x) - 1)
    
    segment = x[safe_index]
    time_axis = base_index / fs / 1.5
    window = 0.5 * np.cos(np.pi * time_axis * current_f0) + 0.5
    window /= np.sqrt(np.sum(window ** 2))
    waveform = segment * window - window * np.mean(segment * window) / np.mean(window)
    return waveform

def linear_smoothing(power_spectrum, f0, fs, fft_size):
    double_frequency_axis = np.arange(2 * fft_size) / fft_size * fs - fs
    double_spectrum = np.concatenate((power_spectrum, power_spectrum))
    double_segment = np.cumsum(double_spectrum * (fs / fft_size))
    
    center_frequency = np.arange(fft_size // 2 + 1) / fft_size * fs
    low_levels = interp1h(double_frequency_axis + fs / fft_size / 2, double_segment, center_frequency - f0 / 3)
    high_levels = interp1h(double_frequency_axis + fs / fft_size / 2, double_segment, center_frequency + f0 / 3)
    
    smoothed_spectrum = (high_levels - low_levels) * 1.5 / f0
    smoothed_spectrum += np.abs(np.random.randn(len(smoothed_spectrum))) * np.finfo(float).eps
    return smoothed_spectrum

def interp1h(x, y, xi):
    xi = np.clip(xi, x[0], x[-1])
    interp_func = interp1d(x, y, kind='linear', fill_value='extrapolate')
    return interp_func(xi)

def smoothing_with_recovery(smoothed_spectrum, f0, fs, fft_size, q1):
    quefrency_axis = np.arange(fft_size) / fs
    smoothing_lifter = np.sinc(f0 * quefrency_axis)
    smoothing_lifter[fft_size // 2 + 1:] = smoothing_lifter[fft_size // 2:0:-1]
    
    compensation_lifter = (1 - 2 * q1) + 2 * q1 * np.cos(2 * np.pi * quefrency_axis * f0)
    compensation_lifter[fft_size // 2 + 1:] = compensation_lifter[fft_size // 2:0:-1]
    
    tandem_cepstrum = fft(np.log(smoothed_spectrum))
    tmp_spectral_envelope = np.exp(np.real(ifft(tandem_cepstrum * smoothing_lifter * compensation_lifter)))
    return tmp_spectral_envelope[:fft_size // 2 + 1]
#从语音信号中提取谱包络