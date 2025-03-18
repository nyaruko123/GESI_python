import numpy as np
from scipy.signal import decimate
from scipy.signal import lfilter


def Harvest(x, fs, option=None):
    # F0 估计通过 Harvest
    # f0_parameter = Harvest(x, fs, option);
    # f0_parameter = Harvest(x, fs);
    #
    # 输入
    #   x  : 输入信号
    #   fs : 采样频率
    #   option : 用户可以设置 f0_floor (Hz), f0_ceil (Hz) 和 frame_period (ms).
    #
    # 输出
    #   f0_parameter : f0 信息

    # 设置参数
    if option is None:
        f0_floor, f0_ceil, frame_period = SetDefaultParameters([])
    else:
        f0_floor, f0_ceil, frame_period = SetDefaultParameters(option)

    # 可以控制帧周期。但估计是以 1 ms 的帧移进行的。
    basic_frame_period = 1
    target_fs = 8000

    basic_temporal_positions = np.arange(0, len(x) / fs, basic_frame_period / 1000)

    channels_in_octave = 40
    # 扩展频率范围以准确获得 F0 候选
    adjusted_f0_floor = f0_floor * 0.9
    adjusted_f0_ceil = f0_ceil * 1.1
    boundary_f0_list = adjusted_f0_floor * 2.0 ** (
        np.arange(1, np.ceil(np.log2(adjusted_f0_ceil / adjusted_f0_floor) * channels_in_octave) + 1) / channels_in_octave
    )

    # 下采样到 target_fs Hz
    y, actual_fs = GetDownsampledSignal(x, fs, target_fs)
    fft_size = 2 ** np.ceil(np.log2(len(y) + 5 + 2 * np.floor(actual_fs / boundary_f0_list[0] * 2)))
    y_spectrum = np.fft.fft(y, int(fft_size))

    raw_f0_candidates = GetRawF0Candidates(
        len(basic_temporal_positions), boundary_f0_list, len(y),
        basic_temporal_positions, actual_fs, y_spectrum, f0_floor, f0_ceil
    )

    f0_candidates, number_of_candidates = DetectOfficialF0Candidates(raw_f0_candidates)
    f0_candidates = OverlapF0Candidates(f0_candidates, number_of_candidates)

    f0_candidates, f0_scores = RefineCandidates(y, actual_fs,
        basic_temporal_positions, f0_candidates, f0_floor, f0_ceil)

    f0_candidates, f0_scores = RemoveUnreliableCandidates(f0_candidates, f0_scores)

    connected_f0, vuv = FixF0Contour(f0_candidates, f0_scores)
    smoothed_f0 = SmoothF0Contour(connected_f0)

    temporal_positions = np.arange(0, len(x) / fs, frame_period / 1000)
    f0_parameter = {}
    f0_parameter['temporal_positions'] = temporal_positions

    # 将 temporal_positions 转换为索引：
    # 这里假设索引 = np.round(temporal_positions * 1000).astype(int) + 1
    # 为避免超出范围，用 np.clip 限制在 0 到 len(smoothed_f0)-1 之间
    indices = np.clip(np.round(temporal_positions * 1000).astype(int) + 1, 0, len(smoothed_f0) - 1)
    f0_parameter['f0'] = smoothed_f0[indices]
    f0_parameter['vuv'] = vuv[indices]
    f0_parameter['f0_candidates'] = f0_candidates
    # f0_parameter['f0_candidates_score'] = f0_candidates_score;

    return f0_parameter

def RemoveUnreliableCandidates(f0_candidates, f0_scores):
    # 移除不可靠的 F0 候选
    new_f0_candidates = f0_candidates.copy()
    new_f0_scores = f0_scores.copy()
    threshold = 0.05

    f0_length = new_f0_candidates.shape[1]  # F0 候选的长度
    number_of_candidates = new_f0_candidates.shape[0]  # 候选的数量

    for i in range(1, f0_length - 1):  # 从第二个到倒数第二个
        for j in range(number_of_candidates):  # 遍历所有候选
            reference_f0 = new_f0_candidates[j, i]
            if reference_f0 == 0:
                continue  # 如果参考 F0 为 0，则跳过

            # 选择最佳 F0
            _, min_error1 = SelectBestF0(reference_f0, new_f0_candidates[:, i + 1], 1)
            _, min_error2 = SelectBestF0(reference_f0, new_f0_candidates[:, i - 1], 1)
            min_error = min(min_error1, min_error2)

            if min_error > threshold:
                new_f0_candidates[j, i] = 0  # 标记为不可靠
                new_f0_scores[j, i] = 0  # 标记为不可靠

    return new_f0_candidates, new_f0_scores


def OverlapF0Candidates(f0_candidates, max_candidates):
    """
    对 F0 候选进行重叠处理。
    这里采用简单的策略：将 f0_candidates 沿行方向重复 number_of_candidates 次，
    其中 number_of_candidates = n * 2 + 1，n 固定为 3。
    """
    n = 3  # 优化参数
    number_of_candidates = n * 2 + 1  # 7
    # f0_candidates 假设形状为 (max_candidates, L)
    new_f0_candidates = np.tile(f0_candidates, (number_of_candidates, 1))
    return new_f0_candidates


def GetDownsampledSignal(x, fs, target_fs):
    # 获取下采样信号
    decimation_ratio = round(fs / target_fs)

    # 确保 x 是二维数组（列向量）
    if x.ndim == 1:
        x = x.reshape(-1, 1)

    if fs <= target_fs:
        y = x[:, 0]  # 如果目标采样率大于等于原采样率，直接返回原信号
        actual_fs = fs
    else:
        offset = np.ceil(140 / decimation_ratio) * decimation_ratio
        # 延伸信号：用二维数组保证维度匹配
        x = np.concatenate((np.ones((int(offset), 1)) * x[0], x, np.ones((int(offset), 1)) * x[-1]), axis=0)
        # 对 x 的第一列进行下采样
        y0 = decimate(x[:, 0], decimation_ratio, ftype='fir')  # 使用 FIR 滤波器进行下采样
        actual_fs = fs / decimation_ratio
        y = y0[int(offset / decimation_ratio): -int(offset / decimation_ratio)]

    y = y - np.mean(y)  # 去除均值
    return y, actual_fs


def GetRawF0Candidates(number_of_frames, boundary_f0_list, y_length, temporal_positions, actual_fs, y_spectrum, f0_floor, f0_ceil):
    # 获取原始 F0 候选
    raw_f0_candidates = np.zeros((len(boundary_f0_list), number_of_frames))  # 初始化原始 F0 候选矩阵

    for i in range(len(boundary_f0_list)):
        raw_f0_candidates[i, :] = GetF0CandidateFromRawEvent(boundary_f0_list[i], actual_fs, y_spectrum, y_length, temporal_positions, f0_floor, f0_ceil)

    return raw_f0_candidates



def DetectOfficialF0Candidates(raw_f0_candidates):
    # 检测官方 F0 候选
    number_of_channels, number_of_frames = raw_f0_candidates.shape
    f0_candidates = np.zeros((round(number_of_channels / 10), number_of_frames))  # 初始化 F0 候选矩阵

    number_of_candidates = 0
    threshold = 10

    for i in range(number_of_frames):
        tmp = raw_f0_candidates[:, i]
        tmp[tmp > 0] = 1  # 将大于 0 的值设为 1
        tmp[0] = 0  # 第一个元素设为 0
        tmp[-1] = 0  # 最后一个元素设为 0
        tmp = np.diff(tmp)  # 计算差分
        st = np.where(tmp == 1)[0]  # 找到开始位置
        ed = np.where(tmp == -1)[0]  # 找到结束位置
        count = 1

        for j in range(len(st)):
            dif = ed[j] - st[j]
            if dif >= threshold:  # 如果持续时间大于阈值
                tmp_f0 = raw_f0_candidates[st[j] + 1: ed[j], i]  # 获取 F0 候选
                f0_candidates[count - 1, i] = np.mean(tmp_f0)  # 计算均值
                count += 1

        number_of_candidates = max(number_of_candidates, count - 1)  # 更新候选数量

    return f0_candidates, number_of_candidates



def GetF0CandidateFromRawEvent(boundary_f0, fs, y_spectrum, y_length, temporal_positions, f0_floor, f0_ceil):
    # 从原始事件获取 F0 候选
    filter_length_half = round(fs / boundary_f0 * 2)  # 计算滤波器长度的一半
    band_pass_filter_base = nuttall(filter_length_half * 2 + 1)  # Nuttall 窗函数
    shifter = np.cos(2 * np.pi * boundary_f0 * np.arange(-filter_length_half, filter_length_half + 1) / fs)  # 频移
    band_pass_filter = band_pass_filter_base[:, np.newaxis] * shifter[:, np.newaxis]  # 带通滤波器

    index_bias = filter_length_half + 1
    spectrum_low_pass_filter = np.fft.fft(band_pass_filter, n=len(y_spectrum))  # 计算低通滤波器的 FFT

    filtered_signal = np.real(np.fft.ifft(spectrum_low_pass_filter * y_spectrum))  # 反 FFT 得到滤波信号

    # 构造要截取的索引数组：从 index_bias 开始，长度为 y_length
    indices = index_bias + np.arange(y_length)  # 生成索引数组
    # 如果最后一个索引超出数组边界，则只取有效部分
    if indices[-1] >= filtered_signal.shape[0]:
        indices = indices[indices < filtered_signal.shape[0]]
    filtered_signal = filtered_signal[indices]

    # 计算四种事件
    negative_zero_cross = ZeroCrossingEngine(filtered_signal, fs)  # 负零交叉
    positive_zero_cross = ZeroCrossingEngine(-filtered_signal, fs)  # 正零交叉
    peak = ZeroCrossingEngine(np.diff(filtered_signal), fs)  # 峰值
    dip = ZeroCrossingEngine(-np.diff(filtered_signal), fs)  # 谷值

    f0_candidate = GetF0CandidateContour(negative_zero_cross, positive_zero_cross, peak, dip, temporal_positions)  # 获取 F0 候选轮廓

    # 移除不可靠的候选
    f0_candidate[f0_candidate > boundary_f0 * 1.1] = 0
    f0_candidate[f0_candidate < boundary_f0 * 0.9] = 0
    f0_candidate[f0_candidate > f0_ceil] = 0
    f0_candidate[f0_candidate < f0_floor] = 0

    return f0_candidate




def GetF0CandidateContour(negative_zero_cross, positive_zero_cross, peak, dip, temporal_positions):
    # 获取 F0 候选轮廓
    usable_channel = (
        max(0, len(negative_zero_cross["interval_locations"]) - 2) *
        max(0, len(positive_zero_cross["interval_locations"]) - 2) *
        max(0, len(peak["interval_locations"]) - 2) *
        max(0, len(dip["interval_locations"]) - 2)
    )

    if usable_channel > 0:
        interpolated_f0_list = np.zeros((4, len(temporal_positions)))  # 初始化插值 F0 列表

        # 线性插值
        interpolated_f0_list[0, :] = np.interp(
            temporal_positions,
            negative_zero_cross["interval_locations"],
            negative_zero_cross["interval_based_f0"],
            left=np.nan,
            right=np.nan
        )  # 负零交叉的插值

        interpolated_f0_list[1, :] = np.interp(
            temporal_positions,
            positive_zero_cross["interval_locations"],
            positive_zero_cross["interval_based_f0"],
            left=np.nan,
            right=np.nan
        )  # 正零交叉的插值

        interpolated_f0_list[2, :] = np.interp(
            temporal_positions,
            peak["interval_locations"],
            peak["interval_based_f0"],
            left=np.nan,
            right=np.nan
        )  # 峰值的插值

        interpolated_f0_list[3, :] = np.interp(
            temporal_positions,
            dip["interval_locations"],
            dip["interval_based_f0"],
            left=np.nan,
            right=np.nan
        )  # 谷值的插值

        f0_candidate = np.mean(interpolated_f0_list, axis=0)  # 计算均值
    else:
        f0_candidate = temporal_positions * 0  # 如果没有可用通道，返回零

    return f0_candidate


def ZeroCrossingEngine(x, fs):
    # 如果 x 是多维，则扁平化为一维数组
    if x.ndim > 1:
        x = x.flatten()
        
    # 计算相邻样本的比较结果
    comp = np.concatenate((x[1:], [x[-1]]))
    condition = (comp * x < 0) * (comp < x)
    
    # 生成1D索引数组（长度为 len(x)），并与条件逐元素相乘
    indices_all = np.arange(1, len(x) + 1)
    negative_going_points = indices_all * condition  # 1D 数组
    
    # 过滤出有效的交叉点
    edge_list = negative_going_points[negative_going_points > 0]
    
    # 精细化交叉点计算：利用线性插值修正交叉位置
    indices = edge_list.astype(int)
    # 为防止越界，确保最后一个索引小于 len(x)-1
    valid = indices < (len(x) - 1)
    indices = indices[valid]
    
    # 计算相邻采样值之差
    diff = x[indices + 1] - x[indices]
    # 替换 diff 中为0的部分，避免除零
    diff[diff == 0] = 1e-8
    fine_edge_list = edge_list[valid] - x[indices] / diff
    
    event_struct = {}
    if fine_edge_list.size < 2:
        event_struct['interval_locations'] = np.array([])
        event_struct['interval_based_f0'] = np.array([])
    else:
        event_struct['interval_locations'] = (fine_edge_list[:-1] + fine_edge_list[1:]) / 2 / fs
        # 计算相邻交叉点的间隔
        delta = np.diff(fine_edge_list)
        # 替换 delta 中为0的部分，避免除零
        delta[delta == 0] = 1e-8
        event_struct['interval_based_f0'] = fs / delta
    
    return event_struct


def FixF0Contour(f0_candidates, f0_scores):
    # 修正 F0 轮廓
    f0_base = SearchF0Base(f0_candidates, f0_scores)  # 搜索 F0 基础
    f0_step1 = FixStep1(f0_base, 0.008)  # 第一步修正，优化参数为 0.008
    f0_step2 = FixStep2(f0_step1, 6)  # 第二步修正，优化参数为 6
    f0_step3 = FixStep3(f0_step2, f0_candidates, 0.18, f0_scores)  # 第三步修正，优化参数为 0.18
    f0 = FixStep4(f0_step3, 9)  # 第四步修正，优化参数为 9

    vuv = f0  # 将 F0 赋值给 vuv
    vuv[vuv != 0] = 1  # 将非零值设为 1

    return f0, vuv  # 返回 F0 和 vuv



def SearchF0Base(f0_candidates, f0_scores):
    # 选择得分最高的 F0 作为基本 F0 轮廓
    f0_base = np.zeros(f0_candidates.shape[1])  # 初始化 f0_base 为零数组
    for i in range(len(f0_base)):
        # 找到得分最高的索引
        max_index = np.argmax(f0_scores[:, i])
        f0_base[i] = f0_candidates[max_index, i]  # 选择对应的 F0 值

    return f0_base  # 返回基本 F0 轮廓



def FixStep1(f0_base, allowed_range):
    # 步骤 1：快速变化的 F0 轮廓被替换为 0
    f0_step1 = f0_base.copy()  # 复制 f0_base
    f0_step1[0] = 0  # 将第一个元素设为 0
    f0_step1[1] = 0  # 将第二个元素设为 0

    for i in range(2, len(f0_base)):  # 从第三个元素开始遍历
        if f0_base[i] == 0:
            continue  # 如果当前 F0 为 0，跳过

        # 计算参考 F0
        reference_f0 = f0_base[i - 1] * 2 - f0_base[i - 2]

        # 检查变化是否超过允许范围
        if (abs((f0_base[i] - reference_f0) / reference_f0) > allowed_range and
                abs((f0_base[i] - f0_base[i - 1]) / f0_base[i - 1]) > allowed_range):
            f0_step1[i] = 0  # 替换为 0

    return f0_step1  # 返回修正后的 F0 轮廓



def FixStep2(f0_step1, voice_range_minimum):
    # 步骤 2：去除周期较短的有声部分
    f0_step2 = f0_step1.copy()  # 复制 f0_step1
    boundary_list = GetBoundaryList(f0_step1)  # 获取边界列表

    for i in range(len(boundary_list) // 2):  # 遍历边界列表
        distance = boundary_list[2 * i + 1] - boundary_list[2 * i]  # 计算距离
        if distance < voice_range_minimum:
            # 将短于最小范围的部分设为 0
            f0_step2[boundary_list[2 * i]:boundary_list[2 * i + 1] + 1] = 0

    return f0_step2  # 返回修正后的 F0 轮廓


def FixStep3(f0_step2, f0_candidates, allowed_range, f0_scores):
    # 步骤 3：根据 F0 轮廓的连续性扩展有声部分
    f0_step3 = f0_step2.copy()  # 复制 f0_step2
    boundary_list = GetBoundaryList(f0_step2)  # 获取边界列表
    multi_channel_f0 = GetMultiChannelF0(f0_step2, boundary_list)  # 获取多通道 F0
    range_ = np.zeros((len(boundary_list) // 2, 2))  # 初始化范围数组
    threshold1 = 100  # 阈值 1
    threshold2 = 2200  # 阈值 2

    count = 0
    for i in range(len(boundary_list) // 2):  # 遍历边界列表
        # 扩展 F0
        extended_f0, tmp_range = ExtendF0(multi_channel_f0[i, :],
                                          boundary_list[i * 2],
                                          min(len(f0_step2) - 1,
                                              boundary_list[i * 2] + threshold1),
                                          1, f0_candidates, allowed_range)

        tmp_f0_sequence, tmp_range[0] = ExtendF0(extended_f0,
                                                 boundary_list[(i * 2) - 1],
                                                 max(2, boundary_list[(i * 2) - 1] - threshold1),
                                                 -1, f0_candidates, allowed_range)

        mean_f0 = np.mean(tmp_f0_sequence[int(tmp_range[0]):int(tmp_range[1])])  # 计算均值 F0
        if threshold2 / mean_f0 < tmp_range[1] - tmp_range[0]:  # 检查条件
            count += 1
            multi_channel_f0[count - 1, :] = tmp_f0_sequence  # 更新多通道 F0
            range_[count - 1, :] = tmp_range  # 更新范围

    multi_channel_f0 = multi_channel_f0[:count, :]  # 截取有效的多通道 F0
    range_ = range_[:count, :]  # 截取有效的范围

    if count > 0:
        f0_step3 = MergeF0(multi_channel_f0, range_, f0_candidates, f0_scores)  # 合并 F0

    return f0_step3  # 返回修正后的 F0 轮廓



def GetMultiChannelF0(f0, boundary_list):
    # 获取多通道 F0
    multi_channel_f0 = np.zeros((len(boundary_list) // 2, len(f0)))  # 初始化多通道 F0 数组
    for i in range(len(boundary_list) // 2):  # 遍历边界列表
        # 将 F0 的相应部分赋值到多通道 F0
        multi_channel_f0[i, boundary_list[2 * i]:boundary_list[2 * i + 1] + 1] = \
            f0[boundary_list[2 * i]:boundary_list[2 * i + 1] + 1]

    return multi_channel_f0  # 返回多通道 F0


def MergeF0(multi_channel_f0, range_, f0_candidates, f0_scores):
    # 合并 F0
    number_of_channels = multi_channel_f0.shape[0]  # 获取通道数量
    sorted_order = np.argsort(range_[:, 0])  # 根据范围的起始值排序
    f0 = multi_channel_f0[sorted_order[0], :]  # 初始化 f0

    for i in range(1, number_of_channels):  # 从第二个通道开始遍历
        # 无重叠情况
        if range_[sorted_order[i], 0] - range_[sorted_order[0], 1] > 0:
            f0[range_[sorted_order[i], 0]:range_[sorted_order[i], 1] + 1] = \
                multi_channel_f0[sorted_order[i], range_[sorted_order[i], 0]:range_[sorted_order[i], 1] + 1]
            range_[sorted_order[0], 0] = range_[sorted_order[i], 0]
            range_[sorted_order[0], 1] = range_[sorted_order[i], 1]
        else:  # 有重叠情况
            f0, range_[sorted_order[0], 1] = \
                MergeF0Sub(f0, range_[sorted_order[0], 0], range_[sorted_order[0], 1],
                            multi_channel_f0[sorted_order[i], :], range_[sorted_order[i], 0],
                            range_[sorted_order[i], 1], f0_candidates, f0_scores)

    return f0  # 返回合并后的 F0


def MergeF0Sub(f0_1, st1, ed1, f0_2, st2, ed2, f0_candidates, f0_scores):
    # 合并两个 F0 轮廓
    merged_f0 = f0_1.copy()  # 复制 f0_1

    # 完全重叠的部分
    if st1 <= st2 and ed1 >= ed2:
        new_ed = ed1
        return merged_f0, new_ed  # 返回合并后的 F0 和新的结束位置

    new_ed = ed2  # 更新新的结束位置

    score1 = 0
    score2 = 0
    for i in range(st2, ed1 + 1):  # 遍历重叠区间
        score1 += SerachScore(f0_1[i], f0_candidates[:, i], f0_scores[:, i])  # 计算 f0_1 的得分
        score2 += SerachScore(f0_2[i], f0_candidates[:, i], f0_scores[:, i])  # 计算 f0_2 的得分

    if score1 > score2:
        merged_f0[ed1:ed2 + 1] = f0_2[ed1:ed2 + 1]  # 用 f0_2 的值替换
    else:
        merged_f0[st2:ed2 + 1] = f0_2[st2:ed2 + 1]  # 用 f0_2 的值替换

    return merged_f0, new_ed  # 返回合并后的 F0 和新的结束位置

def SerachScore(f0, f0_candidates, f0_scores):
    # 计算得分
    score = 0
    for i in range(len(f0_candidates)):  # 遍历所有候选 F0
        if f0 == f0_candidates[i] and score < f0_scores[i]:  # 如果匹配且得分更高
            score = f0_scores[i]  # 更新得分
    return score  # 返回最终得分


def ExtendF0(f0, origin, last_point, shift, f0_candidates, allowed_range):
    # 扩展 F0 值
    threshold = 4  # 阈值
    extended_f0 = f0.copy()  # 复制原始 F0
    tmp_f0 = extended_f0[origin]  # 获取初始值
    shifted_origin = origin  # 初始化偏移原点

    count = 0  # 计数器
    for i in range(origin, last_point + 1, shift):  # 遍历从 origin 到 last_point 的范围
        extended_f0[i + shift] = SelectBestF0(tmp_f0, f0_candidates[:, i + shift], allowed_range)  # 选择最佳 F0
        if extended_f0[i + shift] != 0:  # 如果选择的 F0 不为 0
            tmp_f0 = extended_f0[i + shift]  # 更新临时 F0
            count = 0  # 重置计数器
            shifted_origin = i + shift  # 更新偏移原点
        else:
            count += 1  # 计数器加一
        if count == threshold:  # 如果达到阈值，退出循环
            break

    return extended_f0, shifted_origin  # 返回扩展后的 F0 和偏移原点

def SelectBestF0(reference_f0, f0_candidates, allowed_range):
    # 选择最佳 F0 值
    best_f0 = 0  # 初始化最佳 F0
    best_error = allowed_range  # 初始化最佳误差

    for i in range(len(f0_candidates)):  # 遍历所有候选 F0
        tmp = abs(reference_f0 - f0_candidates[i]) / reference_f0  # 计算误差比例
        if tmp > best_error:  # 如果误差超过允许范围，继续下一个候选
            continue
        best_f0 = f0_candidates[i]  # 更新最佳 F0
        best_error = tmp  # 更新最佳误差

    return best_f0, best_error  # 返回最佳 F0 和最佳误差

def FixStep4(f0_step3, threshold):
    # 修正步骤 4：在短的无声段中伪造 F0
    f0_step4 = f0_step3.copy()  # 复制 F0 数据
    boundary_list = GetBoundaryList(f0_step3)  # 获取边界列表

    for i in range(int(len(boundary_list) / 2) - 1):  # 遍历边界列表
        distance = boundary_list[(2 * i) + 1] - boundary_list[2 * i] - 1  # 计算距离
        if distance >= threshold:  # 如果距离大于等于阈值，继续下一个
            continue
        boundary0 = f0_step3[boundary_list[2 * i]] + 1  # 从 0 到 1 的无声段
        boundary1 = f0_step3[boundary_list[(2 * i) + 1]] - 1  # 从 1 到 0 的无声段
        coefficient = (boundary1 - boundary0) / (distance + 1)  # 计算系数
        count = 1
        for j in range(boundary_list[2 * i] + 1, boundary_list[(2 * i) + 1]):  # 填充 F0
            f0_step4[j] = boundary0 + coefficient * count  # 更新 F0 值
            count += 1  # 计数器加一

    return f0_step4  # 返回修正后的 F0

def RefineCandidates(x, fs, temporal_positions, f0_candidates, f0_floor, f0_ceil):
    # 精炼 F0 候选值
    new_f0_candidates = f0_candidates.copy()  # 复制 F0 候选值
    f0_scores = f0_candidates * 0  # 初始化 F0 分数

    for i in range(len(temporal_positions)):  # 遍历时间位置
        for j in range(f0_candidates.shape[0]):  # 遍历 F0 候选值
            tmp_f0 = f0_candidates[j, i]  # 获取当前 F0 候选值
            if tmp_f0 == 0:  # 如果 F0 候选值为 0，继续下一个
                continue
            new_f0_candidates[j, i], f0_scores[j, i] = GetRefinedF0(x, fs, temporal_positions[i], tmp_f0, f0_floor, f0_ceil)  # 调用未定义函数

    return new_f0_candidates, f0_scores  # 返回新的 F0 候选值和 F0 分数



def GetRefinedF0(x, fs, current_time, current_f0, f0_floor, f0_ceil):
    # 获取精炼的 F0 和分数
    half_window_length = int(np.ceil(3 * fs / current_f0 / 2))  # 半窗口长度
    window_length_in_time = (2 * half_window_length + 1) / fs  # 窗口长度（时间）
    base_time = np.arange(-half_window_length, half_window_length + 1) / fs  # 基础时间
    fft_size = 2 ** int(np.ceil(np.log2((half_window_length * 2 + 1) + 1)))  # FFT 大小
    frequency_axis = np.arange(fft_size) / fft_size * fs  # 频率轴

    # 急救处理：计算基于时间的索引（注意：将索引转换为 0 基准）
    base_index = np.round((current_time + base_time) * fs + 0.001).astype(int)
    # 将 base_index 从 1-based 转为 0-based：clip 时上界为 len(x)，减 1 后有效索引范围 0 ~ len(x)-1
    safe_index = np.clip(base_index, 1, len(x)) - 1

    index_time = safe_index / fs  # 使用安全索引转换为时间（已经是 0 基准）
    window_time = index_time - current_time  # 窗口时间
    main_window = (
        0.42 + 0.5 * np.cos(2 * np.pi * window_time / window_length_in_time) +
        0.08 * np.cos(4 * np.pi * window_time / window_length_in_time)
    )  # 主窗口
    diff_window = -(np.diff(np.concatenate(([0], main_window))) + np.diff(np.concatenate((main_window, [0])))) / 2  # 差分窗口

    # 计算频谱时，取 x[safe_index]
    spectrum = np.fft.fft(x[safe_index] * main_window, fft_size)  # 频谱
    diff_spectrum = np.fft.fft(x[safe_index] * diff_window, fft_size)  # 差分频谱
    numerator_i = np.real(spectrum) * np.imag(diff_spectrum) - np.imag(spectrum) * np.real(diff_spectrum)  # 分子
    power_spectrum = np.abs(spectrum) ** 2  # 功率谱
    # 注意：如果 power_spectrum 中有 0，这里会产生除零问题
    instantaneous_frequency = frequency_axis + numerator_i / power_spectrum * fs / (2 * np.pi)  # 瞬时频率

    number_of_harmonics = min(np.floor(fs / (2 * current_f0)), 6)  # 确保安全
    harmonics_index = np.arange(1, number_of_harmonics + 1)  # 谐波索引
    # 修改：不再加1，因为索引从0开始
    index_list = np.round(current_f0 * fft_size / fs * harmonics_index).astype(int)
    # 确保 index_list 中的值不超出频谱范围
    index_list = np.clip(index_list, 0, len(instantaneous_frequency)-1)
    
    instantaneous_frequency_list = instantaneous_frequency[index_list]  # 瞬时频率列表
    amplitude_list = np.sqrt(power_spectrum[index_list])  # 振幅列表
    refined_f0 = np.sum(amplitude_list * instantaneous_frequency_list) / np.sum(amplitude_list * harmonics_index)  # 精炼的 F0

    variation = np.abs((instantaneous_frequency_list / harmonics_index - current_f0) / current_f0)  # 变化
    refined_score = 1 / (1e-12 + np.mean(variation))  # 精炼分数
    if refined_f0 < f0_floor or refined_f0 > f0_ceil or refined_score < 2.5:  # 检查条件
        refined_f0 = 0
        refined_score = 0

    return refined_f0, refined_score


def SmoothF0Contour(f0):
    # 平滑 F0 轮廓
    b = [0.0078202080334971724, 0.015640416066994345, 0.0078202080334971724]  # 滤波器系数 b
    a = [1.0, -1.7347257688092754, 0.76600660094326412]  # 滤波器系数 a

    smoothed_f0 = np.concatenate((np.zeros(300), f0, np.zeros(300)))  # 在 F0 前后添加零
    boundary_list = GetBoundaryList(smoothed_f0)  # 获取边界列表
    multi_channel_f0 = GetMultiChannelF0(smoothed_f0, boundary_list)  # 获取多通道 F0

    for i in range(len(boundary_list) // 2):  # 遍历边界列表
        tmp_f0_contour = FilterF0Contour(multi_channel_f0[i, :],  # 过滤 F0 轮廓
                                          boundary_list[(i * 2) - 1],
                                          boundary_list[i * 2],
                                          b, a)
        smoothed_f0[boundary_list[(i * 2) - 1]: boundary_list[i * 2]] = \
            tmp_f0_contour[boundary_list[(i * 2) - 1]: boundary_list[i * 2]]

    smoothed_f0 = smoothed_f0[300: -300]  # 去掉前后的零

    return smoothed_f0  # 返回平滑后的 F0



def FilterF0Contour(f0, st, ed, b, a):
    # 过滤 F0 轮廓
    smoothed_f0 = f0.copy()  # 复制 F0
    smoothed_f0[0:st - 1] = smoothed_f0[st]  # 用 st 位置的值填充前面的部分
    smoothed_f0[ed:] = smoothed_f0[ed]  # 用 ed 位置的值填充后面的部分

    aaa = lfilter(b, a, smoothed_f0)  # 应用滤波器
    bbb = lfilter(b, a, aaa[::-1])  # 反向滤波
    smoothed_f0 = bbb[::-1]  # 再次反向

    smoothed_f0[0:st - 1] = 0  # 将前面的部分置为 0
    smoothed_f0[ed:] = 0  # 将后面的部分置为 0

    return smoothed_f0  # 返回平滑后的 F0


def nuttall(N):
    # 生成 Nuttall 窗函数
    t = np.arange(N).reshape(-1, 1) * 2 * np.pi / (N - 1)  # 生成时间向量
    coefs = np.array([0.355768, -0.487396, 0.144232, -0.012604])  # 窗函数系数
    window = np.cos(t * np.array([0, 1, 2, 3])) @ coefs  # 计算窗函数

    return window  # 返回 Nuttall 窗函数


def GetBoundaryList(f0):
    # 获取边界列表
    vuv = f0.copy()  # 复制 F0
    vuv[vuv != 0] = 1  # 非零值设为 1
    vuv[0] = 0  # 第一个元素设为 0
    vuv[-1] = 0  # 最后一个元素设为 0
    diff_vuv = np.diff(vuv)  # 计算差分
    boundary_list = np.where(diff_vuv != 0)[0]  # 找到边界
    boundary_list[0::2] += 1  # 偶数索引加 1

    return boundary_list  # 返回边界列表

def SetDefaultParameters(option):
    # 设置默认参数
    f0_floor = 71  # F0 下限
    f0_ceil = 800  # F0 上限
    frame_period = 5  # 帧周期

    if option is not None:  # 检查 option 是否为空
        if 'f0_floor' in option:  # 检查是否有 f0_floor 字段
            f0_floor = option['f0_floor']  # 更新 f0_floor
        if 'f0_ceil' in option:  # 检查是否有 f0_ceil 字段
            f0_ceil = option['f0_ceil']  # 更新 f0_ceil
        if 'frame_period' in option:  # 检查是否有 frame_period 字段
            frame_period = option['frame_period']  # 更新 frame_period

    return f0_floor, f0_ceil, frame_period  # 返回参数
