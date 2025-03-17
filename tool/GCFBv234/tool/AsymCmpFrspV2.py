import numpy as np
from scipy.interpolate import interp1d
from .Freq2ERB import Freq2ERB
from .MakeAsymCmpFiltersV2 import MakeAsymCmpFiltersV2

def interpolate_rows(data, new_length):
    """
    将 data 的每一行插值到 new_length 列。
    data.shape => (row, old_length)
    插值后 => (row, new_length)
    """
    row, old_length = data.shape
    x = np.linspace(0, 1, old_length)
    x_new = np.linspace(0, 1, new_length)
    out = np.zeros((row, new_length), dtype=data.dtype)
    for i in range(row):
        f = interp1d(x, data[i, :], kind='linear', fill_value='extrapolate')
        out[i, :] = f(x_new)
    return out

def force_2d_shape(mat, num_ch, n_frq_rsl, do_interpolate=True):
    """
    强制将 mat 变为 (num_ch, n_frq_rsl)：
    1) 若 mat.shape[0] != num_ch, 则截断或扩展(简单 tile)到 num_ch 行。
    2) 若 mat.shape[1] != n_frq_rsl, 则:
       - 若只有1列 => tile到 n_frq_rsl;
       - 若 do_interpolate=True => 线性插值到 n_frq_rsl 列;
       - 否则做简单截断或 tile 到 n_frq_rsl 列。
    """
    mat = np.array(mat, ndmin=2)
    # 先行数修正
    if mat.shape[0] > num_ch:
        mat = mat[:num_ch, :]
    elif mat.shape[0] < num_ch:
        repeat_r = (num_ch // mat.shape[0]) + 1
        mat = np.tile(mat, (repeat_r, 1))[:num_ch, :]

    # 再列数修正
    old_cols = mat.shape[1]
    if old_cols == n_frq_rsl:
        return mat
    elif old_cols == 1:
        # 若仅1列，则简单 tile
        mat = np.tile(mat, (1, n_frq_rsl))
        return mat
    else:
        if do_interpolate and old_cols > 1:
            # 线性插值到 n_frq_rsl
            mat = interpolate_rows(mat, n_frq_rsl)
        else:
            # 否则简单截断/扩展
            if old_cols > n_frq_rsl:
                mat = mat[:, :n_frq_rsl]
            else:
                rep_c = (n_frq_rsl // old_cols) + 1
                mat = np.tile(mat, (1, rep_c))[:, :n_frq_rsl]
        return mat

def AsymCmpFrspV2(Frs, fs, b, c, NfrqRsl=1024, NumFilt=4, SwCoef=1):
    """
    Compute the frequency response of an Asymmetric Compensation Filter.

    * 强制所有相关变量都对齐为 (NumCh, NfrqRsl), 以避免 broadcasting 错误.
    """

    # 0) 基本检查
    if NfrqRsl is None or NfrqRsl == 0:
        NfrqRsl = 1024
    if NumFilt != 4:
        raise ValueError('NumFilter should be 4.')

    # 1) 转为2D
    Frs = np.array(Frs, ndmin=2)
    b   = np.array(b, ndmin=2)
    c   = np.array(c, ndmin=2)
    NumCh = Frs.shape[0]  # 行数 = 通道数

    # 2) 构造 freq
    # 当 NfrqRsl >= 64，按 0 ~ fs/2 均匀生成
    if NfrqRsl >= 64:
        freq_1d = np.linspace(0, fs/2, NfrqRsl)  # shape (NfrqRsl,)
    else:
        # 若有其它逻辑,请自行调整;演示只保留
        freq_1d = Frs[0,:]  # shape(?)
    freq_2d = freq_1d[np.newaxis, :]  # 变为 (1, NfrqRsl)
    # 对齐成 (NumCh, NfrqRsl)
    freq_2d = force_2d_shape(freq_2d, NumCh, NfrqRsl, do_interpolate=False)

    # 3) 对齐 Frs => (NumCh, NfrqRsl)
    Frs = force_2d_shape(Frs, NumCh, NfrqRsl, do_interpolate=True)

    print(f'[DEBUG] freq shape: {freq_2d.shape}, Frs shape: {Frs.shape}')

    # 4) 计算 ERBw => (NumCh, NfrqRsl)
    _, ERBw = Freq2ERB(Frs)
    ERBw = force_2d_shape(ERBw, NumCh, NfrqRsl, do_interpolate=True)

    # 5) 若 SwCoef=1 => 计算滤波器系数; 否则走硬编码
    if SwCoef == 1:
        ACFcoef, _ = MakeAsymCmpFiltersV2(fs, Frs, b, c)

    # 初始化
    ACFFrsp = np.ones((NumCh, NfrqRsl), dtype=float)

    # 6) 循环计算4个二阶滤波器
    for nf in range(1, NumFilt+1):
        if SwCoef == 0:
            # 硬编码
            p0 = 2
            p1 = 1.7818*(1 - 0.0791*b)*(1 - 0.1655*np.abs(c))
            p2 = 0.5689*(1 - 0.1620*b)*(1 - 0.0857*np.abs(c))
            r = np.exp(-p1*((p0/1.0724)**(nf-1))*2*np.pi*b*ERBw/fs)
            delfr = ((p0*1.0724)**(nf-1)) * p2*c*b*ERBw
            phi = 2*np.pi*np.maximum(Frs + delfr, 0)/fs
            psy = 2*np.pi*np.maximum(Frs - delfr, 0)/fs

            # 生成 ap, bz => (NumCh, 3)
            # 这里要注意 r, phi 均是(NumCh, NfrqRsl), 但你只想要一个(NumCh,)?
            # 常见做法: 只挑选最后一个点 or 中心点 or ...
            # 此处为了演示：直接用 r.ravel()[0:NumCh] => 仅取第一个频率点
            # 你也可自己决定如何选 freqIndex
            # ================
            # DEMO: 全部取第一个频率点(列=0)
            r0   = r[:, 0]
            phi0 = phi[:, 0]
            psy0 = psy[:, 0]
            # ================
            ap = np.column_stack((
                np.ones(NumCh),
                -2 * r0 * np.cos(phi0),
                r0**2
            ))
            bz = np.column_stack((
                np.ones(NumCh),
                -2 * r0 * np.cos(psy0),
                r0**2
            ))
        else:
            # 来自 ACFcoef
            # shape => (ChCoef, 3, 4)? 具体看 MakeAsymCmpFiltersV2 的返回
            ap = ACFcoef['ap'][:, :, nf-1]  # => (something, 3)
            bz = ACFcoef['bz'][:, :, nf-1]  # => (something, 3)

        print(f'[DEBUG] bz shape: {bz.shape}, ap shape: {ap.shape}')

        # 若 shape(4096,3) 而 NumCh=2 => 截断
        if bz.shape[0] > NumCh:
            bz = bz[:NumCh, :]
            ap = ap[:NumCh, :]
        elif bz.shape[0] < NumCh:
            # 也可报错 or tile
            raise ValueError('Filter channel < actual channel!')

        # 7) 计算滤波器幅度响应
        cs1 = np.cos(2*np.pi*freq_2d/fs)  # => (NumCh, NfrqRsl)
        cs2 = np.cos(4*np.pi*freq_2d/fs)  # => (NumCh, NfrqRsl)

        # 先分子
        hb0 = (bz[:,0]**2 + bz[:,1]**2 + bz[:,2]**2)[:, np.newaxis]
        hb1 = (2*bz[:,1]*(bz[:,0]+bz[:,2]))[:, np.newaxis]*cs1
        hb2 = (2*bz[:,0]*bz[:,2])[:, np.newaxis]*cs2
        hb  = hb0 + hb1 + hb2  # => (NumCh, NfrqRsl)

        # 分母
        ha0 = (ap[:,0]**2 + ap[:,1]**2 + ap[:,2]**2)[:, np.newaxis]
        ha1 = (2*ap[:,1]*(ap[:,0]+ap[:,2]))[:, np.newaxis]*cs1
        ha2 = (2*ap[:,0]*ap[:,2])[:, np.newaxis]*cs2
        ha  = ha0 + ha1 + ha2

        H = np.sqrt(hb/ha)          # => (NumCh, NfrqRsl)
        Hnorm = H[:, -1]            # => (NumCh,)
        ACFFrsp *= H / Hnorm[:, np.newaxis]

    # 8) 计算 AsymFunc
    fd = freq_2d - Frs  # => (NumCh, NfrqRsl)

    # 强制 b, c => (NumCh, NfrqRsl)
    b = force_2d_shape(b, NumCh, NfrqRsl, do_interpolate=True)
    c = force_2d_shape(c, NumCh, NfrqRsl, do_interpolate=True)

    # be => (NumCh, NfrqRsl)
    be = b * ERBw  # 保证二者都 (NumCh, NfrqRsl)

    AsymFunc = np.exp(c * np.arctan2(fd, be))

    # freq_2d[0,:] 当作 1D 频率向量返回
    return ACFFrsp, freq_2d[0, :], AsymFunc
