import numpy as np

def GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, CompressionHealth, PindB):
    """
    Compute asymmetric function output in dB scale.

    Parameters:
        GCparam: dict
            Parameters for the gammatone filter bank.
        GCresp: dict
            Response data of the filter bank.
        Fr1query: float
            Frequency query specified by Fr1, usually used in specifying filter bank frequency (not Fp1).
        CompressionHealth: float
            Compression health parameter.
        PindB: float
            Input level in dB.

    Returns:
        AFoutdB: float
            Output of the asymmetric function in dB.
        IOfuncdB: float
            Input-output function in dB.
        GCparam: dict
            Updated GCparam with the normalization parameter.
    """
    # 新增维度验证
    assert 'ERBw2' in GCresp, "ERBw2 not found in GCresp"
    assert GCresp['ERBw2'].shape == GCresp['Fp1'].shape, \
        f"ERBw2维度不匹配 {GCresp['ERBw2'].shape} vs Fp1 {GCresp['Fp1'].shape}"

    GCparam['AsymFunc_NormdB'] = 100

    # 增强类型检查
    if isinstance(PindB, (list, np.ndarray)):
        if PindB.size == 1:
            PindB = float(PindB)
        else:
            raise ValueError("PindB必须是标量值")
    elif not isinstance(PindB, (int, float)):
        raise TypeError("PindB必须是int/float类型")

    # 新增参数维度调试
    print(f"[DEBUG] ERBw2维度: {GCresp['ERBw2'].shape}")
    print(f"[DEBUG] Fp1维度: {GCresp['Fp1'].shape}")

    # 计算非对称函数输出
    AFoutLin = CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, PindB)
    AFoutLinNorm = CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, GCparam['AsymFunc_NormdB'])

    # 转换到dB标度
    AFoutdB = 20 * np.log10(AFoutLin / AFoutLinNorm)
    IOfuncdB = AFoutdB + PindB

    return AFoutdB, IOfuncdB, GCparam

def CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, PindB):
    """
    Compute the asymmetric function in linear scale.

    Parameters:
        GCparam: dict
            Parameters for the gammatone filter bank.
        GCresp: dict
            Response data of the filter bank.
        Fr1query: float
            Frequency query specified by Fr1.
        CompressionHealth: float
            Compression health parameter.
        PindB: float
            Input level in dB.

    Returns:
        AFoutLin: float
            Output of the asymmetric function in linear scale.
    """
    # 强制转换为NumPy数组
    Fr1_array = np.asarray(GCparam['Fr1'])
    assert Fr1_array.ndim == 1, "Fr1应为1维数组"

    # 新增维度统一检查
    n_channels = Fr1_array.size
    assert GCresp['PcHPAF'].shape[1] == n_channels, \
        f"PcHPAF维度不匹配 {GCresp['PcHPAF'].shape} vs 通道数 {n_channels}"

    # 查找最近频率索引
    nch = np.argmin(np.abs(Fr1_array - Fr1query))
    
    # 索引安全约束（关键修改）
    nch = np.clip(nch, 0, n_channels - 1).item()  # 确保转为Python整数
    
    # 调试输出（增强版）
    print(f"[DEBUG] 当前索引nch={nch} (允许范围0-{n_channels-1})")
    print(f"[DEBUG] PcHPAF有效索引范围: 0-{GCresp['PcHPAF'].shape[1]-1}")
    print(f"[DEBUG] frat0Pc有效索引范围: 0-{GCresp['frat0Pc'].shape[1]-1}")

    # 参数提取（带维度断言）
    assert GCresp['PcHPAF'].shape[0] == 2, "PcHPAF应为二维数组(2, N)"
    PcHPAF = GCresp['PcHPAF'][0, nch].item()  # 转换为Python标量

    assert GCresp['frat0Pc'].shape == (2, n_channels), "frat0Pc维度错误"
    frat0Pc = GCresp['frat0Pc'][0, nch].item()

    assert GCresp['frat1val'].ndim == 1, "frat1val应为1维数组"
    frat1val = GCresp['frat1val'][nch].item()

    # 输入验证
    if not isinstance(PindB, (int, float)):
        raise TypeError(f"PindB必须为标量，当前类型: {type(PindB)}")

    # 计算动态参数
    frat = frat0Pc + frat1val * (PindB - PcHPAF)
    Fr2 = frat * GCresp['Fp1'][0, nch].item()  # Fp1维度已验证

    # 获取ERB参数（关键修改）
    assert GCresp['ERBw2'].shape == (2, n_channels), "ERBw2维度错误"
    ERBw2 = GCresp['ERBw2'][0, nch].item()

    # 计算压缩参数
    assert GCresp['b2val'].shape == (2, n_channels), "b2val维度错误"
    b2E = GCresp['b2val'][0, nch].item() * ERBw2

    assert GCresp['c2val'].shape == (2, n_channels), "c2val维度错误"
    c2CH = CompressionHealth * GCresp['c2val'][0, nch].item()

    # 计算最终输出（添加异常捕获）
    try:
        AFoutLin = np.exp(c2CH * np.arctan2(
            GCresp['Fp1'][0, nch].item() - Fr2,
            b2E
        ))
    except FloatingPointError as e:
        print(f"数学计算错误: {e}")
        print(f"调试参数: Fr2={Fr2}, b2E={b2E}")
        AFoutLin = 1.0  # 安全回退值

    return AFoutLin
