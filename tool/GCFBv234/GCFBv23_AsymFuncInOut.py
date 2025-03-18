import numpy as np

def GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, CompressionHealth, PindB):
    """
    Compute asymmetric function output in dB scale.
    
    如果 PindB 传入的是数组或列表，将取其第一个元素，并给出警告。
    """
    # 确保 PindB 为 NumPy 数组
    PindB = np.asarray(PindB)
    if PindB.size == 1:
        PindB = float(PindB)  # 转换为标量
    else:
        # 修改：当 PindB 不是标量时，给出警告并取第一个元素
        #print("Warning: PindB 不是标量，自动取第一个元素。")
        PindB = float(PindB.flat[0])
    
    # 维度验证
    assert 'ERBw2' in GCresp, "ERBw2 not found in GCresp"
    assert GCresp['ERBw2'].shape == GCresp['Fp1'].shape, \
        f"ERBw2 维度不匹配 {GCresp['ERBw2'].shape} vs Fp1 {GCresp['Fp1'].shape}"

    GCparam['AsymFunc_NormdB'] = 100

    # 计算非对称函数输出（线性尺度）
    AFoutLin = CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, PindB)
    AFoutLinNorm = CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, GCparam['AsymFunc_NormdB'])

    # 转换到 dB 标度
    AFoutdB = 20 * np.log10(AFoutLin / AFoutLinNorm)
    IOfuncdB = AFoutdB + PindB

    return AFoutdB, IOfuncdB, GCparam

def CalAsymFunc(GCparam, GCresp, Fr1query, CompressionHealth, PindB):
    """
    Compute the asymmetric function in linear scale.
    """
    # 强制转换为 NumPy 数组
    Fr1_array = np.asarray(GCparam['Fr1'])
    assert Fr1_array.ndim == 1, "Fr1 应为 1 维数组"

    # 维度检查
    n_channels = Fr1_array.size
    if 'PcHPAF' not in GCresp or GCresp['PcHPAF'] is None:
        raise ValueError("GCresp['PcHPAF'] 为空或未初始化，请检查 GCFBv23_SetParam 是否正确生成 GCresp")

    if GCresp['PcHPAF'].ndim < 2 or GCresp['PcHPAF'].shape[1] != n_channels:
        raise ValueError(f"PcHPAF 维度错误: {GCresp['PcHPAF'].shape}, 需要至少是 (2, N)")

    # 查找最近频率索引
    nch = np.argmin(np.abs(Fr1_array - Fr1query))
    nch = np.clip(nch, 0, n_channels - 1).item()  # 转换为 Python 整数

    # 提取参数
    assert GCresp['PcHPAF'].shape[0] == 2, "PcHPAF 应为二维数组(2, N)"
    PcHPAF = GCresp['PcHPAF'][0, nch].item()

    assert GCresp['frat0Pc'].shape == (2, n_channels), "frat0Pc 维度错误"
    frat0Pc = GCresp['frat0Pc'][0, nch].item()

    assert GCresp['frat1val'].ndim == 1, "frat1val 应为 1 维数组"
    frat1val = GCresp['frat1val'][nch].item()

    if not isinstance(PindB, (int, float)):
        raise TypeError(f"PindB 必须为标量，当前类型: {type(PindB)}")

    # 计算动态参数
    frat = frat0Pc + frat1val * (PindB - PcHPAF)
    Fr2 = frat * GCresp['Fp1'][0, nch].item()

    # 获取 ERB 参数
    assert GCresp['ERBw2'].shape == (2, n_channels), "ERBw2 维度错误"
    ERBw2 = GCresp['ERBw2'][0, nch].item()

    assert GCresp['b2val'].shape == (2, n_channels), "b2val 维度错误"
    b2E = GCresp['b2val'][0, nch].item() * ERBw2

    assert GCresp['c2val'].shape == (2, n_channels), "c2val 维度错误"
    c2CH = CompressionHealth * GCresp['c2val'][0, nch].item()

    # 计算最终输出
    try:
        AFoutLin = np.exp(c2CH * np.arctan2(
            GCresp['Fp1'][0, nch].item() - Fr2,
            b2E
        ))
    except FloatingPointError as e:
        print(f"数学计算错误: {e}")
        print(f"调试参数: Fr2={Fr2}, b2E={b2E}")
        AFoutLin = 1.0  # 设定安全回退值

    return AFoutLin
