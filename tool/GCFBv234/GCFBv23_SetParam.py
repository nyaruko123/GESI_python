import numpy as np
from .tool.EqualFreqScale import EqualFreqScale
from tool.GCFBv234.tool.Freq2ERB import Freq2ERB
from .tool.Fr2Fpeak import Fr2Fpeak
from .GCFBv23_HearingLoss import GCFBv23_HearingLoss

def GCFBv23_SetParam(GCparam):
    # ================== 参数初始化阶段 ==================
    # 基础参数设置
    GCparam.setdefault('fs', 48000)
    GCparam.setdefault('OutMidCrct', 'ELC')
    GCparam.setdefault('NumCh', 100)
    GCparam.setdefault('FRange', [100, 6000])
    
    # 参数校验增强
    if not isinstance(GCparam['NumCh'], int) or GCparam['NumCh'] < 1:
        raise ValueError("NumCh必须是正整数")
    if len(GCparam['FRange']) != 2 or GCparam['FRange'][0] >= GCparam['FRange'][1]:
        raise ValueError("FRange参数格式错误，应为[min_freq, max_freq]")
    
    # 采样率校验
    if GCparam['FRange'][1] * 3 > GCparam['fs']:
        print('警告: GCFB可能无法正常工作，因 max(FreqRange)*3 > fs')
        print(GCparam)
        input('按回车键继续...')
    
    # ================== Gammachirp参数配置 ==================
    GCparam.setdefault('n', 4)
    GCparam.setdefault('b1', [1.81, 0])
    GCparam.setdefault('c1', [-2.96, 0])
    GCparam.setdefault('frat', [[0.466, 0], [0.0109, 0]])
    GCparam.setdefault('b2', [[2.17, 0], [0, 0]])
    GCparam.setdefault('c2', [[2.20, 0], [0, 0]])
    GCparam.setdefault('Ctrl', 'dynamic')
    GCparam.setdefault('GainCmpnstdB', -1)
    
    # 控制模式处理
    if GCparam['Ctrl'].startswith('fix'):
        GCparam['Ctrl'] = 'static'
    if GCparam['Ctrl'].startswith('tim'):
        GCparam['Ctrl'] = 'dynamic'
    if not GCparam['Ctrl'].startswith(('sta', 'dyn', 'lev')):
        raise ValueError('控制模式错误: 应为"static", "dynamic"或"level-estimation"')
    
    # ================== 频率生成核心逻辑 ==================
    # 生成基础频率参数
    freq_result = EqualFreqScale('ERB', GCparam['NumCh'], GCparam['FRange'])
    if not (isinstance(freq_result, (tuple, list)) and len(freq_result) >= 2):
        raise ValueError("EqualFreqScale返回结构异常，应返回(Fr1, ERBrate1)")
    
    Fr1_raw, ERBrate1 = freq_result[0], freq_result[1]
    
    # 维度处理增强
    Fr1 = np.asarray(Fr1_raw, dtype=np.float32).squeeze()
    if Fr1.ndim > 1:
        Fr1 = np.mean(Fr1, axis=1) if Fr1.shape[1] == 2 else Fr1[:,0]
    Fr1 = Fr1[:GCparam['NumCh']]  # 强制截断
    
    # 异常处理
    if len(Fr1) != GCparam['NumCh']:
        print(f"[维度修正] Fr1长度自动调整 {len(Fr1)}→{GCparam['NumCh']}")
        Fr1 = np.linspace(GCparam['FRange'][0], GCparam['FRange'][1], GCparam['NumCh'])
        
    GCparam['Fr1'] = np.clip(Fr1, 50, 12000)
    GCresp = {'Fr1': GCparam['Fr1'], 'ERBspace1': np.mean(np.diff(ERBrate1))}
    
    # ================== ERB计算模块 ==================
    try:
        ERBrate, ERBw = Freq2ERB(Fr1)
        ERBrate1kHz, _ = Freq2ERB(1000)
        GCresp['Ef'] = ERBrate / ERBrate1kHz - 1
    except Exception as e:
        raise RuntimeError(f"Freq2ERB计算失败: {str(e)}")
    
    # 安全计算处理
    with np.errstate(invalid='ignore', divide='ignore'):
        Fr1_safe = Fr1.astype(np.float64)
        GCresp['ERBw2'] = np.tile(24.7 * (4.37 * (Fr1_safe / 1000.0) + 1.0), (2, 1))
        GCresp['ERBw2'] = np.nan_to_num(GCresp['ERBw2'], nan=24.7)
    
    # ================== 修复 KeyError: 'LvlEst' ==================
    GCparam.setdefault('LvlEst', {})
    GCparam['LvlEst'].setdefault('ExpDecayVal', 0.5)  # decay exponential 默认值
    GCparam['LvlEst'].setdefault('LctERB', 1.5)
    GCparam['LvlEst'].setdefault('DecayHL', 0.5)
    GCparam['LvlEst'].setdefault('b2', GCparam['b2'][0][0])
    GCparam['LvlEst'].setdefault('c2', GCparam['c2'][0][0])
    GCparam['LvlEst'].setdefault('frat', 1.08)
    GCparam['LvlEst'].setdefault('RMStoSPLdB', 30)
    GCparam['LvlEst'].setdefault('Weight', 0.5)
    GCparam['LvlEst'].setdefault('RefdB', 50)
    GCparam['LvlEst'].setdefault('Pwr', [1.5, 0.5])
    GCparam['LvlEst'].setdefault('LvlLinRef', 1)
    GCparam['LvlEst'].setdefault('LvlLinMinLim', 1e-8)
    
    # ------------------ 新增默认：为每个通道设定级别估计通道索引 ------------------
    GCparam['LvlEst'].setdefault('NchLvlEst', np.arange(GCparam['NumCh']))
    
    # ================== 修复 KeyError: 'DynHPAF' ==================
    GCparam.setdefault('DynHPAF', {})  # 确保 DynHPAF 存在
    
    # 设定默认的帧处理参数
    GCparam['DynHPAF'].setdefault('StrPrc', 'sample-by-sample')  # 默认逐样本处理
    GCparam['DynHPAF'].setdefault('Tframe', 0.001)  # 默认帧长 1ms
    GCparam['DynHPAF'].setdefault('Tshift', 0.0005)  # 默认帧移 0.5ms
    GCparam['DynHPAF'].setdefault('LenFrame', int(GCparam['DynHPAF']['Tframe'] * GCparam['fs']))  # 帧长度
    GCparam['DynHPAF'].setdefault('LenShift', int(GCparam['DynHPAF']['Tshift'] * GCparam['fs']))  # 帧移长度
    GCparam['DynHPAF'].setdefault('fs', 1 / GCparam['DynHPAF']['Tshift'])  # 帧移采样率
    GCparam['DynHPAF'].setdefault('NameWin', 'hanning')  # 窗函数类型
    GCparam['DynHPAF'].setdefault('ValWin', np.hanning(GCparam['DynHPAF']['LenFrame']))  # 默认窗口
    
    # 确保窗口值归一化，防止 NaN
    if GCparam['DynHPAF']['ValWin'].sum() != 0:
        GCparam['DynHPAF']['ValWin'] /= np.sum(GCparam['DynHPAF']['ValWin'])
    
    # ================== 修复 KeyError: 'HLoss' ==================
    GCparam.setdefault('HLoss', {})  # 确保 HLoss 这个字典存在
    # 修改：FB_CompressionHealth默认值设为NumCh个1的数组
    GCparam['HLoss'].setdefault('FB_CompressionHealth', np.ones(GCparam['NumCh']))
    
    # ================== 关键参数生成 ==================
    OneVec = np.ones(GCparam['NumCh'])
    GCresp['b1val'] = GCparam['b1'][0] * OneVec + GCparam['b1'][1] * GCresp['Ef']
    GCresp['c1val'] = GCparam['c1'][0] * OneVec + GCparam['c1'][1] * GCresp['Ef']
    
    # ================== 新增核心修复：b2val/c2val计算 ==================
    GCresp['b2val'] = np.array([
        GCparam['b2'][0][0] * OneVec + GCparam['b2'][0][1] * GCresp['Ef'],
        GCparam['b2'][1][0] * OneVec + GCparam['b2'][1][1] * GCresp['Ef']
    ])
    GCresp['c2val'] = np.array([
        GCparam['c2'][0][0] * OneVec + GCparam['c2'][0][1] * GCresp['Ef'],
        GCparam['c2'][1][0] * OneVec + GCparam['c2'][1][1] * GCresp['Ef']
    ])
    assert GCresp['b2val'].shape == (2, GCparam['NumCh']), f"b2val维度错误: {GCresp['b2val'].shape}"
    assert GCresp['c2val'].shape == (2, GCparam['NumCh']), f"c2val维度错误: {GCresp['c2val'].shape}"
    
    # ================== 动态参数计算修复部分 ==================
    GCresp['frat0val'] = GCparam['frat'][0][0] * OneVec + GCparam['frat'][0][1] * GCresp['Ef']
    GCresp['frat1val'] = GCparam['frat'][1][0] * OneVec + GCparam['frat'][1][1] * GCresp['Ef']
    
    # PcHPAF计算及维度扩展
    GCresp['PcHPAF'] = (1 - GCresp['frat0val']) / GCresp['frat1val']
    GCresp['PcHPAF'] = np.tile(GCresp['PcHPAF'], (2, 1))  # 维度扩展为(2, NumCh)
    
    # frat0Pc生成
    GCresp['frat0Pc'] = GCresp['frat0val'] + GCresp['frat1val'] * GCresp['PcHPAF']
    
    # ================== 维度对齐修复核心 ==================
    try:
        fp1_low, fp1_high = Fr2Fpeak(GCparam['n'], GCresp['b1val'], GCresp['c1val'], Fr1)
        GCresp['Fp1'] = np.vstack([fp1_low, fp1_high])
        assert GCresp['Fp1'].shape == (2, GCparam['NumCh']), f"Fp1维度异常: {GCresp['Fp1'].shape}"
    except Exception as e:
        raise RuntimeError(f"Fp1生成失败: {str(e)}")
    
    LvldB = GCparam.get('LeveldBscGCFB', 50)
    print(f"[维度校验] frat0Pc形状: {GCresp['frat0Pc'].shape}")
    print(f"[维度校验] frat1val形状: {GCresp['frat1val'].shape}")
    print(f"[维度校验] PcHPAF形状: {GCresp['PcHPAF'].shape}")
    
    # 计算 fratVal
    fratVal = GCresp['frat0Pc'] + GCresp['frat1val'] * (LvldB - GCresp['PcHPAF'])
    
    if fratVal.ndim == 1:
        fratVal = np.tile(fratVal, (2, 1))
    elif fratVal.shape[0] != 2:
        fratVal = fratVal.T
    
    assert fratVal.shape == GCresp['Fp1'].shape, f"维度不匹配: fratVal{fratVal.shape} vs Fp1{GCresp['Fp1'].shape}"
    
    # ================== 最终计算阶段 ==================
    Fr2_temp = fratVal * GCresp['Fp1']
    # 只取第一行
    Fp1_single = Fr2_temp[0, :]  # shape => (NumCh,)
    GCresp['Fr2'] = Fp1_single[:, np.newaxis]
    GCresp['Fr2'] = np.asarray(GCresp['Fr2'], dtype=np.float64)
    
    print(f"[维度追踪] fratVal形状: {fratVal.shape} → Fp1形状: {GCresp['Fp1'].shape}")
    print(f"[最终输出] Fr2维度: {GCresp['Fr2'].shape}")
    
    # ================== 后处理阶段 ==================
    GCparam = GCFBv23_HearingLoss(GCparam, GCresp)
    
    print(f"[系统状态] 成功生成 {GCparam['NumCh']} 通道滤波器")
    print(f"          频率范围: {Fr1.min():.1f}Hz - {Fr1.max():.1f}Hz")
    print(f"          ERB空间均值: {GCresp['ERBspace1']:.2f} ERB")
    
    # --- 新增：将 frat 值存入 GCresp，便于后续调用使用 ---
    GCresp['frat'] = GCparam['LvlEst']['frat']
    
    return GCparam, GCresp
