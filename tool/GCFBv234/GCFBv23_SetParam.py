import numpy as np
from .tool.EqualFreqScale import EqualFreqScale
from tool.GCFBv234.tool.Freq2ERB import Freq2ERB
from .tool.Fr2Fpeak import Fr2Fpeak
from .GCFBv23_HearingLoss import GCFBv23_HearingLoss

def GCFBv23_SetParam(GCparam):
    # Handling Input Parameters
    GCparam.setdefault('fs', 48000)
    GCparam.setdefault('OutMidCrct', 'ELC')
    GCparam.setdefault('NumCh', 100)
    GCparam.setdefault('FRange', [100, 6000])

    # ========== 新增参数校验 ==========
    if not isinstance(GCparam['NumCh'], int) or GCparam['NumCh'] < 1:
        raise ValueError("NumCh必须是正整数")
    if len(GCparam['FRange']) != 2 or GCparam['FRange'][0] >= GCparam['FRange'][1]:
        raise ValueError("FRange参数格式错误，应为[min_freq, max_freq]")

    if GCparam['FRange'][1] * 3 > GCparam['fs']:
        print(GCparam)
        print('GCFB may not work properly when max(FreqRange)*3 > fs.')
        print('---> Set fs properly. OR If you wish to continue as is, press RETURN > ')
        input()

    # Gammachirp parameters
    GCparam.setdefault('n', 4)
    GCparam.setdefault('b1', [1.81, 0])
    GCparam.setdefault('c1', [-2.96, 0])
    GCparam.setdefault('frat', [[0.466, 0], [0.0109, 0]])
    GCparam.setdefault('b2', [[2.17, 0], [0, 0]])
    GCparam.setdefault('c2', [[2.20, 0], [0, 0]])
    GCparam.setdefault('Ctrl', 'dynamic')
    GCparam.setdefault('GainCmpnstdB', -1)

    if GCparam['Ctrl'].startswith('fix'):
        GCparam['Ctrl'] = 'static'
    if GCparam['Ctrl'].startswith('tim'):
        GCparam['Ctrl'] = 'dynamic'

    if not (GCparam['Ctrl'].startswith('sta') or GCparam['Ctrl'].startswith('dyn') or GCparam['Ctrl'].startswith(
            'lev')):
        raise ValueError('Specify GCparam.Ctrl: "static", "dynamic", or "level(-estimation)"')

    # Parameters for level estimation
    if 'PpgcRef' in GCparam or 'LvlRefdB' in GCparam:
        print('The parameter "GCparam.PpgcRef" is obsolete.')
        print('The parameter "GCparam.LvlRefdB" is obsolete.')
        raise ValueError('Please change it to GCparam.GainRefdB.')

    GCparam.setdefault('GainRefdB', 'NormIOfunc')
    GCparam.setdefault('LeveldBscGCFB', 50)
    GCparam.setdefault('LvlEst', {})
    GCparam['LvlEst'].setdefault('LctERB', 1.5)
    GCparam['LvlEst'].setdefault('DecayHL', 0.5)
    GCparam['LvlEst'].setdefault('b2', GCparam['b2'][0][0])
    GCparam['LvlEst'].setdefault('c2', GCparam['c2'][0][0])
    GCparam['LvlEst'].setdefault('frat', 1.08)
    GCparam['LvlEst'].setdefault('RMStoSPLdB', 30)
    GCparam['LvlEst'].setdefault('Weight', 0.5)
    GCparam['LvlEst'].setdefault('RefdB', 50)
    GCparam['LvlEst'].setdefault('Pwr', [1.5, 0.5])
    GCparam.setdefault('NumUpdateAsymCmp', 1)

    # Sample-by-sample or Frame-base processing
    GCparam.setdefault('DynHPAF', {})
    GCparam['DynHPAF'].setdefault('StrPrc', 'sample-by-sample')

    if GCparam['DynHPAF']['StrPrc'].startswith('frame'):
        GCparam['DynHPAF']['Tframe'] = 0.001
        GCparam['DynHPAF']['Tshift'] = 0.0005
        GCparam['DynHPAF']['LenFrame'] = int(GCparam['DynHPAF']['Tframe'] * GCparam['fs'])
        GCparam['DynHPAF']['LenShift'] = int(GCparam['DynHPAF']['Tshift'] * GCparam['fs'])
        GCparam['DynHPAF']['Tframe'] = GCparam['DynHPAF']['LenFrame'] / GCparam['fs']
        GCparam['DynHPAF']['Tshift'] = GCparam['DynHPAF']['LenShift'] / GCparam['fs']
        GCparam['DynHPAF']['fs'] = 1 / GCparam['DynHPAF']['Tshift']
        GCparam['DynHPAF']['NameWin'] = 'hanning'
        GCparam['DynHPAF']['ValWin'] = np.hanning(GCparam['DynHPAF']['LenFrame'])
        GCparam['DynHPAF']['ValWin'] /= np.sum(GCparam['DynHPAF']['ValWin'])

    # ========== 核心修改：Fr1维度处理 ==========
    # Generate frequency parameters
    freq_result = EqualFreqScale('ERB', GCparam['NumCh'], GCparam['FRange'])
    
    # 处理不同版本的返回结构
    if isinstance(freq_result, (tuple, list)) and len(freq_result) >= 2:
        Fr1_raw, ERBrate1 = freq_result[0], freq_result[1]
        
        # 处理二维数组情况（假设返回低/高频率边界）
        if isinstance(Fr1_raw, np.ndarray) and Fr1_raw.ndim == 2:
            if Fr1_raw.shape[0] == GCparam['NumCh'] and Fr1_raw.shape[1] == 2:
                Fr1 = np.mean(Fr1_raw, axis=1)  # 取中心频率
            else:
                Fr1 = Fr1_raw.flatten()[:GCparam['NumCh']]  # 安全截断
        else:
            Fr1 = np.array(Fr1_raw).flatten()
    else:
        raise ValueError("EqualFreqScale返回结构异常，应返回(Fr1, ERBrate1)")
    
    # 强制类型转换和维度验证
    Fr1 = np.asarray(Fr1, dtype=np.float32).squeeze()  # 去除单维度
    if Fr1.ndim != 1:
        Fr1 = Fr1.reshape(-1)[:GCparam['NumCh']]  # 强制转为1维并截断
        
    # 最终维度校验
    if len(Fr1) != GCparam['NumCh']:
        print(f"[WARN] Fr1长度自动修正 {len(Fr1)}→{GCparam['NumCh']}")
        Fr1 = np.linspace(GCparam['FRange'][0], GCparam['FRange'][1], GCparam['NumCh'])  # 生成备用频率
        
    # 数值范围保护
    Fr1 = np.clip(Fr1, 50, 12000)
    
    GCparam['Fr1'] = Fr1
    GCresp = {}
    GCresp['Fr1'] = Fr1
    GCresp['ERBspace1'] = np.mean(np.diff(ERBrate1))
    
    # ========== 参数计算保护 ==========
    try:
        ERBrate, ERBw = Freq2ERB(Fr1)
        ERBrate1kHz, ERBw1kHz = Freq2ERB(1000)
        GCresp['Ef'] = ERBrate / ERBrate1kHz - 1
    except Exception as e:
        raise RuntimeError(f"Freq2ERB计算失败: {str(e)}")

    # 安全计算ERBw2
    with np.errstate(invalid='ignore', divide='ignore'):
        Fr1_safe = Fr1.astype(np.float64)
        GCresp['ERBw2'] = 24.7 * (4.37 * (Fr1_safe/1000.0) + 1.0)
        GCresp['ERBw2'] = np.nan_to_num(GCresp['ERBw2'], nan=24.7)

    # ========== 后续参数生成 ==========
    OneVec = np.ones(GCparam['NumCh'])
    GCresp['b1val'] = GCparam['b1'][0] * OneVec + GCparam['b1'][1] * GCresp['Ef']
    GCresp['c1val'] = GCparam['c1'][0] * OneVec + GCparam['c1'][1] * GCresp['Ef']
    
    # Fr2Fpeak维度保护
    try:
        fp1_low, fp1_high = Fr2Fpeak(GCparam['n'], GCresp['b1val'], GCresp['c1val'], Fr1)
    except ValueError as e:
        raise RuntimeError(f"Fr2Fpeak计算错误: {str(e)}，请检查输入参数")
    
    GCresp['Fp1'] = np.vstack([fp1_low, fp1_high])

    # 调试信息
    print(f"[GCFBv23] 成功生成 {GCparam['NumCh']} 通道滤波器")
    print(f"          Fr1范围: {Fr1.min():.1f}Hz - {Fr1.max():.1f}Hz")
    print(f"          ERB空间均值: {GCresp['ERBspace1']:.2f} ERB")

    return GCparam, GCresp
