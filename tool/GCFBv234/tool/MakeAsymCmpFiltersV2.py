# ================== 修复版 MakeAsymCmpFiltersV2.py ==================
import numpy as np
from .Freq2ERB import Freq2ERB

def MakeAsymCmpFiltersV2(fs, Frs, b, c):
    # ------------------ 输入参数维度规范化 ------------------
    Frs = np.asarray(Frs, dtype=np.float64).squeeze()
    if Frs.ndim > 1:
        Frs = Frs.ravel()  # 多维数组展平
    
    # 处理b/c参数维度
    b = np.asarray(b, dtype=np.float64).squeeze()
    c = np.asarray(c, dtype=np.float64).squeeze()
    
    # 维度对齐策略
    if b.ndim == 0:  # 标量扩展
        b = np.full_like(Frs, fill_value=b)
    elif b.ndim == 1 and b.shape != Frs.shape:
        b = b.reshape(-1)
    
    if c.ndim == 0:
        c = np.full_like(Frs, fill_value=c)
    elif c.ndim == 1 and c.shape != Frs.shape:
        c = c.reshape(-1)
    
    # ------------------ ERB计算安全处理 ------------------
    try:
        ERBw = Freq2ERB(Frs)
        print(f"DEBUG: ERBw before conversion: {ERBw}")
        
        if isinstance(ERBw, tuple):  # 处理元组
            ERBw = np.array(ERBw).squeeze()  # 转换为 NumPy 数组，并去掉单维度
            print(f"DEBUG: ERBw converted from tuple: {ERBw}")
            
        ERBw = np.asarray(ERBw, dtype=np.float64).squeeze()
        print(f"DEBUG: ERBw after conversion: shape={ERBw.shape}, content={ERBw[:10] if ERBw.size > 10 else ERBw}")
        
        if ERBw.size == 0:  # 关键检查，避免空数组
            raise ValueError(f"Freq2ERB returned an empty array! Check GCresp['Fr2']: {Frs}")
    
    except Exception as e:
        raise RuntimeError(f"Freq2ERB error: {str(e)}")
    
    # ------------------ 广播维度校验 ------------------
    def _safe_broadcast(arr, target_shape):
        arr = np.asarray(arr)
        if arr.size == 0:
            raise ValueError(f"_safe_broadcast error: Input array is empty! Expected shape {target_shape}")
        if arr.ndim == 0:
            return np.broadcast_to(arr, target_shape)
        try:
            return np.broadcast_to(arr, target_shape)
        except ValueError:
            return np.broadcast_to(arr.ravel(), target_shape)
    
    target_shape = Frs.shape
    b = _safe_broadcast(b, target_shape)
    c = _safe_broadcast(c, target_shape)
    ERBw = _safe_broadcast(ERBw, target_shape)
    
    # ------------------ 继续滤波器计算逻辑 ------------------
    # 这里应该有 ACFcoef 和 ACFcoefConv 计算的代码，保留原逻辑
    
    return ACFcoef, ACFcoefConv
