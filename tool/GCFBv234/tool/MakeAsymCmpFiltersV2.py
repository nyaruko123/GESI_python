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
        ERBw = np.asarray(ERBw, dtype=np.float64).squeeze()
        if ERBw.size == 0:
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
            arr = arr.ravel()
            if arr.shape[0] > target_shape[0]:
                arr = arr[:target_shape[0]]
            elif arr.shape[0] < target_shape[0]:
                arr = np.pad(arr, (0, target_shape[0] - arr.shape[0]), mode='edge')
            return np.broadcast_to(arr, target_shape)

    target_shape = Frs.shape
    b = _safe_broadcast(b, target_shape)
    c = _safe_broadcast(c, target_shape)
    ERBw = _safe_broadcast(ERBw, target_shape)
    
    # ------------------ 计算 ACFcoef ------------------
    NumCh = len(Frs)
    NumFilt = 4
    p0 = 2
    p1 = 1.7818 * (1 - 0.0791 * b) * (1 - 0.1655 * np.abs(c))
    p2 = 0.5689 * (1 - 0.1620 * b) * (1 - 0.0857 * np.abs(c))
    p3 = 0.2523 * (1 - 0.0244 * b) * (1 + 0.0574 * np.abs(c))
    p4 = 1.0724

    ACFcoef = {'fs': fs, 'bz': np.zeros((NumCh, 3, NumFilt)), 'ap': np.zeros((NumCh, 3, NumFilt))}
    
    for Nfilt in range(NumFilt):
        r = np.exp(-p1 * (p0 / p4) ** Nfilt * 2 * np.pi * b * ERBw / fs)
        delFrs = (p0 * p4) ** Nfilt * p2 * c * b * ERBw
        phi = 2 * np.pi * np.maximum(Frs + delFrs, 0) / fs
        psy = 2 * np.pi * np.maximum(Frs - delFrs, 0) / fs
        
        ap = np.column_stack((np.ones_like(r), -2 * r * np.cos(phi), r ** 2))
        bz = np.column_stack((np.ones_like(r), -2 * r * np.cos(psy), r ** 2))
        
        vwr = np.exp(1j * 2 * np.pi * Frs / fs)
        vwrs = np.column_stack((np.ones_like(vwr), vwr, vwr ** 2))
        nrm = np.abs(np.sum(vwrs * ap, axis=1) / np.sum(vwrs * bz, axis=1))
        bz = bz * nrm[:, np.newaxis]
        
        ACFcoef['ap'][:, :, Nfilt] = ap
        ACFcoef['bz'][:, :, Nfilt] = bz

    # ------------------ 计算 ACFcoefConv ------------------
    ACFcoefConv = {'ap': np.zeros((NumCh, 9)), 'bz': np.zeros((NumCh, 9))}
    for nch in range(NumCh):
        ap1 = [1]
        bz1 = [1]
        for Nfilt in range(NumFilt):
            ap1 = np.convolve(ACFcoef['ap'][nch, :, Nfilt], ap1)
            bz1 = np.convolve(ACFcoef['bz'][nch, :, Nfilt], bz1)
        
        ACFcoefConv['ap'][nch, :9] = ap1[:9]
        ACFcoefConv['bz'][nch, :9] = bz1[:9]
    
    return ACFcoef, ACFcoefConv
