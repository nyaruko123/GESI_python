import numpy as np
from .Fr2Fpeak import Fr2Fpeak
from .GammaChirpFrsp import GammaChirpFrsp
from .AsymCmpFrspV2 import AsymCmpFrspV2

def CmprsGCFrsp(Fr1, fs=48000, n=4, b1=1.81, c1=-2.96, frat=1, b2=2.17, c2=2.20, NfrqRsl=1024):
    """
    Compute compressive GammaChirp frequency response, using only one line of Fp1 
    if Fr2Fpeak() returns two arrays.

    1) Detect whether Fr2Fpeak returns a single array or a tuple of two arrays.
    2) If it's a tuple (fp1_low, fp1_high), stack them by np.vstack => shape (2,NumCh).
       Then pick [0,:] to keep only the first row (NumCh,).
    3) If it's already a single array, proceed as usual.
    """

    # ---------------------- (可选)调试打印 ----------------------
    print("DEBUG in CmprsGCFrsp: shapes of input arguments:")
    print("  Fr1.shape =", np.shape(Fr1))
    print("  n.shape   =", np.shape(n))
    print("  b1.shape  =", np.shape(b1))
    print("  c1.shape  =", np.shape(c1))
    print("  frat.shape=", np.shape(frat))
    print("  b2.shape  =", np.shape(b2))
    print("  c2.shape  =", np.shape(c2))
    print(f"  fs={fs}, NfrqRsl={NfrqRsl}")
    print("---------------------------------------------------")

    # ---------------------- 参数扩展 ----------------------
    Fr1 = np.array(Fr1, ndmin=1)
    NumCh = len(Fr1)

    if np.isscalar(n):
        n = np.full(NumCh, n)
    if np.isscalar(b1):
        b1 = np.full(NumCh, b1)
    if np.isscalar(c1):
        c1 = np.full(NumCh, c1)
    if np.isscalar(frat):
        frat = np.full(NumCh, frat)
    if np.isscalar(b2):
        b2 = np.full(NumCh, b2)
    if np.isscalar(c2):
        c2 = np.full(NumCh, c2)

    # ---------------------- 被动GammaChirp频率响应 ----------------------
    pGCFrsp, freq, _, _, _ = GammaChirpFrsp(Fr1, fs, n, b1, c1, 0, NfrqRsl)

    # ---------------------- 获取 Fp1 ----------------------
    fp1_return = Fr2Fpeak(n, b1, c1, Fr1)
    # 如果是 tuple: e.g. (fp1_low, fp1_high) => 合并
    if isinstance(fp1_return, tuple):
        # 假设只有2个
        # fp1_return[0] shape => (NumCh,)
        # fp1_return[1] shape => (NumCh,)
        print("DEBUG: Fr2Fpeak returned a tuple. vstack => shape (2,NumCh)")
        Fp1_2d = np.vstack(fp1_return)  # => shape (2, NumCh)
    else:
        # 否则它可能是单个数组
        Fp1_2d = np.asarray(fp1_return)
        # 若它是一维(NumCh,) => shape
        # 若它是(2,NumCh) => 保持原样
        # 也可能是其他情况

    # 若 Fp1_2d 的第一维是 2，说明 2 行 => 只取[0,:]
    if Fp1_2d.ndim == 2 and Fp1_2d.shape[0] == 2:
        print("DEBUG: Detected Fp1_2d has shape (2, NumCh) => picking row[0,:] only.")
        Fp1 = Fp1_2d[0, :]  # shape => (NumCh,)
    else:
        # 如果是单行或(NumCh,)等，直接用
        Fp1 = Fp1_2d

    # ---------------------- 计算 Fr2 并传给AsymCmpFrspV2 ----------------------
    Fr2 = frat * Fp1  # => shape (NumCh,) or (?)
    ACFFrsp, freq, AsymFunc = AsymCmpFrspV2(Fr2, fs, b2, c2, NfrqRsl)

    # ---------------------- 合成压缩输出 ----------------------
    cGCFrsp = pGCFrsp * AsymFunc  # => shape (NumCh, NfrqRsl)

    # 计算峰值，用于归一化
    ValFp2 = np.max(cGCFrsp, axis=1)  # => (NumCh,)
    NormFactFp2 = 1.0 / ValFp2       # => (NumCh,)

    cGCresp = {
        'Fr1': Fr1,
        'n': n,
        'b1': b1,
        'c1': c1,
        'frat': frat,
        'b2': b2,
        'c2': c2,
        'NfrqRsl': NfrqRsl,
        'pGCFrsp': pGCFrsp,
        'cGCFrsp': cGCFrsp,
        'cGCNrmFrsp': cGCFrsp * (NormFactFp2[:, np.newaxis]),
        'ACFFrsp': ACFFrsp,
        'AsymFunc': AsymFunc,
        'Fp1': Fp1,
        'Fr2': Fr2,
        'Fp2': freq[np.argmax(cGCFrsp, axis=1)],
        'ValFp2': ValFp2,
        'NormFctFp2': NormFactFp2,
        'freq': freq
    }

    return cGCresp
