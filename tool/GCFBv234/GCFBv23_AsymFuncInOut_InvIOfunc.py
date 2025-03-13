import numpy as np
from .GCFBv23_AsymFuncInOut import GCFBv23_AsymFuncInOut

def GCFBv23_AsymFuncInOut_InvIOfunc(GCparam, GCresp, Fr1query, CmprsHlthQuery, IOfuncdB):
    # 定义 PindB 值范围 (-120dB 到 150dB，步长 0.1)
    PindBList = np.arange(-120.0, 150.1, 0.1)  # 用 NumPy 数组更高效

    # 遍历 PindBList，分别计算 IOfuncdBlist
    IOfuncdBlist = np.array([
        GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, CmprsHlthQuery, float(PindB))[1]
        for PindB in PindBList
    ])

    # 插值计算 PindB
    PindB = np.interp(IOfuncdB, IOfuncdBlist, PindBList)

    return PindB
