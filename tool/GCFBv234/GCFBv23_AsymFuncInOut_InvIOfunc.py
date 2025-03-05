import numpy as np
from .GCFBv23_AsymFuncInOut import GCFBv23_AsymFuncInOut

def GCFBv23_AsymFuncInOut_InvIOfunc(GCparam, GCresp, Fr1query, CmprsHlthQuery, IOfuncdB):
    # Define a wide range of PindB values
    PindBList = [x * 0.1 for x in range(-1200, 1501)]  # -120 to 150 with step 0.1

    # Call the undefined function GCFBv23_AsymFuncInOut
    _, IOfuncdBlist = GCFBv23_AsymFuncInOut(GCparam, GCresp, Fr1query, CmprsHlthQuery, PindBList)

    # Interpolate to find PindB
    PindB = np.interp(IOfuncdB, IOfuncdBlist, PindBList)

    return PindB




