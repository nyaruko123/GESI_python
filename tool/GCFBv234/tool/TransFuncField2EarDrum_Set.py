import numpy as np


def TransFuncField2EarDrum_Set(StrCrct, FreqList=None):
    """
    Transfer function from field to ear drum various set.

    Parameters:
        StrCrct: Type of transfer function (e.g., 'FreeField', 'DiffuseField', 'ITU').
        FreqList: Optional list of frequencies to select from the table.

    Returns:
        FreqTbl: Frequency table.
        FrspdBTbl: Corresponding response values.
        Param: Dictionary containing the type of field to ear drum.
    """
    Param = {}

    if StrCrct.startswith('FreeField') or StrCrct.upper() == 'FF':
        Param['TypeField2EarDrum'] = 'FreeField'
        FreqTbl, FrspdBTbl = TransFuncFreeField2EarDrum_Moore16()
    elif StrCrct.startswith('DiffuseField') or StrCrct.upper() == 'DF':
        Param['TypeField2EarDrum'] = 'DiffuseField'
        FreqTbl, FrspdBTbl = TransFuncDiffuseField2EarDrum_Moore16()
    elif StrCrct.startswith('ITU'):
        Param['TypeField2EarDrum'] = 'ITU'
        FreqTbl, FrspdBTbl = TransFuncField2EarDrum_ITU()
    else:
        raise ValueError("Specify: FreeField (FF) / DiffuseField (DF) / ITU")

    if FreqList is None:
        return FreqTbl, FrspdBTbl, Param

    # Selection of FreqList
    for freq in FreqList:
        NumFreq = np.where(FreqTbl == freq)[0]
        if NumFreq.size > 0:
            FreqTbl = np.append(FreqTbl, FreqTbl[NumFreq])
            FrspdBTbl = np.append(FrspdBTbl, FrspdBTbl[NumFreq])
        else:
            raise ValueError(f"Freq {freq} is not listed on the table.")

    return FreqTbl, FrspdBTbl, Param


def TransFuncFreeField2EarDrum_Moore16():
    """
    Transfer function from Free Field to Ear Drum.

    Returns:
        FreqTbl: Frequency table.
        FrspdBTbl: Corresponding response values.
    """
    table = np.array([
        [20.0, 0.0],
        [25.0, 0.0],
        [31.5, 0.0],
        [40.0, 0.0],
        [50.0, 0.0],
        [63.0, 0.0],
        [80.0, 0.0],
        [100.0, 0.0],
        [125.0, 0.1],
        [160.0, 0.3],
        [200.0, 0.5],
        [250.0, 0.9],
        [315.0, 1.4],
        [400.0, 1.6],
        [500.0, 1.7],
        [630.0, 2.5],
        [750.0, 2.7],
        [800.0, 2.6],
        [1000.0, 2.6],
        [1250.0, 3.2],
        [1500.0, 5.2],
        [1600.0, 6.6],
        [2000.0, 12.0],
        [2500.0, 16.8],
        [3000.0, 15.3],
        [3150.0, 15.2],
        [4000.0, 14.2],
        [5000.0, 10.7],
        [6000.0, 7.1],
        [6300.0, 6.4],
        [8000.0, 1.8],
        [9000.0, -0.9],
        [10000.0, -1.6],
        [11200.0, 1.9],
        [12500.0, 4.9],
        [14000.0, 2.0],
        [15000.0, -2.0],
        [16000.0, 2.5],
    ])

    FreqTbl = table[:, 0]
    FrspdBTbl = table[:, 1]
    return FreqTbl, FrspdBTbl


def TransFuncDiffuseField2EarDrum_Moore16():
    """
    Transfer function from Diffuse Field to Ear Drum.

    Returns:
        FreqTbl: Frequency table.
        FrspdBTbl: Corresponding response values.
    """
    table = np.array([
        [20.0, 0.0],
        [25.0, 0.0],
        [31.5, 0.0],
        [40.0, 0.0],
        [50.0, 0.0],
        [63.0, 0.0],
        [80.0, 0.0],
        [100.0, 0.0],
        [125.0, 0.1],
        [160.0, 0.3],
        [200.0, 0.4],
        [250.0, 0.5],
        [315.0, 1.0],
        [400.0, 1.6],
        [500.0, 1.7],
        [630.0, 2.2],
        [750.0, 2.7],
        [800.0, 2.9],
        [1000.0, 3.8],
        [1250.0, 5.3],
        [1500.0, 6.8],
        [1600.0, 7.2],
        [2000.0, 10.2],
        [2500.0, 14.9],
        [3000.0, 14.5],
        [3150.0, 14.4],
        [4000.0, 12.7],
        [5000.0, 10.8],
        [6000.0, 8.9],
        [6300.0, 8.7],
        [8000.0, 8.5],
        [9000.0, 6.2],
        [10000.0, 5.0],
        [11200.0, 4.5],
        [12500.0, 4.0],
        [14000.0, 3.3],
        [15000.0, 2.6],
        [16000.0, 2.0],
    ])

    FreqTbl = table[:, 0]
    FrspdBTbl = table[:, 1]
    return FreqTbl, FrspdBTbl


def TransFuncField2EarDrum_ITU():
    """
    ITU transfer function from field to ear drum.

    Returns:
        FreqTbl: Frequency table.
        FrspdBTbl: Corresponding response values.
    """
    FreqTbl = np.array([
        0, 100, 125, 160, 200, 250, 315, 400, 500, 630,
        800, 1000, 1250, 1600, 2000, 2500, 3150, 4000,
        5000, 6300, 8000, 10000
    ])

    FrspdBTbl = np.array([
        0.0, 0.0, 0.0, 0.0, 0.3, 0.2, 0.5, 0.6, 0.7, 1.1,
        1.7, 2.6, 4.2, 6.5, 9.4, 10.3, 6.6, 3.2, 3.3, 16.0, 14.4
    ])

    return FreqTbl, FrspdBTbl
