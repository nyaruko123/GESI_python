import numpy as np

def TransFuncMiddleEar_Moore16(FreqList=None):
    table = [
        (20.0, -39.6), (25.0, -32.0), (31.5, -25.85), (40.0, -21.4),
        (50.0, -18.5), (63.0, -15.9), (80.0, -14.1), (100.0, -12.4),
        (125.0, -11.0), (160.0, -9.6), (200.0, -8.3), (250.0, -7.4),
        (315.0, -6.2), (400.0, -4.8), (500.0, -3.8), (630.0, -3.3),
        (750.0, -2.9), (800.0, -2.6), (1000.0, -2.6), (1250.0, -4.5),
        (1500.0, -5.4), (1600.0, -6.1), (2000.0, -8.5), (2500.0, -10.4),
        (3000.0, -7.3), (3150.0, -7.0), (4000.0, -6.6), (5000.0, -7.0),
        (6000.0, -9.2), (6300.0, -10.2), (8000.0, -12.2), (9000.0, -10.8),
        (10000.0, -10.1), (11200.0, -12.7), (12500.0, -15.0), (14000.0, -18.2),
        (15000.0, -23.8), (16000.0, -32.3), (18000.0, -45.5), (20000.0, -50.0)
    ]

    if FreqList is None:
        FreqTbl = [pair[0] for pair in table]
        FrspdBTbl = [pair[1] for pair in table]
        return FreqTbl, FrspdBTbl

    #  **确保 FreqList 是列表**
    FreqList = np.atleast_1d(FreqList)  # **转换为 1D 数组，防止整数输入**
    
    FreqTbl = []
    FrspdBTbl = []
    for freq in FreqList:
        matches = [pair for pair in table if pair[0] == freq]
        if matches:
            FreqTbl.append(matches[0][0])
            FrspdBTbl.append(matches[0][1])
        else:
            raise ValueError(f'Freq {freq} is not listed on the table.')

    return FreqTbl, FrspdBTbl
