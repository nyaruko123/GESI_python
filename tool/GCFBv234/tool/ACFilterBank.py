import numpy as np

def ACFilterBank(ACFcoef, ACFstatus, SigIn, SwOrdr=0):
    if ACFstatus is None or len(ACFstatus) == 0:
        NumCh, Lbz, NumFilt = ACFcoef['bz'].shape
        NumCh, Lap, NumFilt = ACFcoef['ap'].shape
        if Lbz != 3 or Lap != 3:
            raise ValueError('No guarantee for usual IIR filters except for AsymCmpFilter.')

        ACFstatus = {
            'NumCh': NumCh,
            'NumFilt': NumFilt,
            'Lbz': Lbz,
            'Lap': Lap,
            'SigInPrev': np.zeros((NumCh, Lbz)),
            'SigOutPrev': np.zeros((NumCh, Lap, NumFilt)),
            'Count': 0
        }
        print('ACFilterBank: Initialization of ACFstatus')
        return np.array([]), ACFstatus

    NumChSig, LenSig = SigIn.shape
    if LenSig != 1:
        raise ValueError('Input Signal should be NumCh*1 vector (1 sample time-slice)')
    if NumChSig != ACFstatus['NumCh']:
        raise ValueError(f'NumChSig ({NumChSig}) != ACFstatus.NumCh ({ACFstatus["NumCh"]})')

    if 'verbose' in ACFcoef and ACFcoef['verbose'] == 1:
        Tdisp = 50  # ms
        Tcnt = ACFstatus['Count'] / (ACFcoef['fs'] // 1000)  # ms
        if ACFstatus['Count'] == 0:
            print('ACFilterBank: Start processing')
        elif Tcnt % Tdisp == 0:
            print(f'ACFilterBank: Processed {int(Tcnt)} (ms).')

    ACFstatus['Count'] += 1

    ACFstatus['SigInPrev'] = np.hstack((ACFstatus['SigInPrev'][:, 1:], SigIn))

    x = ACFstatus['SigInPrev']
    NfiltList = range(ACFstatus['NumFilt'])
    if SwOrdr == 1:
        NfiltList = reversed(NfiltList)

    for Nfilt in NfiltList:
        forward = ACFcoef['bz'][:, ACFstatus['Lbz']-1::-1, Nfilt] * x
        feedback = ACFcoef['ap'][:, ACFstatus['Lap']-1:1:-1, Nfilt] * ACFstatus['SigOutPrev'][:, 1:, Nfilt]

        fwdSum = np.sum(forward, axis=1)
        fbkSum = np.sum(feedback, axis=1)

        y = (fwdSum - fbkSum) / ACFcoef['ap'][:, 0, Nfilt]
        ACFstatus['SigOutPrev'][:, :, Nfilt] = np.hstack((ACFstatus['SigOutPrev'][:, 1:, Nfilt], y[:, np.newaxis]))
        x = ACFstatus['SigOutPrev'][:, :, Nfilt]

    SigOut = y
    return SigOut, ACFstatus

# 未定义的函数或变量: 无
