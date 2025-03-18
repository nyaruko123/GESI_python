"""
INPUT:
    Metric:
    ParamSigmoid:
        Case STOI, ESTOI, GECI objective measures
            ParamSigmoid:
            length(ParamSigmoid) == 2 --- param for sigmoid [STOIparam.a, STOIparam.b]
            length(ParamSigmoid) == 3 --- [STOIparam.a, STOIparam.b, MaxPcorrect]
        Case HASPI or multiple metrics
        MaxPcorrect (option) was added to the end of the parameters 23 May 2022

OUTPUT:
    Pcorrect (%)

Note:
    STOI type
        ParamSigmoid = [STOIparam.a, STOIparam.b, MaxPcorrect]
        a = STOIparam.a;
        b = STOIparam.b;
        a = -13.1903; # Taal et al., ICASSP Proc., 2011
        b = 6.5293;
    HASPI_v1
        arg = bias + wgtcep * CepCorr + sum(wgtcov * cov3')
        Intel = 1.0 / (1.0 + exp(-arg)) # Logistic (logsig) function
"""

import numpy as np

def Metric2Pcorrect_Sigmoid(Metric, ParamSigmoid):
    # 如果 Metric 是标量，则转换为列表
    if np.isscalar(Metric):
        Metric = [Metric]
        
    LenMetric = len(Metric)
    LenPS = len(ParamSigmoid)

    if LenPS == LenMetric + 1:  # default
        MaxPcorrect = 100
    elif LenPS == LenMetric + 2:  # Setting MaxPcorrect
        MaxPcorrect = ParamSigmoid[LenPS - 1]
        ParamSigmoid = ParamSigmoid[:LenPS - 1]
    else:
        raise ValueError('Lengths of Metric and ParamSigmoid are inconsistent.')

    if LenMetric == 1:  # STOI, ESTOI type
        Pcorrect = MaxPcorrect / (1 + np.exp(ParamSigmoid[0] * Metric[0] + ParamSigmoid[1]))
    else:  # HASPI type
        ValArg = np.dot(ParamSigmoid, np.concatenate(([1], Metric)))
        Pcorrect = MaxPcorrect / (1 + np.exp(-ValArg))

    return Pcorrect
