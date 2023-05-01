import geatpy as ea
import numpy as np

def Hv(Objv,max_of,o):
    Max = max_of*1.1   # 每一列最大值
    if o ==2:
    # print(np.shape(max_of))
        Min = np.min(Objv,0)  # 每一列最小值
    else:
        Min = np.zeros(len(Objv[0]))
    # print(np.shape(Min))
    b = np.shape(Objv)

    # print(b)
    for i in range(b[0]):  # 归一化  b[0] 为行数
        for j in range(b[1]):  # b[1] 为列数
            Objv[i][j] = (Objv[i][j] - Min[j])/(Max[j] - Min[j])
    bbbb = np.ones((b[0], b[1]))
    # print(np.shape(Objv), np.shape(bbbb))
    # Objv = np.array(Objv)
    # print(type(Objv), type(bbbb))
    # print(Objv)
    HV = ea.indicator.HV(Objv,bbbb)
    # print(np.shape(Objv))
    # print(HV)
    return  HV