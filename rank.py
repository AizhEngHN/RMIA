import numpy as np

def rank(P_non,M):
    aa = np.ones((M,len(P_non)))
    for m in range(M):
        for i in range(len(P_non)):
            aa[m][i] = P_non[i].fv[m]

    for m in range(M):
        b = np.argsort(aa[m])
        for i in range(len(b)):
            P_non[b[i]].rank[m] = i + 1

    # return aa

