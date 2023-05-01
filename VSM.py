import numpy as np
import math

def eucliDist(A,B):
    return math.sqrt(sum([(a - b)**2 for (a,b) in zip(A,B)]))



def Vsm(P,n_min,n_s,n_v,rr):
    V_set = []  # 保存的是被选个体的编号
    N_s = n_s
    ii=[]
    rrr = sorted(rr.keys())
    Vac =[]
    for i in rrr:
        if N_s > len(rr[i]): # 先选sum小的
            N_s = N_s - len(rr[i])
            for m in range(len(rr[i])):
                V_set.append(rr[i][m])
        else:
            for j in range(len(rr[i])): # 在同一 sum里选pro大的
                ii.append(rr[i][j])
            ii.sort(reverse=True)
            for k in range(N_s):
                V_set.append(rr[i][k])
            break
    eD = np.ones((len(V_set),len(P)))  # 每一个被选的点与其他所有点在目标空间的欧氏距离。
    eDR = [] # 每一个被选的点与其他所有点在目标空间欧氏距离 从小到大的排名、保存的的是个体序号

    for i in range(len(V_set)):   # 每一个被选的点与其他所有点在目标空间的欧氏距离。
        for j in range(len(P)):
            eD[i][j] = eucliDist(P[V_set[i]].fv,P[j].fv)
        eDR.append(np.argsort(eD[i]))

    # 每一组是 以eDR为准， eDR的每一行的编号输入在 V_set中得到备选个体的编号， eDR每一行中的向量代表 对应备选个体欧氏距离从大到小排名的序号。
    G =[]  # 疫苗小组集合
    for i in range(len(V_set)):
        G.append(eDR[i][:n_v])

    for i in range(len(V_set)):
        vac = {}  # 疫苗片段 是一个字典
        for d in range(P[0].D_multitask):
            vvac = []  # 中间变量用来提取疫苗片段
            SUM = 0
            for j in range(n_v):
                vvac.append(P[eDR[i][j]].Nd[d])
            maxnum = vvac.count(max(vvac, key=vvac.count))  # 出现最多元素的 个数
            maxTimes = max(vvac, key=vvac.count)   # 哪个元素最多

            if maxnum >= n_min:
                vac[d] = maxTimes

        Vac.append(vac)
    # Vac每一个元素是一个字典，每一个字典的key代表维度，value 代表 邻域，元素编号 从1到n_s（等于V_set个数）
    # V_set 每一个存的是 被选个体的编号。————> 所以疫苗片段是，P[V_set[i]].Nd[Vac[i][j](key的值是d)] = Vac[i][j](value)
    return Vac










