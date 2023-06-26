import numpy as np
import geatpy as ea
from Initialization import Individual
import objectivefunction
import math
import scipy.stats
import matplotlib.pyplot as plt
from rank import rank
from VSM import Vsm
from DAEM import Daem
from Hv import Hv
from Mms import MMS

""" 改于 6.20 晚"""

n = 210  # 个体数  m=5 n=210 10-440   15-695
Maxiter = 300  # 最大迭代次数
m = 5  # 目标数

K = 2 * (m - 1)  # WFG
L = 20
n_min = 4  # 相同片段最少满足值
Dim = L + K   # 维度

mum = 5  # 突变参数
# OF = objectivefunction.MaF1(m, K, Dim)  # objective function
# P_pf = []  # The Optimal solution set
# population = np.asarray([Individual(_, Dim, m, n)  # 初始化个体
#                          for _ in range(n)])


Finall_mean_30times = []
Finall_std_30times = []
rank_Finall_mean_30times = []
rank_Finall_std_30times = []

MAXT = 5
of = [1, 2, 3, 4, 5, 6, 7, 8, 9]
# of =[1]
for o in of:
    if o == 1:
        Dim = L + K
        OF = objectivefunction.WFG1(m, K, Dim)
    elif o == 2:
        Dim = L + K
        OF = objectivefunction.WFG2(m, K, Dim)
    elif o == 3:
        Dim = L + K
        OF = objectivefunction.WFG3(m, K, Dim)
    elif o == 4:
        Dim = L + K
        OF = objectivefunction.WFG4(m, K, Dim)
    elif o == 5:
        Dim = L + K
        OF = objectivefunction.WFG5(m, K, Dim)
    elif o == 6:
        Dim = L + K
        OF = objectivefunction.WFG6(m, K, Dim)
    elif o == 7:
        Dim = L + K
        OF = objectivefunction.WFG7(m, K, Dim)

    elif o == 8:
        Dim = L + K
        OF = objectivefunction.WFG8(m, K , Dim)

    elif o == 9:
        Dim = L + K
        OF = objectivefunction.WFG9(m, K , Dim)

    Frut = []
    Frut1 = []
    for T in range(MAXT):  # 重复次数

        population = np.asarray([Individual(_, Dim, m, n)  # 初始化个体
                                 for _ in range(n)])
        P_pff = []  # 所有非控制解目标向量
        NoN_I = []  # 所有非控制解个体编码
        P_pf = []
        t = 0
        max_of = np.zeros(m)
        while t < Maxiter:
            Vac = []  # the set of vaccine
            ObjV = []  # 目标函数矩阵
            P_non = []
            al = (math.cos((t / Maxiter) * (math.pi / 2)) + 1) / 10  # alpha
            # Calculating the population fitness
            for i in range(n):
                population[i].fv = OF.fnc(population[i].code)
                ObjV.append(population[i].fv)
            ObjV = np.array(ObjV)
            for i in range(n):
                for j in range(m):
                    if ObjV[i][j] > max_of[j]:
                        max_of[j] = ObjV[i][j]

            # print(ObjV)
            # 基于ESS的非支配快速排序 levels 是第几层 levels =1 时 在PF上
            sort = [levels, criLevel] = ea.ndsortESS(ObjV=ObjV, needLevel=1)
            # print(sort[0])
            # Algorithm lin 7
            num_P_non = 0  # |Pnon|
            for i in range(len(sort[0])):
                if sort[0][i] == 1:
                    population[i].PF = 1

                    P_non.append(population[i])
                    num_P_non = num_P_non + 1
            # P_pf.append(P_non)
            n_s = math.ceil(num_P_non * al)
            n_v = math.ceil((Dim / m) * (Maxiter - t) / Maxiter + n_min)
            el = math.ceil(((Maxiter - t) / Maxiter * 0.3 + 0.1) * len(population))
            rank(P_non, m)
            # Rsum = []
            rr = {}  # 字典，key为sum，value为 编号
            for i in range(len(P_non)):
                P_non[i].Rsum = np.sum(P_non[i].rank)
                # Rsum.append(P_non[i].Rsum)
                P_non[i].Rpro = np.prod(P_non[i].rank)
                # print(P_non[i].rank, P_non[i].Rsum, P_non[i].Rpro)
                if P_non[i].Rsum in rr:
                    rr[P_non[i].Rsum].append(P_non[i].i)
                else:
                    rr[P_non[i].Rsum] = []
                    rr[P_non[i].Rsum].append(P_non[i].i)

            # 把当代非支配解放入非支配解集
            for i in range(len(P_non)):
                P_pff.append(P_non[i].fv.copy())
                NoN_I.append(P_non[i].code.copy())
            # np.savetxt('m={},OF={}个体编码3.txt'.format(m, o), NoN_I, fmt='%f')
            # np.savetxt('m={},OF={}个体目标向量3.txt'.format(m, o), P_pff, fmt='%f')

            Vac = Vsm(population, n_min, n_s, n_v, rr)  # 疫苗集
            Daem(population, el)
            MMS(population, mum, Vac)

            if t + 1 != Maxiter:
                for i in population:
                    # print(i.code)
                    i.PF = 0
                    i.ep = 0
                    i.Rsum = 0
                    i.Rpro = 0
                    i.fv = []

            t = t + 1
            # 第 t 代 进化结束

        # 对于每一代非支配解，提取最终的 非支配解集合
        # for i in range(len(ObjV)):
        #     P_pf.append(ObjV[i])
        # np.savetxt('m={},OF={}个体编码2.txt'.format(m, o), NoN_I, fmt='%f')
        # np.savetxt('m={},OF={}个体目标向量2.txt'.format(m, o), P_pff, fmt='%f')
        P_pff = np.array(P_pff, dtype=float)
        NoN_I = np.array(NoN_I, dtype=float)
        sort1 = [levels, criLevel] = ea.ndsortESS(ObjV=P_pff, needLevel=1)
        P_pf1 = P_pff.copy()
        NoN_I1 = NoN_I.copy()
        # np.savetxt('m={},OF={}个体编码1.txt'.format(m, o), NoN_I, fmt='%f')
        # np.savetxt('m={},OF={}个体目标向量1.txt'.format(m, o), P_pff, fmt='%f')
        NoN_I = []
        P_pff = []

        for i in range(len(sort1[0])):
            if sort1[0][i] == 1:
                P_pff.append(P_pf1[i])  # P_pff 是所有非支配解 的目标向量集合
                NoN_I.append(NoN_I1[i])  # NoN_I 是所有非支配解 的个体集合
        # np.savetxt('m={},OF={}个体编码0.txt'.format(m, o), NoN_I, fmt='%f')
        # np.savetxt('m={},OF={}个体目标向量0.txt'.format(m, o), P_pff, fmt='%f')
        P_pff = np.array(P_pff, dtype=float)
        NoN_I = np.array(NoN_I, dtype=float)
        # np.savetxt('m={},OF={}个体编码0.txt'.format(m, o), NoN_I, fmt='%f')
        # np.savetxt('m={},OF={}个体目标向量0.txt'.format(m, o), P_pff, fmt='%f')
        # print(np.shape(P_pf))
        print('P_pff----HV:')
        # print(len((P_pff)))
        rut1 = Hv(P_pff, max_of, o)
        # np.savetxt('m={},OF={}个体编码3.txt'.format(m, o), NoN_I, fmt='%f')
        # np.savetxt('m={},OF={}个体目标向量3.txt'.format(m, o), P_pff, fmt='%f')
        print(rut1)
        aaa = P_pff.copy()
        aa = aaa.T
        # np.savetxt('m={},OF={}个体编码3.txt'.format(m, o), NoN_I, fmt='%f')
        # np.savetxt('m={},OF={}个体目标向量3.txt'.format(m, o), P_pff, fmt='%f')

        P_rank = np.ones((len(P_pff), m))
        # print(np.shape(P_rank))
        # print(np.shape(aa))
        # print(P_pf)
        for mm in range(m):
            b = np.argsort(aa[mm])
            for i in range(len(b)):
                P_rank[b[i]][mm] = i + 1
        # print(P_rank)
        P_rank_sum = np.ones((len(P_pff), 1))
        P_rank_pro = np.ones((len(P_pff), 1))
        for i in range(len(P_pff)):
            P_rank_sum[i] = np.sum(P_rank[i])
            P_rank_pro[i] = np.prod(P_rank[i])

        max_num = len(P_pff)
        bb = np.reshape(P_rank_sum, -1)  # bb 每个个体的排名总和
        cc = np.reshape(P_rank_pro, -1)  # cc 每个个体的排名乘积
        b = np.argsort(bb)  # b排名索引
        # print(b)
        # c = bb.rank(method='dense', ascending=False)

        P_pfpf = []
        I_pfpf = []
        # print(b)
        rut = 0
        if max_num > n:
            Psum = bb[b[0]]
            P_pfpf.append(P_pff[b[0]])
            I_pfpf.append(NoN_I[b[0]])
            for i in range(1, max_num):
                if Psum == bb[b[i]]:
                    if cc[b[i - 1]] < cc[b[i]]:
                        P_pfpf[-1] = P_pff[b[i]]
                        I_pfpf[-1] = NoN_I[b[i]]

                    elif cc[b[i - 1]] > cc[b[i]]:  # 乘积小的不选
                        continue
                    else:  # 乘积一样的不选
                        continue
                        # P_pfpf.append(P_pff[b[i]])
                        # I_pfpf.append(NoN_I[b[i]])
                else:
                    Psum = bb[b[i]]
                    P_pfpf.append(P_pff[b[i]])
                    I_pfpf.append(NoN_I[b[i]])

                if len(P_pfpf) == n:
                    break

                # P_pfpf.append(P_pff[b[i]])
                # I_pfpf.append(NoN_I[b[i]])
            P_pfpf = np.array(P_pfpf)
            # print(np.shape(P_pfpf))
            print('P_pfpf----HV:')
            rut = Hv(P_pfpf, max_of, o)
            print(rut)
        # if rut>rut1:
        Frut.append(rut)  # 用rank的
        # else:
        Frut1.append(rut1)  # 没rank 的
        # P_pfpf = np.reshape(P_pfpf,(len(b),1))
        # print(P_pf,np.shape(P_pf),type(P_pf))

        # for i in range(len()):
        #     pass

        # print(P_rank_sum)
        # print(P_rank_pro)
        # print('done')
        np.savetxt('m={},OF={},popcode.txt'.format(m, o), NoN_I, fmt='%f')
        np.savetxt('m={},OF={},objVector.txt'.format(m, o), P_pff, fmt='%f')
        if len(I_pfpf) > 0:
            np.savetxt('m={},OF={}popcode(rank).txt'.format(m, o), I_pfpf, fmt='%f')
            np.savetxt('m={},OF={}objVector(rank).txt'.format(m, o), P_pfpf, fmt='%f')

    Frut = np.array(Frut)
    Frut1 = np.array(Frut1)
    # print("第{}个任务30次HV的平均值：{}，标准差：{}".format(o,np.mean(Frut),np.std(Frut)))
    # var = [np.mean(Frut),np.std(Frut)]
    print("(rank)massion:{},mean:{}std:{}".format(o, np.mean(Frut), np.std(Frut)))
    print("massion:{},mean:{}std:{}".format(o, np.mean(Frut1), np.std(Frut1)))
    # aaaaa = np.reshape(np.mean(Frut),-1)
    Finall_mean_30times.append(np.mean(Frut1))
    Finall_std_30times.append(np.std(Frut1))
    rank_Finall_mean_30times.append(np.mean(Frut))
    rank_Finall_std_30times.append(np.std(Frut))

np.savetxt('mean(rank),m={},OF={}.txt'.format(m, o), rank_Finall_mean_30times, fmt='%f')
np.savetxt('std(rank),m={},OF={}.txt'.format(m, o), rank_Finall_std_30times, fmt='%f')
np.savetxt('mean,m={},OF={}.txt'.format(m, o), Finall_mean_30times, fmt='%f')
np.savetxt('std,m={},OF={}.txt'.format(m, o), Finall_std_30times, fmt='%f')
