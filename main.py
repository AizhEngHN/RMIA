#coding:utf-8
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

n = 210  # 个体数
Maxiter = 300
m = 5  # 目标数

K = 10  # Maf1-MaF6 K =10 ， MaF7 K=20
n_min = 4
Dim = m + K - 1

mum = 5
OF = objectivefunction.MaF1(m, K, Dim)  # objective function
# P_pf = []  # The Optimal solution set

MAXT = 8
of = [1,2,3,4,5,6,7]
# of = [2]
for o in of:
    if o == 1:
        OF = objectivefunction.MaF1(m, K, Dim)
    elif o == 2:
        OF = objectivefunction.MaF2(m, K, Dim)
    elif o == 3:
        OF = objectivefunction.MaF3(m, K, Dim)
    elif o == 4:
        OF = objectivefunction.MaF4(m, K, Dim)
    elif o == 5:
        OF = objectivefunction.MaF5(m, K, Dim)
    elif o == 6:
        OF = objectivefunction.MaF6(m, K, Dim)
    elif o == 7:
        OF = objectivefunction.MaF7(m, K + 10, Dim)
    Frut = []
    Frut1 = []
    for T in range(MAXT):
        population = np.asarray([Individual(_, Dim, m, n)  # 初始化个体
                                 for _ in range(n)])
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

            Vac = Vsm(population, n_min, n_s, n_v, rr)  # 疫苗集
            Daem(population, el)
            MMS(population, mum, Vac)

            for i in population:
                # print(i.code)
                i.PF = 0
                i.ep = 0
                i.Rsum = 0
                i.Rpro = 0
                i.fv = []

            t = t + 1

            for i in range(len(ObjV)):
                P_pf.append(ObjV[i])
            P_pf = np.array(P_pf)
            sort1 = [levels, criLevel] = ea.ndsortESS(ObjV=P_pf, needLevel=1)
            P_pf1 = P_pf.copy()
            P_pf = []
            for i in range(len(sort1[0])):
                if sort1[0][i] == 1:
                    P_pf.append(P_pf1[i])

        P_pf = np.array(P_pf)
        # print(np.shape(P_pf))
        print('P_pf----HV:')
        rut1 = Hv(P_pf,max_of,o)
        print(rut1)
        aaa = P_pf.copy()
        aa = aaa.T

        P_rank = np.ones((len(P_pf), m))

        for mm in range(m):
            b = np.argsort(aa[mm])
            for i in range(len(b)):
                P_rank[b[i]][mm] = i + 1
        P_rank_sum = np.ones((len(P_pf), 1))
        P_rank_pro = np.ones((len(P_pf), 1))
        for i in range(len(P_pf)):
            P_rank_sum[i] = np.sum(P_rank[i])
            P_rank_pro[i] = np.prod(P_rank[i])

        max_num = len(P_pf)
        bb = np.reshape(P_rank_sum, -1)
        b = np.argsort(bb)
        P_pfpf = []
        rut = 0
        if max_num > n:
            for i in range(n):
                P_pfpf.append(P_pf[b[i]])
            P_pfpf = np.array(P_pfpf)
            # print(np.shape(P_pfpf))
            print('P_pfpf----HV:')
            rut = Hv(P_pfpf, max_of,o)
            print(rut)
        
        Frut.append(rut)
        
        Frut1.append(rut1)

    Frut = np.array(Frut)
    if len(Frut) > 0:
        print("rank_massion:{},mean:{}std:{}".format(o, np.mean(Frut), np.std(Frut)))
    else:
        print('done')
    
    print("massion:{},mean:{}std:{}".format(o, np.mean(Frut1), np.std(Frut1)))
    np.savetxt('m={},OF={}.txt'.format(m, o), Frut, fmt='%f')
