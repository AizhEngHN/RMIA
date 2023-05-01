import math

import numpy as np
import scipy.optimize as op


def Neigh_divi(code, n):

    each_n = 1 / n  # 每个分区的值
    Nd = []
    for i in range(len(code)):
        for j in range(n):
            if j * each_n <= code[i] < (j + 1) * each_n:
                Nd.append(j)
    return Nd


class Individual(object):  # 生成个体

    def __init__(self, i, D_multitask, M, n):
        self.i = i
        self.D_multitask = D_multitask  # 维度
        self.PF = 0  # 1 代表在 PF中
        self.rank = np.reshape(np.ones((1, M)), -1)
        self.Rsum = 0
        self.Rpro = 0
        self.fv = []  # 目标函数适应度值
        self.age = []
        self.ep = 0
        self.n = 10  # self.n:分区数量    有n个分区 ， n越大 分区越多 每个分区越小   1 / self.n 是一个分区的大小
        # self.a = np.random.uniform(-1,1,size=(1,D_multitask))
        self.code = np.reshape(np.random.uniform(0, 1, size=(1, D_multitask)), -1)
        self.Nd = Neigh_divi(self.code, self.n)  # 邻域划分 n:分区数量

    def Ndd(self,code, n):
        each_n = 1 / n
        Nd = []
        for i in range(len(code)):
            for j in range(n):
                if j * each_n <= code[i] < (j + 1) * each_n:
                    Nd.append(j)
        return Nd
