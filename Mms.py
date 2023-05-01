# Multimodal mutation strategy
import numpy as np
import math


def GM(P, mum):
    for i in range(len(P)):
        if np.random.uniform() < (1 / len(P)):
            u = np.random.uniform()
            if u <= 0.5:
                delta = (2 * u) ** (1 / (1 + mum)) - 1
                P[i] = P[i] + delta * P[i]
                if P[i] > 1:
                    P[i] = 1 - 1e-10
                if P[i] < 0:
                    P[i] = 0
            else:
                delta = 1 - (2 * (1 - u)) ** (1 / (1 + mum))
                P[i] = P[i] + delta * (1 - P[i])
                if P[i] > 1:
                    P[i] = 1 - 1e-10
                if P[i] < 0:
                    P[i] = 0


def MMS(P, mum, Vac):
    b = int(len(Vac) / 2)  # 超参
    for i in P:
        if i.PF == 1 and i.code != []:
            GM(i.code, mum)
            i.Nd = i.Ndd(i.code, i.n)

        if i.PF == 0 and i.code != []:
            u = np.random.randint(len(Vac))  # 随机一个 疫苗
            for key, value in Vac[u].items():
                i.code[key] = (1 / i.n) * value + (1 / i.n) * np.random.uniform()
                if i.code[key] > 1:
                    i.code[key] = 1 - 1e-10
                if i.code[key] < 0:
                    i.code[key] = 0
            i.Nd = i.Ndd(i.code, i.n)

        if i.code == []:
            i.code = np.reshape(np.random.uniform(0, 1, size=(1, i.D_multitask)), -1)
            for j in range(b):
                u = np.random.randint(len(Vac))
                for key, value in Vac[u].items():

                    i.code[key] = 0
                    i.code[key] = ((1 / i.n) * value + (1 / i.n) * np.random.uniform()) / b

            i.Nd = i.Ndd(i.code, i.n)
            i.age = []
