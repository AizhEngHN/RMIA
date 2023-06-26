import numpy as np
import math


class MaF1():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.D = D

    def fnc(self, var):
        g = 0
        f = []

        for i in range(self.m - 1, self.D):  # 从0开始
            g = g + (var[i] - 0.5) * (var[i] - 0.5)

        f1 = (1 + g) * (1 - np.prod(var[:self.m - 1]))

        f.append(f1)

        for j in range(1, self.m - 1):
            f11 = (1 + g) * (1 - np.prod(var[:self.m - j - 1]) * (1 - var[self.m - 1 - j]))
            f.append(f11)
        fM = (1 + g) * var[0]
        f.append(fM)
        return f


class MaF2():  #
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.D = D

    def fnc(self, var):
        var_o = var.copy()
        var_c = var.copy()
        var_s = var.copy()
        g = []
        f = []
        for i in range(self.m):
            var_o[i] = math.pi / 2 * (var[i] / 2 + 0.25)
            var_c[i] = math.cos(math.pi / 2 * (var[i] / 2 + 0.25))
            var_s[i] = math.sin(math.pi / 2 * (var[i] / 2 + 0.25))

        for i in range(1, self.m + 1):
            g1 = 0
            if i == self.m:
                up = self.D
            else:
                up = self.m + (i) * math.floor((self.D - self.m + 1) / self.m) - 1
            down = self.m + (i - 1) * math.floor((self.D - self.m + 1) / self.m)
            # xx = up - down + 1
            for j in range(down, up + 1):
                # if j > len(Xm):
                #     j = len(Xm)
                g1 = g1 + ((var[j - 1] / 2 + 0.25) - 0.5) ** 2
            g.append(g1)

        f1 = np.prod(var_c[:self.m - 1]) * (1 + g[0])
        f.append(f1)

        for i in range(1, self.m - 1):
            f11 = (1 + g[i]) * (np.prod(var_c[:self.m - i - 1]) * var_s[self.m - 1 - i])
            f.append(f11)
        fM = var_s[0] * (1 + g[self.m - 1])
        f.append(fM)
        return f


class MaF3():  #
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.D = D

    def fnc(self, var):
        g = 0
        f = []
        var_c = var.copy()
        var_s = var.copy()

        for i in range(self.m):
            var_c[i] = math.cos(math.pi / 2 * (var[i]))
            var_s[i] = math.sin(math.pi / 2 * (var[i]))

        for i in range(self.m - 1, self.D):  # 从0开始
            g = g + (var[i] - 0.5) * (var[i] - 0.5) - math.cos(20 * math.pi * (var[i] - 0.5))
        g = 100 * (self.K + g)

        f1 = (np.prod(var_c[:self.m - 1]) * (1 + g)) ** 4
        f.append(f1)

        for i in range(1, self.m - 1):
            f11 = ((1 + g) * (np.prod(var_c[:self.m - i - 1]) * var_s[self.m - 1 - i])) ** 4
            f.append(f11)

        fM = (var_s[0] * (1 + g)) ** 2
        f.append(fM)
        return f


class MaF4():  #
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.D = D

    def fnc(self, var):
        a = 2
        g = 0
        f = []
        var_c = var.copy()
        var_s = var.copy()

        for i in range(self.m):
            var_c[i] = math.cos(math.pi / 2 * (var[i]))
            var_s[i] = math.sin(math.pi / 2 * (var[i]))

        for i in range(self.m - 1, self.D):  # 从0开始
            g = g + (var[i] - 0.5) * (var[i] - 0.5) - math.cos(20 * math.pi * (var[i] - 0.5))
        g = 100 * (self.K + g)

        f1 = a * (1 - np.prod(var_c[:self.m - 1])) * (1 + g)
        f.append(f1)

        for i in range(1, self.m - 1):
            f11 = a ** (i + 1) * ((1 + g) * (1 - np.prod(var_c[:self.m - i - 1]) * var_s[self.m - 1 - i]))
            f.append(f11)

        fM = a ** (self.m) * (1 - var_s[0])* (1 + g)
        f.append(fM)
        return f


class MaF5():  #  没有4次方 暂时与PlatEMO 相同
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.D = D

    def fnc(self, var):
        a = 2
        b = 100
        g = 0
        f = []
        var_c = var.copy()
        var_s = var.copy()

        for i in range(self.m):
            var_c[i] = math.cos(math.pi / 2 * (var[i]**b))
            var_s[i] = math.sin(math.pi / 2 * (var[i]**b))

        for i in range(self.m - 1, self.D):  # 从0开始
            g = g + (var[i] - 0.5) * (var[i] - 0.5)

        f1 = (a ** (self.m)) * (((np.prod(var_c[:self.m - 1])) * (1 + g)))
        f.append(f1)

        for i in range(1, self.m - 1):
            f11 = a ** (self.m - i) * (((1 + g) * (np.prod(var_c[:self.m - i - 1]) * var_s[self.m - 1 - i])))
            f.append(f11)

        fM = a  * ((var_s[0])* (1 + g))
        f.append(fM)
        return f


class MaF6():  #
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.D = D

    def fnc(self, var):
        g = 0
        f=[]
        var_o = var.copy()
        var_c = var.copy()
        var_s = var.copy()
        for i in range(self.m - 1, self.D):  # 从0开始
            g = g + (var[i] - 0.5) * (var[i] - 0.5)

        for i in range(self.m):
            if i ==0:
                var_o[i] = (math.pi / 2) * var[i]
            else:
                var_o[i] = (math.pi / 2) * (1/(2*(1+g)))*(1+2*(g*var[i]))   # 于原论文不符 ，但是于PlatEMO 中的MF6一致。
                # var_o[i] =  (math.pi / 2) * var[i]
            var_c[i] = math.cos(var_o[i])
            var_s[i] = math.sin(var_o[i])

        f1 = ((np.prod(var_c[:self.m - 1])) * (1 + 100*g))
        f.append(f1)

        for i in range(1, self.m - 1):
            f11 =((1 + 100*g) * (np.prod(var_c[:self.m - i - 1]) * (var_s[self.m - i - 1])))
            f.append(f11)

        fM = (var_s[0])* (1 + 100*g)
        f.append(fM)
        return f

class MaF7():  #
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.D = D

    def fnc(self, var):
        g = 0
        f=[]
        h=0
        for i in range(self.m - 1, self.D):  # 从0开始
            g = g + var[i]
        g = 1+ (9/self.K)*g

        for i in range(self.m -1):
            f1 = var[i]
            f.append(f1)
            h = h+(f1/(1+g))*(1+math.sin(3*math.pi*f1))
        h = self.m - h

        fM = h * (1+g)
        f.append(fM)
        return f


class WFG1():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D



    def fnc(self, Vars):  # 目标函数
        # N, Lind = Vars.shape
        Vars = np.reshape(Vars, (1, self.Lind))

        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        D = 1
        L = Lind - K
        S = np.array([list(range(2, 2 * M + 1, 2))])

        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub

        A = np.ones((1, M - 1))
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, 2 * Lind + 1, 2)])
        t1 = np.zeros((N, K + L))
        t1[:, :K] = Z[:, :K]
        t1[:, K:] = s_linear(Z[:, K:], 0.35)
        t2 = np.zeros((N, K + L))
        t2[:, :K] = t1[:, :K]
        t2[:, K:] = b_flat(t1[:, K:], 0.8, 0.75, 0.85)
        t3 = b_poly(t2, 0.02)
        t4 = np.zeros((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))
        for i in range(1, M):
            t4[:, i - 1] = r_sum(
                t3[:,
                list(range(
                    (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))],
                list(
                    range(2 * ((i - 1) * K_divide_M_sub_1 + 1),
                          # 原论文里 +1 在 乘2的括号外  2 * (i - 1) * K_divide_M_sub_1 + 1  ， Plat复现错误
                          2 * i * K_divide_M_sub_1 + 1,
                          2)))
        t4[:, M - 1] = r_sum(t3[:, K:K + L],
                             list(range(2 * (K + 1), 2 * (K + L) + 1, 2)))
        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t4[:, [M - 1]], A[:, [i - 1]]],
                                   0) * (t4[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t4[:, [M - 1]]
        h = convex(x)
        h[:, [M - 1]] = mixed(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f,-1)
        return f




class WFG2():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D

    def fnc(self, Vars):  # 目标函数
        Vars = np.reshape(Vars, (1, self.Lind))
        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        # L = self.L
        S = np.array([list(range(2, 2 * M + 1, 2))])
        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub
        D = 1
        L = Lind - K
        A = np.ones((1, M - 1))
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, Lind * 2 + 1, 2)])
        t1 = np.zeros((N, K + L))
        t1[:, :K] = Z[:, :K]
        t1[:, K:] = s_linear(Z[:, K:], 0.35)
        t2 = np.zeros((N, int(K + L / 2)))
        t2[:, :K] = t1[:, :K]
        t2[:,
        K:int(K
              + L / 2)] = (t1[:, K::2] + t1[:, K + 1::2]
                           + 2 * np.abs(t1[:, K::2] - t1[:, K + 1::2])) / 3
        t3 = np.ones((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))
        for i in range(1, M):
            t3[:, i - 1] = r_sum(
                t2[:,
                list(range(
                    (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))],
                np.ones((1, K_divide_M_sub_1)))
        t3[:, M - 1] = r_sum(t2[:, K:int(K + L / 2)], np.ones((1, int(L / 2))))
        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t3[:, [M - 1]], A[:, [i - 1]]],
                                   0) * (t3[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t3[:, [M - 1]]
        h = convex(x)
        h[:, [M - 1]] = disc(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f, -1)
        return f


class WFG3():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D

    def fnc(self, Vars):  # 目标函数
        Vars = np.reshape(Vars, (1, self.Lind))
        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        S = np.array([list(range(2, 2 * M + 1, 2))])
        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub
        D = 1
        L = Lind - K
        A = [1] + [0] * (M - 2)
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, Lind * 2 + 1, 2)])
        t1 = np.zeros((N, K + L))
        t1[:, :K] = Z[:, :K]
        t1[:, K:] = s_linear(Z[:, K:], 0.35)
        t2 = np.zeros((N, int(K + L / 2)))
        t2[:, :K] = t1[:, :K]
        t2[:,
        K:int(K
              + L / 2)] = (t1[:, K::2] + t1[:, K + 1::2]
                           + 2 * np.abs(t1[:, K::2] - t1[:, K + 1::2])) / 3
        t3 = np.ones((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))
        for i in range(1, M):
            t3[:, i - 1] = r_sum(
                t2[:,
                list(range(
                    (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))],
                np.ones((1, K_divide_M_sub_1)))
        t3[:, M - 1] = r_sum(t2[:, K:int(K + L / 2)], np.ones((1, int(L / 2))))
        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t3[:, [M - 1]], A[:, [i - 1]]],
                                   0) * (t3[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t3[:, [M - 1]]
        h = linear(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f, -1)
        return f


class WFG4():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D

    def fnc(self, Vars):  # 目标函数
        Vars = np.reshape(Vars, (1, self.Lind))
        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        S = np.array([list(range(2, 2 * M + 1, 2))])
        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub
        D = 1
        L = Lind - K
        A = np.ones((1, M - 1))
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, Lind * 2 + 1, 2)])
        t1 = s_multi(Z, 30, 10, 0.35)
        t2 = np.zeros((N, int(K + L / 2)))
        t2[:, :K] = t1[:, :K]
        t2[:,
        K:int(K
              + L / 2)] = (t1[:, K::2] + t1[:, K + 1::2]
                           + 2 * np.abs(t1[:, K::2] - t1[:, K + 1::2])) / 3
        t2 = np.ones((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))
        for i in range(1, M):
            t2[:, i - 1] = r_sum(
                t1[:,
                list(range(
                    (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))],
                np.ones((1, K_divide_M_sub_1)))
        t2[:, M - 1] = r_sum(t1[:, K:K + L], np.ones((1, L)))
        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t2[:, [M - 1]], A[:, [i - 1]]],
                                   0) * (t2[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t2[:, [M - 1]]
        h = concave(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f, -1)
        return f


class WFG5():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D



    def fnc(self, Vars):  # 目标函数
        Vars = np.reshape(Vars, (1, self.Lind))
        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        S = np.array([list(range(2, 2 * M + 1, 2))])
        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub
        D = 1
        L = Lind - K
        A = np.ones((1, M - 1))
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, Lind * 2 + 1, 2)])
        t1 = s_decept(Z, 0.35, 0.001, 0.05)
        t2 = np.zeros((N, int(K + L / 2)))
        t2[:, :K] = t1[:, :K]
        t2[:,
        K:int(K
              + L / 2)] = (t1[:, K::2] + t1[:, K + 1::2]
                           + 2 * np.abs(t1[:, K::2] - t1[:, K + 1::2])) / 3
        t2 = np.ones((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))
        for i in range(1, M):
            t2[:, i - 1] = r_sum(
                t1[:,
                list(range(
                    (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))],
                np.ones((1, K_divide_M_sub_1)))
        t2[:, M - 1] = r_sum(t1[:, K:K + L], np.ones((1, L)))
        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t2[:, [M - 1]], A[:, [i - 1]]],
                                   0) * (t2[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t2[:, [M - 1]]
        h = concave(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f, -1)
        return f




class WFG6():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D

    def fnc(self, Vars):  # 目标函数
        Vars = np.reshape(Vars, (1, self.Lind))
        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        S = np.array([list(range(2, 2 * M + 1, 2))])
        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub
        D = 1
        L = Lind - K
        A = np.ones((1, M - 1))
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, Lind * 2 + 1, 2)])

        t1 = np.zeros((N, K + L))
        t1[:, :K] = Z[:, :K]
        t1[:, K:] = s_linear(Z[:, K:], 0.35)
        t2 = np.zeros((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))
        for i in range(1, M):
            t2[:, i - 1] = r_nonsep(t1[:,
                                    list(range(
                                        (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))], K_divide_M_sub_1)

        SUM = np.zeros((N, 1))
        for i in range(K + 1, K + L):
            for j in range(i + 1, K + L + 1):
                SUM = SUM + np.abs(t1[:, i - 1] - t1[:, j - 1])

        t2[:, M - 1] = (np.sum(t1[:, K:]) + SUM * 2) / np.ceil(L / 2) / (1 + 2 * L - 2 * np.ceil(L / 2))

        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t1[:, [M - 1]], A[:, [i - 1]]], 0) * (t2[:, [i - 1]] - 0.5) + 0.5

        x[:, [M - 1]] = t2[:, [M - 1]]

        h = concave(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f, -1)
        return f




class WFG7():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D

    def fnc(self, Vars):  # 目标函数
        Vars = np.reshape(Vars, (1, self.Lind))
        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        S = np.array([list(range(2, 2 * M + 1, 2))])
        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub
        D = 1
        L = Lind - K
        A = np.ones((1, M - 1))
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, Lind * 2 + 1, 2)])
        t1 = np.zeros((N, K + L))
        Y = (np.fliplr(np.cumsum(np.fliplr(Z)) - np.zeros((N, D))) - Z) / np.array([range(K + L - 1, -1, -1)])
        t1[:, : K] = np.power(Z[:, : K],
                              (0.02 + (50 - 0.02) * (0.98 / 49.98 - (1 - 2 * Y[:, :K]) * np.abs(np.floor(0.5 - Y[
                                                                                                               :,
                                                                                                               :K]) + 0.98 / 49.98))))
        # Z[:, 1: K]^ (0.02 + (50 - 0.02) * (0.98 / 49.98 - (1 - 2 * Y[:, 1:K])* np.abs(np.floor(0.5 - Y[
        #                           :, 1:K])+0.98 / 49.98)))
        t1[:, K:] = Z[:, K:]

        t2 = np.zeros((N, K + L))
        t2[:, :K] = t1[:, :K]
        t2[:, K:] = s_linear(t1[:, K:], 0.35)
        t3 = np.zeros((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))
        for i in range(1, M):
            t3[:, i - 1] = r_sum(
                t2[:,
                list(range(
                    (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))],
                np.ones((1, K_divide_M_sub_1)))
        t3[:, M - 1] = r_sum(t2[:, K:K + L], np.ones((1, L)))

        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t3[:, [M - 1]], A[:, [i - 1]]],
                                   0) * (t3[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t3[:, [M - 1]]
        h = concave(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f, -1)
        return f



class WFG8():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D

    def fnc(self, Vars):  # 目标函数
        Vars = np.reshape(Vars, (1, self.Lind))
        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        S = np.array([list(range(2, 2 * M + 1, 2))])
        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub
        D = 1
        L = Lind - K
        A = np.ones((1, M - 1))
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, Lind * 2 + 1, 2)])
        t1 = np.zeros((N, K + L))
        t1[:, : K] = Z[:, : K]
        Y = (np.cumsum(Z) - np.zeros((N, D)) - Z) / np.array([range(0, K + L, 1)])
        t1[:, K: K + L] = np.power(Z[:, K: K + L], (
                    0.02 + (50 - 0.02) * (0.98 / 49.98 - (1 - 2 * Y[:, K: K + L]) * np.abs(np.floor(0.5 - Y[
                                                                                                          :,
                                                                                                          K: K + L]) + 0.98 / 49.98))))
        # Z[:, 1: K]^ (0.02 + (50 - 0.02) * (0.98 / 49.98 - (1 - 2 * Y[:, 1:K])* np.abs(np.floor(0.5 - Y[
        #                           :, 1:K])+0.98 / 49.98)))

        t2 = np.zeros((N, K + L))
        t2[:, :K] = t1[:, :K]
        t2[:, K:] = s_linear(t1[:, K:], 0.35)
        t3 = np.zeros((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))
        for i in range(1, M):
            t3[:, i - 1] = r_sum(
                t2[:,
                list(range(
                    (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))],
                np.ones((1, K_divide_M_sub_1)))
        t3[:, M - 1] = r_sum(t2[:, K:K + L], np.ones((1, L)))

        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t3[:, [M - 1]], A[:, [i - 1]]],
                                   0) * (t3[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t3[:, [M - 1]]
        h = concave(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f, -1)
        return f




class WFG9():
    def __init__(self, m, K, D):
        self.m = m
        self.K = K
        self.Lind = D

    def fnc(self, Vars):  # 目标函数
        Vars = np.reshape(Vars, (1, self.Lind))
        N = 1
        Lind = self.Lind
        M = self.m
        K = self.K
        S = np.array([list(range(2, 2 * M + 1, 2))])
        ub = list(range(2, 2 * Lind + 1, 2))
        ub = np.array(ub)
        Vars = Vars * ub
        D = 1
        L = Lind - K
        A = np.ones((1, M - 1))
        A = np.tile(A, (N, 1))
        Z = Vars / np.array([range(2, Lind * 2 + 1, 2)])
        t1 = np.zeros((N, K + L))
        Y = (np.fliplr(np.cumsum(np.fliplr(Z)) - np.zeros((N, D))) - Z) / np.array([range(K + L - 1, -1, -1)])

        t1[:, : K + L - 1] = np.power(Z[:, : K + L - 1], (
                    0.02 + (50 - 0.02) * (0.98 / 49.98 - (1 - 2 * Y[:, : K + L - 1]) * np.abs(np.floor(0.5 - Y[
                                                                                                             :,
                                                                                                             : K + L - 1]) + 0.98 / 49.98))))

        t1[:, -1] = Z[:, -1]
        t2 = np.zeros((N, K + L))
        t2[:, :K] = s_decept(t1[:, :K], 0.35, 0.001, 0.05)
        t2[:, K:] = s_multi(t1[:, K:], 30, 95, 0.35)
        t3 = np.zeros((N, M))
        K_divide_M_sub_1 = int(K / (M - 1))

        for i in range(1, M):
            t3[:, i - 1] = r_nonsep(t2[:,
                                    list(range(
                                        (i - 1) * K_divide_M_sub_1, i * K_divide_M_sub_1))], K_divide_M_sub_1)


        SUM = np.zeros((N, 1))
        for i in range(K + 1, K + L):  # WFG6检查
            for j in range(i + 1, K + L + 1):
                SUM = SUM + np.abs(t2[:, i - 1] - t2[:, j - 1])



        t3[:, M - 1] = (np.sum(t2[:, K:]) + SUM * 2) / np.ceil(L / 2) / (1 + 2 * L - 2 * np.ceil(L / 2))

        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t3[:, [M - 1]], A[:, [i - 1]]],
                                   0) * (t3[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t3[:, [M - 1]]
        h = concave(x)
        f = D * x[:, [M - 1]] + S * h
        f = np.reshape(f, -1)
        return f

def b_param(y,Y,A,B,C):
    return y^(B+(C-B)*(A-(1-2*Y)* np.abs(np.floor(0.5-Y)+A)))

def s_multi(x, A, B, C):
    return (1 + np.cos((4 * A + 2) * np.pi * (0.5 - np.abs(x - C) / 2 /
                                              (np.floor(C - x) + C))) + 4 * B *
            (np.abs(x - C) / 2 / (np.floor(C - x) + C))**2) / (B + 2)

def disc(x):
    return 1 - x[:, [0]] * (np.cos(5 * np.pi * x[:, [0]]))**2


def concave(x):
    return np.fliplr(
        np.cumprod(
            np.hstack(
                [np.ones((x.shape[0], 1)), np.sin(x[:, :-1] * np.pi / 2)]),
            1)) * np.hstack([
                np.ones((x.shape[0], 1)),
                np.cos(
                    x[:, list(range(x.shape[1] - 1 - 1, -1, -1))] * np.pi / 2)
            ])



def convex(x):
    return np.fliplr(
        np.cumprod(
            np.hstack(
                [np.ones((x.shape[0], 1)), 1 - np.cos(x[:, :-1] * np.pi / 2)]),
            1)) * np.hstack([
        np.ones((x.shape[0], 1)),
        1 - np.sin(
            x[:, list(range(x.shape[1] - 1 - 1, -1, -1))] * np.pi / 2)
    ])


def mixed(x):
    return 1 - x[:,
               [0]] - np.cos(10 * np.pi * x[:, [0]] + np.pi / 2) / 10 / np.pi


def s_linear(x, A):
    return np.abs(x - A) / np.abs(np.floor(A - x) + A)


def b_flat(x, A, B, C):
    Output = A + np.min([0 * np.floor(x - B), np.floor(x - B)], 0) * A * (
            B - x) / B - np.min([0 * np.floor(C - x), np.floor(C - x)],
                                0) * (1 - A) * (x - C) / (1 - C)
    return np.round(Output, 6)


def b_poly(x, a):
    return np.sign(x) * np.abs(x) ** a


def r_sum(x, w):
    Output = np.sum(x * w, 1) / np.sum(w)
    return Output


def linear(x):
    return np.fliplr(
        np.cumprod(np.hstack([np.ones(
            (x.shape[0], 1)), x[:, :-1]]), 1)) * np.hstack([
                np.ones((x.shape[0], 1)),
                1 - x[:, list(range(x.shape[1] - 1 - 1, -1, -1))]
            ])



def s_decept(x, A, B, C):
    return 1 + (np.abs(x - A) - B) * (np.floor(x - A + B) * (1 - C +
                                                             (A - B) / B) /
                                      (A - B) + np.floor(A + B - x) *
                                      (1 - C +
                                       (1 - A - B) / B) / (1 - A - B) + 1 / B)


def r_nonsep(y,A):
    Output = np.zeros((np.shape(y)[0],1))

    for j in range(1,np.shape(y)[1]+1):
        # print(np.shape(y)[1]+1)
        Temp = np.zeros((np.shape(y)[0],1))
        # print(y)
        for k in range(A-1):
            # print(y[:, j])
            # print(((j + k)%(np.shape(y)[1])))

            Temp = Temp + np.abs(y[:, j-1]-y[:, 1+((j + k)%(np.shape(y)[1])) -1])
        Output = Output + y[:, j-1]+Temp
    Output = Output/ (np.shape(y)[1] / A) / np.ceil(A / 2) / (1 + 2 * A - 2 * np.ceil(A / 2))
    return Output
