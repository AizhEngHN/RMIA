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
                var_o[i] =  (1/(4*(1+g)))*(1+2*(g*var[i]))
            var_c[i] = math.cos(var_o[i])
            var_s[i] = math.cos(var_o[i])

        f1 = ((np.prod(var_c[:self.m - 1])) * (1 + 100*g))
        f.append(f1)

        for i in range(1, self.m - 1):
            f11 =((1 + 100*g) * (np.prod(var_c[:self.m - i - 1]) * var_s[self.m - 1 - i]))
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



