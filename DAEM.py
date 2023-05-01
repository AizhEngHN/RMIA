# Dynamic age-based elimination mechanism
import numpy as np
import math

def Daem(Pop,el):
    a = 2/3
    b = 1/3
    sum = 0
    for i in Pop:
        if i.PF ==1:
            i.age.append(a)
        else:
            i.age.append(b)

    for i in Pop:
        i.ep = (1 - np.prod(i.age[:]))
        sum = sum + i.ep

    p=[]
    for i in range(el):
        p.append(np.random.uniform())
    p.sort()
    m = 0
    fenzi = 0

    for i in range(len(Pop)):
        fenzi1 = fenzi +Pop[i].ep
        if fenzi/sum <=p[m] <fenzi1/sum:
            Pop[i].code = []
            Pop[i].age = []
            m = m + 1
        if m >= el:
            break

        fenzi = fenzi1






