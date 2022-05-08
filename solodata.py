import numpy as np
import os


def readSoloArc(arcFileName):
    path = os.path.expanduser(arcFileName)
    raw = [line.rstrip() for line in open(path)]
    vol0 = []
    m0 = []
    x0 = []
    v0 = []
    for line in raw:
        l = line.split()
        if l[0] == 'L':
            vol0.append(float(l[1]))
        elif l[0] == '1':
            m0.append(float(l[2]))
            x0.append(float(l[4]))
            v0.append(float(l[6]))
    return np.array(x0), np.array(v0), np.array(m0), np.array(vol0)


def readSoloVbar(vbarFileName):
    path = os.path.expanduser(vbarFileName)
    raw = [line.rstrip() for line in open(path)]
    m0 = []
    vbar0 = []
    for line in raw:
        l = line.split()
        m0.append(float(l[1]))
        vbar0.append(float(l[3]))
    return np.array(vbar0), np.array(m0)
