import numpy as np

def getHist(array, bins=40, density=True, half=False):
    a0 = array
    if half:
        a0 = []
        for a in array:
            a0.append(a)
            a0.append(-a)
    h, edges = np.histogram(a0, bins=bins, density=density)
    l = len(h)
    edg0 = []
    for i in range(l):
        edg0.append(0.5*(edges[i]+edges[i+1]))
    return h, np.array(edg0)
