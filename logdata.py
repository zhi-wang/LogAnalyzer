from math import pi, cos, sin, sqrt
from matplotlib.pyplot import step
import numpy as np
import os


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


class Constant:
    avogadro = 6.02214076e23
    gasconst = 1.987204259e-3


class StaticMethod:
    def volume(a, b, c, alpha, beta, gamma):
        radian = 180.0 / pi
        al, be, ga = alpha/radian, beta/radian, gamma/radian
        ca, cb, cg, sg = cos(al), cos(be), cos(ga), sin(ga)
        vol = a*b*c*sqrt(sg**2 - ca**2 - cb**2 + 2.0*ca*cb*cg)
        return vol

    def density(mass, volume):
        avogadro = Constant.avogadro
        density = (mass/volume)*(1.e24/avogadro)
        return density

    def temperature(ekin, dof):
        gasconst = Constant.gasconst
        temper = 2.0*ekin/(dof*gasconst)
        return temper


class GLogData:
    def __init__(self, mass, dof, moldof, file, dir):
        path = os.path.join(os.path.expanduser(dir), file)
        self.rawFile = [line.rstrip() for line in open(path)]
        self.path = path

        self.list0 = ['nframe', 'dof', 'moldof', 'performance', 'mass']
        self.performance = 0.0
        self.nframe = 0
        self.mass = mass
        self.dof = dof
        self.moldof = moldof

    def fmtstr(self, p) -> str:
        shared = '\n{:20s}{:12.3f}'
        blank = '            '
        if p == 'performance':
            return shared + blank + ' ns/day'
        elif p == 'mass':
            return shared + blank + ' a.u.'
        elif p == 'nframe':
            return '\n{:20s}{:12d}'
        elif p == 'dof' or p == 'moldof':
            return '\n{:20s}{:12d}'
        else:
            return shared

    def brief(self) -> str:
        out = '{}'.format(self.path)
        for k in self.list0:
            foo = getattr(self, k)
            if foo != 0:
                out = out + self.fmtstr(k).format(k, getattr(self, k))
        return out


class ALogData:
    def __init__(self):
        self.list1 = ['time', 'potential', 'kinetic',
            'xbox', 'ybox', 'zbox', 'alpha', 'beta', 'gamma']
        self.list2 = ['volume', 'density', 'temperature']
        self.list3 = ['pressure', 'stress']
        self.list0 = self.list1 + self.list2 + self.list3
        for k in self.list0:
            setattr(self, k, [])

    def keep(self, start: int, stop: int, step: int):
        for k in self.list1 + self.list2:
            p = getattr(self, k)
            if len(p):
                setattr(self, k, np.array(p[start:stop:step]))

    def keep2(self, start: int, stop: int, step: int):
        for k in self.list3:
            p = getattr(self, k)
            lenp = len(p)
            if stop == -1:
                stop2 = lenp
            else:
                stop2 = stop
            if lenp:
                setattr(self, k, np.array(p[start:stop2:step]))

    def fmtstr(self, p) -> str:
        shared = '\n{:20s}{:12.3f}{:12.3f}'
        if p == 'potential':
            return shared + ' kcal/mol'
        elif p == 'temperature':
            return shared + ' Kelvin'
        else:
            return shared

    def brief(self) -> str:
        out = ''
        len1 = len(self.potential)
        out = out + '\n{:20s}{:12d}'.format('ndata', len1)
        for k in self.list0:
            if k == 'time':
                continue

            foo = getattr(self, k)
            if len(foo) == 0:
                continue
            elif k == 'stress':
                avg, std = foo.mean(axis=0), foo.std(axis=0)
                fmt = '\n{:20s}{:12.3f}{:12.3f}{:12.3f}'
                out = out + '\n'
                out = out + fmt.format('stress-x', avg[0], avg[1], avg[2])
                out = out + fmt.format('stress-y', avg[3], avg[4], avg[5])
                out = out + fmt.format('stress-z', avg[6], avg[7], avg[8])
                out = out + fmt.format('stress-std-x', std[0], std[1], std[2])
                out = out + fmt.format('stress-std-y', std[3], std[4], std[5])
                out = out + fmt.format('stress-std-z', std[6], std[7], std[8])
            else:
                avg, std = foo.mean(), foo.std()
                out = out + self.fmtstr(k).format(k, avg, std)
        return out

    def convert(self, gdata):
        for i in range(gdata.nframe):
            a, b, c = self.xbox[i], self.ybox[i], self.zbox[i]
            alpha, beta, gamma = self.alpha[i], self.beta[i], self.gamma[i]
            ekin = self.kinetic[i]
            volume = StaticMethod.volume(a, b, c, alpha, beta, gamma)
            density = StaticMethod.density(gdata.mass, volume)
            temper = StaticMethod.temperature(ekin, gdata.dof)
            self.volume.append(volume)
            self.density.append(density)
            self.temperature.append(temper)

        for k in self.list0:
            setattr(self, k, np.array(getattr(self, k)))


class MLogData:
    def __init__(self):
        self.list0 = ['molkin', 'molT', 'molP', 'molStress']
        for k in self.list0:
            setattr(self, k, [])

    def convert(self, gdata):
        for i in range(len(self.molkin)):
            ekin = self.molkin[i]
            temper = StaticMethod.temperature(ekin, gdata.moldof)
            self.molT.append(temper)

        for k in self.list0:
            setattr(self, k, np.array(getattr(self, k)))

    def keep2(self, start: int, stop: int, step: int):
        for k in self.list0:
            p = getattr(self, k)
            lenp = len(p)
            if stop == -1:
                stop2 = lenp
            else:
                stop2 = stop
            if lenp:
                setattr(self, k, np.array(p[start:stop2:step]))

    def brief(self) -> str:
        out = ''
        out = '\n{:20s}{:12d}'.format('ndata/mol', len(self.molkin))
        for k in self.list0:
            foo = getattr(self, k)
            if len(foo) == 0:
                continue
            elif k == 'molStress':
                avg, std = foo.mean(axis=0), foo.std(axis=0)
                fmt = '\n{:20s}{:12.3f}{:12.3f}{:12.3f}'
                out = out + '\n'
                out = out + fmt.format('molStress-x', avg[0], avg[1], avg[2])
                out = out + fmt.format('molStress-y', avg[3], avg[4], avg[5])
                out = out + fmt.format('molStress-z', avg[6], avg[7], avg[8])
                out = out + fmt.format('molStress-std-x', std[0], std[1], std[2])
                out = out + fmt.format('molStress-std-y', std[3], std[4], std[5])
                out = out + fmt.format('molStress-std-z', std[6], std[7], std[8])
            else:
                avg, std = foo.mean(), foo.std()
                out = out + '\n{:20s}{:12.3f}{:12.3f}'.format(k, avg, std)
        return out


class LogData:
    def __init__(self, mass=0.0, dof=0, moldof=0, file='', dir='', keep1=[1, -1, 1], keep2=[1, -1, 1]):
        self.gdata = GLogData(mass, dof, moldof, file, dir)
        self.adata = ALogData()
        self.mdata = None if moldof == 0 else MLogData()

        for line in self.gdata.rawFile:
            if 'Current Time ' in line:
                # ps
                foo = float(line.split()[2])
                self.adata.time.append(foo)
            elif 'Current Potential ' in line:
                # kcal/mol
                foo = float(line.split()[2])
                self.adata.potential.append(foo)
            elif 'Current Kinetic ' in line:
                # kcal/mol
                foo = float(line.split()[2])
                self.adata.kinetic.append(foo)
            elif 'Lattice Lengths ' in line:
                # angstrom
                vs = line.split()
                self.adata.xbox.append(float(vs[2]))
                self.adata.ybox.append(float(vs[3]))
                self.adata.zbox.append(float(vs[4]))
            elif 'Lattice Angles ' in line:
                # degree
                vs = line.split()
                self.adata.alpha.append(float(vs[2]))
                self.adata.beta.append(float(vs[3]))
                self.adata.gamma.append(float(vs[4]))
            elif 'Performance: ' in line:
                # ns/day
                self.gdata.performance = float(line.split()[2])
            elif 'Atomic Pressure Atm ' in line:
                # atm
                vs = line.split()
                self.adata.pressure.append(float(vs[3]))
            elif 'Group Pressure Atm ' in line:
                # atm
                vs = line.split()
                self.mdata.molP.append(float(vs[3]))
            elif 'Group Kinetic ' in line:
                # kcal/mol
                vs = line.split()
                self.mdata.molkin.append(float(vs[2]))
            elif 'Atomic Pressure x- Atm ' in line:
                # atm
                vs = line.split()
                self.adata.stress.append([float(vs[4]), float(vs[5]), float(vs[6])])
            elif 'Atomic Pressure y- Atm ' in line:
                vs = line.split()
                self.adata.stress[-1].append(float(vs[4]))
                self.adata.stress[-1].append(float(vs[5]))
                self.adata.stress[-1].append(float(vs[6]))
            elif 'Atomic Pressure z- Atm ' in line:
                vs = line.split()
                self.adata.stress[-1].append(float(vs[4]))
                self.adata.stress[-1].append(float(vs[5]))
                self.adata.stress[-1].append(float(vs[6]))
            elif 'Group Pressure x- Atm ' in line:
                # atm
                vs = line.split()
                self.mdata.molStress.append([float(vs[4]), float(vs[5]), float(vs[6])])
            elif 'Group Pressure y- Atm ' in line:
                vs = line.split()
                self.mdata.molStress[-1].append(float(vs[4]))
                self.mdata.molStress[-1].append(float(vs[5]))
                self.mdata.molStress[-1].append(float(vs[6]))
            elif 'Group Pressure z- Atm ' in line:
                vs = line.split()
                self.mdata.molStress[-1].append(float(vs[4]))
                self.mdata.molStress[-1].append(float(vs[5]))
                self.mdata.molStress[-1].append(float(vs[6]))

        self.gdata.nframe = len(self.adata.time)
        self.adata.convert(self.gdata)
        if self.mdata:
            self.mdata.convert(self.gdata)

        start1, stop1, step1 = keep1
        start1 = start1 - 1
        if stop1 == -1:
            stop1 = self.gdata.nframe
        start2, stop2, step2 = keep2
        start2 = start2 - 1
        self.adata.keep(start1, stop1, step1)
        self.adata.keep2(start2, stop2, step2)
        if self.mdata:
            self.mdata.keep2(start2, stop2, step2)

    def brief(self) -> str:
        out = ''
        out = out + self.gdata.brief()
        out = out + '\n'
        out = out + self.adata.brief()

        if self.mdata:
            out = out + '\n'
            out = out + self.mdata.brief()
        return out

    def __str__(self) -> str:
        return self.brief()
