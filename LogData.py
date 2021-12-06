from enum import Enum
from math import gamma, pi, cos, sin, sqrt
import numpy as np
import os

from numpy.lib.function_base import select

class Property(Enum):
    time = 1
    potential = 2
    kinetic = 3
    xbox = 4
    ybox = 5
    zbox = 6
    alpha = 7
    beta = 8
    gamma = 9

    volume = 10
    density = 11
    temperature = 12

    nframe = 0

class LogData:
    _avogadro = 6.02214076e23
    _gasconst = 1.987204259e-3

    def _calc_volume(a, b, c, alpha, beta, gamma):
        radian = 180.0 / pi
        al, be, ga = alpha/radian, beta/radian, gamma/radian
        ca, cb, cg = cos(al), cos(be), cos(ga)
        sg = sin(ga)
        vol = a*b*c*sqrt(sg**2 - ca**2 - cb**2 + 2.0*ca*cb*cg)
        return vol

    def _calc_density(amass, volume):
        avogadro = LogData._avogadro
        density = (amass/volume)*(1.e24/avogadro)
        return density

    def _calc_temperature(ekin, dof):
        gasconst = LogData._gasconst
        temper = 2.0*ekin/(dof*gasconst)
        return temper

    def __init__(self, amass=0.0, dof=0.0, file='', dir=''):
        path = os.path.join(dir, file)
        self._rawFile = [line.strip() for line in open(path)]
        self._time, self._potential, self._kinetic = [], [], []
        self._xbox, self._ybox, self._zbox = [], [], []
        self._alpha, self._beta, self._gamma = [], [], []
        for line in self._rawFile:
            if 'Current Time ' in line:
                # ps
                foo = float(line.split()[2])
                self._time.append(foo)
            elif 'Current Potential ' in line:
                # kcal/mol
                foo = float(line.split()[2])
                self._potential.append(foo)
            elif 'Current Kinetic ' in line:
                # kcal/mol
                foo = float(line.split()[2])
                self._kinetic.append(foo)
            elif 'Lattice Lengths ' in line:
                # angstrom
                vs = line.split()
                self._xbox.append(float(vs[2]))
                self._ybox.append(float(vs[3]))
                self._zbox.append(float(vs[4]))
            elif 'Lattice Angles ' in line:
                # degree
                vs = line.split()
                self._alpha.append(float(vs[2]))
                self._beta.append(float(vs[3]))
                self._gamma.append(float(vs[4]))

        self._nframe = len(self._time)
        self.amass = float(amass) # atomic unit
        self.dof = float(dof)

        self._volume, self._density, self._kelvin = [], [], []
        for i in range(self._nframe):
            a, b, c = self._xbox[i], self._ybox[i], self._zbox[i]
            alpha, beta, gamma = self._alpha[i], self._beta[i], self._gamma[i]
            ekin = self._kinetic[i]
            volume = LogData._calc_volume(a, b, c, alpha, beta, gamma)
            density = LogData._calc_density(self.amass, volume)
            temper = LogData._calc_temperature(ekin, self.dof)
            self._volume.append(volume)
            self._density.append(density)
            self._kelvin.append(temper)

        self._properties = {}
        self.keep()

    def keep(self, start=0, stop=-1, step=1):
        if stop == -1:
            stop = self._nframe

        self._properties = {
            Property.time: np.array(self._time[start:stop:step]),
            Property.potential: np.array(self._potential[start:stop:step]),
            Property.kinetic: np.array(self._kinetic[start:stop:step]),
            Property.xbox : np.array(self._xbox[start:stop:step]),
            Property.ybox : np.array(self._ybox[start:stop:step]),
            Property.zbox : np.array(self._zbox[start:stop:step]),
            Property.alpha : np.array(self._alpha[start:stop:step]),
            Property.beta : np.array(self._beta[start:stop:step]),
            Property.gamma : np.array(self._gamma[start:stop:step]),

            Property.volume: np.array(self._volume[start:stop:step]),
            Property.density: np.array(self._density[start:stop:step]),
            Property.temperature: np.array(self._kelvin[start:stop:step]),
        }
        self._properties[Property.nframe] = len(self._properties[Property.time])

    def discard(self, first_n=0):
        self.keep(first_n, -1, 1)

    def get(self, p: Property):
        return self._properties[p]

    def mean(self, p: Property):
        array = self.get(p)
        mean, std = None, None
        if p != Property.nframe:
            mean, std = array.mean(), array.std()
        return mean, std

if __name__ == '__main__':
    pass
