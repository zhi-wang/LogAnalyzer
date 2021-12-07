from enum import IntEnum
from math import pi, cos, sin, sqrt
import numpy as np
import os


class _Property:
    is_array = 0x0001

    nframe = (1<<1)
    mass = (2<<1)
    dof = (3<<1)
    path = (4<<1)
    performance = (5<<1)

    time = (11<<1) | is_array
    potential = (12<<1) | is_array
    kinetic = (13<<1) | is_array
    xbox = (14<<1) | is_array
    ybox = (15<<1) | is_array
    zbox = (16<<1) | is_array
    alpha = (17<<1) | is_array
    beta = (18<<1) | is_array
    gamma = (19<<1) | is_array

    volume = (20<<1) | is_array
    density = (21<<1) | is_array
    temperature = (22<<1) | is_array

class Property(IntEnum):
    nframe = _Property.nframe
    mass = _Property.mass
    dof = _Property.dof
    path = _Property.path
    performance = _Property.performance

    time = _Property.time
    potential = _Property.potential
    kinetic = _Property.kinetic
    xbox = _Property.xbox
    ybox = _Property.ybox
    zbox = _Property.zbox
    alpha = _Property.alpha
    beta = _Property.beta
    gamma = _Property.gamma

    volume = _Property.volume
    density = _Property.density
    temperature = _Property.temperature


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


class LogData:
    def __init__(self, mass=0.0, dof=0.0, file='', dir=''):
        path = os.path.join(os.path.expanduser(dir), file)
        self._rawFile = [line.strip() for line in open(path)]
        self._path = path

        self._time, self._potential, self._kinetic = [], [], []
        self._xbox, self._ybox, self._zbox = [], [], []
        self._alpha, self._beta, self._gamma = [], [], []
        self._performance = 0.0
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
            elif 'Performance: ' in line:
                # ns/day
                self._performance = float(line.split()[2])

        self._nframe = len(self._time)
        # atomic unit
        self._mass = float(mass)
        self._dof = float(dof)

        self._volume, self._density, self._temperature = [], [], []
        for i in range(self._nframe):
            a, b, c = self._xbox[i], self._ybox[i], self._zbox[i]
            alpha, beta, gamma = self._alpha[i], self._beta[i], self._gamma[i]
            ekin = self._kinetic[i]
            volume = StaticMethod.volume(a, b, c, alpha, beta, gamma)
            density = StaticMethod.density(self._mass, volume)
            temper = StaticMethod.temperature(ekin, self._dof)
            self._volume.append(volume)
            self._density.append(density)
            self._temperature.append(temper)

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
            Property.temperature: np.array(self._temperature[start:stop:step]),
        }
        self._properties[Property.nframe] = len(self._properties[Property.time])
        self._properties[Property.mass] = self._mass
        self._properties[Property.dof] = self._dof
        self._properties[Property.path] = self._path
        self._properties[Property.performance] = self._performance

    def discard(self, first_n=0):
        self.keep(first_n, -1, 1)

    def get(self, p: Property):
        return self._properties[p]

    def mean(self, p: Property):
        array = self.get(p)
        mean, std = None, None
        if p & _Property.is_array:
            mean, std = array.mean(), array.std()
        else:
            mean = array
        return mean, std
