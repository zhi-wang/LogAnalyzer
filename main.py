#!/usr/bin/env python3


from logdata import LogData, Property
import sys


if __name__ == '__main__':
    argv = []
    for r in sys.argv:
        argv.append(r)
    for r in sys.stdin:
        argv.append(r.rstrip())
    file, directory = '', ''
    mass, dof = 0.0, 0.0
    if file == '':
        file = argv[1]
        mass, dof = float(argv[2]), float(argv[3])
    d = LogData(mass=mass, dof=dof, file=file, dir=directory)
    print(d)
