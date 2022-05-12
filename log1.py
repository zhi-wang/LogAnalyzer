#!/usr/bin/env python3

import argparse

gasconst = 1.987204259e-3
prescon = 6.85684112e4
avogadro = 6.02214076e23
kcal = 4184.
atm = 1.01325e5
format_string1 = '{:30s}{:12.3e} {:s}'

def isothermal_compressibility_inverse_atm(v_ang3, stdv_ang3, T_K):
    '''Units:
    Isothermal Compressibility: 1/atm
    Volume, StdDevVolume: angstrom**3
    Temperature: Kelvin
    '''
    var = stdv_ang3**2
    return var/(v_ang3*gasconst*T_K*prescon)

def IsothermalCompressibility(lines, args):
    for line in lines:
        vs = line.split()
        if len(vs):
            t = args.temperature
            if vs[0] == 'volume':
                v, s = float(vs[1]), float(vs[2])
                ka = isothermal_compressibility_inverse_atm(v, s, t)
                print(format_string1.format('isothermal compressibility', ka, '1/atm'))
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, required=True)
    parser.add_argument('-p', '--property', type=str, required=True)

    parser.add_argument('-T', '--temperature', type=float, required=False)
    args = parser.parse_args()

    lines = [l.strip() for l in open(args.file)]
    if args.property == 'kappa':
        IsothermalCompressibility(lines, args)
