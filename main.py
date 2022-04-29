#!/usr/bin/env python3


from logdata import LogData
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True)
    parser.add_argument('--dir', type=str, required=False)
    parser.add_argument('--mass', type=float, required=True)
    parser.add_argument('--dof', type=int, required=True)
    parser.add_argument('--moldof', type=int, required=False)
    parser.add_argument('--keep1', type=int, nargs=3, required=False)
    parser.add_argument('--keep2', type=int, nargs=3, required=False)
    args = parser.parse_args()

    mass = args.mass
    dof = args.dof
    moldof = 0
    if args.moldof:
        moldof = args.moldof
    file = args.file
    directory = ''
    if args.dir:
        directory = args.dir
    keep1 = [1, -1, 1]
    if args.keep1:
        keep1 = args.keep1
    keep2 = [1, -1, 1]
    if args.keep2:
        keep2 = args.keep2
    d = LogData(mass=mass, dof=dof, moldof=moldof, file=file, dir=directory, keep1=keep1, keep2=keep2)
    print(d)
