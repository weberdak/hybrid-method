
import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Generate artificial phi/psi backbone dihedral restraints for XPLOR-NIH.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--start', type=int,
        help='First amino acid residue ID'
    )
    parser.add_argument(
        '--stop', type=int,
        help='Last amino acid residue ID'
    )
    parser.add_argument(
        '--outfile', type=str,
        help='Output name of restraint file.'
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    f = open(args.outfile,'w')
    for i in range(args.start,args.stop+1):

        # Write phi angle
        f.write('assign\n')
        f.write('  ( resid {} and name    C)\n'.format(i-1))
        f.write('  ( resid {} and name    N)\n'.format(i))
        f.write('  ( resid {} and name   CA)\n'.format(i))
        f.write('  ( resid {} and name    C)\n'.format(i))
        f.write('  1.0 -63.00 20.00 2\n')

        # Write psi angle
        f.write('assign\n')
        f.write('  ( resid {} and name    N)\n'.format(i))
        f.write('  ( resid {} and name   CA)\n'.format(i))
        f.write('  ( resid {} and name    C)\n'.format(i))
        f.write('  ( resid {} and name    N)\n'.format(i+1))
        f.write('  1.0 -42.00 20.00 2\n')
    f.close()

    
if __name__ == '__main__':
    main()
