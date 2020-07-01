# fakeHelix.py
# ------------
# Write DIHE, NOE and HBDA restraint files for a specified helical segment.
#
# Example Usage
# -------------
# python3 fakeHelix.py --start 6 --stop 28 \
#    --sequence MGINTRELFLNFTIVLITVILMWLLVRSYQY --start_id 1 --out_prefix sln
#
# Will write dihedral (sln.dihe.tbl), noe (sln.hbnoe.tbl) and hydrogen bond (sln.hbda.tbl)
# restraint files. Flags --sequence, --start_id (default 1) and --out_prefix (default fakeHelix)
# are optional. The --sequence input is required if prolines are to be ignored.
#
# Written by Daniel Weber, Veglia Lab
# Last revised June 30 2020

import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Generate artificial phi/psi backbone dihedral restraints for XPLOR-NIH.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--start', type=int,
        help='First amino acid residue ID of helical segment'
    )
    parser.add_argument(
        '--stop', type=int,
        help='Last amino acid residue ID'
    )
    parser.add_argument(
        '--out_prefix', type=str,
        help='Output prefix of restraint files.', default='fakeHelix'
    )
    parser.add_argument(
        '--sequence', type=str,
        help='Amino acid sequence. 1 letter code. Optional.', default=''
    )
    parser.add_argument(
        '--start_id', type=int,
        help='Residue ID of first amino acid in input sequence.', default=1
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    f = open(args.out_prefix+'.dihe.tbl','w')
    f_n = open(args.out_prefix+'.hbnoe.tbl','w')
    f_h = open(args.out_prefix+'.hbda.tbl','w')

    # Convert sequence to dictionary
    if args.sequence:
        seq = [ a for a in args.sequence ]
        rid = [ i for i in range(args.start_id,args.start_id+len(seq)+1) ]
        d = { i: a for i,a in zip(rid,seq) }
    
    for i in range(args.start,args.stop+1):

        # Write phi angle
        f.write('assign\n')
        f.write('  ( resid {} and name    C)\n'.format(i-1))
        f.write('  ( resid {} and name    N)\n'.format(i))
        f.write('  ( resid {} and name   CA)\n'.format(i))
        f.write('  ( resid {} and name    C)\n'.format(i))
        f.write('  1.0 -63.00 30.00 2\n')

        # Write psi angle
        f.write('assign\n')
        f.write('  ( resid {} and name    N)\n'.format(i))
        f.write('  ( resid {} and name   CA)\n'.format(i))
        f.write('  ( resid {} and name    C)\n'.format(i))
        f.write('  ( resid {} and name    N)\n'.format(i+1))
        f.write('  1.0 -42.00 30.00 2\n')

        # Write NOE and HBDA restraints
        if i >= args.start+4:
            # Check sequence for prolines
            if args.sequence:
                if d[i] != 'P':
                    f_n.write('assign (resid {} and name N)(resid {} and name O) 2.00 0.1 0.1\n'.format(i,i-4))
                    f_h.write('assign (resid {} and name N)(resid {} and name HN) (resid {} and name O)\n'.format(i,i,i-4))
                else:
                    print('Skipping residue {}{}'.format(d[i],i))
            # Sequence not specified, don't bother checking
            else:
                f_n.write('assign (resid {} and name N)(resid {} and name O) 2.00 0.1 0.1\n'.format(i,i-4))
                f_h.write('assign (resid {} and name N)(resid {} and name HN) (resid {} and name O)\n'.format(i,i,i-4))

    f.close()
    f_n.close()
    f_h.close()
    
if __name__ == '__main__':
    main()
