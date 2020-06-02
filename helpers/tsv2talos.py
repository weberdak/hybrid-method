import pandas as pd
import numpy as np
import argparse

# Written by D. K. Weber, Veglia Lab. (last revised June 2 2020)
# Script to convert CSV file of chemical shifts into TALOS-N input
# https://spin.niddk.nih.gov/bax/nmrserver/talosn/


def chunks(l, n):
    """Yield successive n-sized chunks from l. From stackoverflow 312443"""
    for i in range(0, len(l), n):
        yield l[i:i + n]

        
def parse_args():
    parser = argparse.ArgumentParser(description='Initial folding for membrane protein structure.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-i', '--infile', type=str,
        help='TSV file of chemical shifts.'
    )
    parser.add_argument(
        '-o', '--outfile', type=str,
        help='Output file (TALOS-N input).'
    )
    parser.add_argument(
        '-s', '--sequence', type=str,
        help='Amino acid sequence. 1 letter code'
    )
    parser.add_argument(
        '-r', '--start_id', type=int,
        help='Residue ID of first amino acid.', default=1
    )
    args = parser.parse_args()
    return args
    

def main():

    # Handle arguments
    args = parse_args()
    
    # Read shifts into dataframe
    df = pd.read_csv(args.infile,header=0,delimiter='\t')

    # Start writing data
    f = open(args.outfile,'w')
    sequence = chunks(args.sequence, 10)
    sequence = ' '.join(sequence)

    # Header
    f.write('REMARK TALOS input file\n\n')
    f.write('DATA FIRST_RESID {}\n\n'.format(args.start_id))
    f.write('DATA SEQUENCE {}\n\n'.format(sequence))
    f.write('VARS   RESID RESNAME ATOMNAME SHIFT\n')
    f.write('FORMAT %4d   %1s     %4s      %8.3f\n\n')

    for index,row in df.iterrows():
        shifts = [s for s in row.keys()]
        shifts.remove('RESID')
        shifts.remove('RESNAME')
        for s in shifts:
            resid = str(row['RESID']).rjust(4,' ')
            resname = str(row['RESNAME'])
            name = str(s).ljust(4,' ')
            if row[s] > 0:
                shift = round(row[s], 3)
                shift = str(shift).rjust(8,' ')
                f.write(resid+' '+resname+' '+name+' '+' '+shift+'\n')
                #print(resid+' '+resname+' '+name+' '+' '+shift)

    f.close()



if __name__ == '__main__':
    main()
