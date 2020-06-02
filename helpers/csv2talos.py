import pandas as pd
import numpy as np

# Written by D. K. Weber
# Script to convert CSV file of chemical shifts to TALOS-N input
# https://spin.niddk.nih.gov/bax/nmrserver/talosn/


def chunks(l, n):
    """Yield successive n-sized chunks from l. From stackoverflow 312443"""
    for i in range(0, len(l), n):
        yield l[i:i + n]
        
df = pd.read_csv('dworf_shifts.dat',header=0,delimiter='\t')
sequence = 'AMAEKAGSTFSHLLVPILLLIGWIVGCIIMIYVVFS'
start = 0

f = open('dworf_shifts.tls','w')
sequence = chunks(sequence, 10)
sequence = ' '.join(sequence)

f.write('REMARK TALOS input file\n\n')
f.write('DATA FIRST_RESID {}\n\n'.format(start))
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


