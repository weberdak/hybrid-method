import argparse
import re
import random
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Divide XPLOR-NIH CSA and DC restraint tables into R-work and R-free.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-i', '--restraint_tables', type=str, nargs='+', 
        help='All DC and CSA restraint tables.'
    )
    parser.add_argument(
        '-r', '--rfree', type=float, default=20.0, 
        help='Percentage of restraints to reserve for R-free. Default: 20.0'
    )
    args = parser.parse_args()
    return args


def unique(list_in):  
    return list(set(list_in))


def main():    
    args = parse_args()

    # Init restraint dictionary and counter
    restraint = dict()
    c = -1
    for tbl in args.restraint_tables:
        fi = open(tbl, 'r')

        for line in fi.readlines():
            
            # Init new restraint
            if re.search('assign', line):
                c += 1
                restraint[c] = dict()
                restraint[c]['lines'] = []
                restraint[c]['file'] = tbl

            # Record lines
            restraint[c]['lines'].append(line)
            
        fi.close()

        
    # Randomize restraints
    n_work = len(restraint.keys()) *  (1 - args.rfree / 100)
    r_work = random.sample(restraint.keys(), int(n_work))
    r_free = [ r for r in restraint.keys() if r not in r_work ]
    r_work = sorted(r_work)
    r_free = sorted(r_free)

    
    # Report
    print('Detected {} total restraints'.format(len(restraint.keys())))
    print('Reserving {} ({}%) for R-free'.format(len(r_free), args.rfree))
    print('R-Work: {}'.format(r_work))
    print('R-free: {}'.format(r_free))

    
    # Write R-work restraints
    # Init output pointers
    fnames = []
    for i in r_work:
        fnames.append(restraint[i]['file'])
    fnames = list(set(fnames))
    fouts = dict()
    for fname in fnames:
        fsplit = os.path.splitext(fname)
        fout_str = fsplit[0]+'.work'+fsplit[1]
        fo = open(fout_str, 'w')
        fouts[fname] = dict()
        fouts[fname]['count'] = 0
        fouts[fname]['fout'] = fo
        fouts[fname]['fout_str'] = fout_str

    # Write restraints
    for i in r_work:
        f = restraint[i]['file']
        for line in restraint[i]['lines']:
            fouts[f]['fout'].write(line)
        fouts[f]['count'] += 1

    # Close outfiles and report
    for f in fnames:
        fouts[f]['fout'].close()
        print('{} restraints written to {}'.format(fouts[f]['count'], fouts[f]['fout_str']))


    # Write R-free restraints
    # Init output pointers
    fnames = []
    for i in r_free:
        fnames.append(restraint[i]['file'])
    fnames = list(set(fnames))
    fouts = dict()
    for fname in fnames:
        fsplit = os.path.splitext(fname)
        fout_str = fsplit[0]+'.free'+fsplit[1]
        fo = open(fout_str, 'w')
        fouts[fname] = dict()
        fouts[fname]['count'] = 0
        fouts[fname]['fout'] = fo
        fouts[fname]['fout_str'] = fout_str

    # Write restraints
    for i in r_free:
        f = restraint[i]['file']
        for line in restraint[i]['lines']:
            fouts[f]['fout'].write(line)
        fouts[f]['count'] += 1

    # Close outfiles and report
    for f in fnames:
        fouts[f]['fout'].close()
        print('{} restraints written to {}'.format(fouts[f]['count'], fouts[f]['fout_str']))

        
if __name__ == '__main__':
    main()

