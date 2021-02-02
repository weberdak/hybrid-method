import glob
import re
import argparse
import shutil
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Summary statistics for OS-ssNMR structures.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--folders', type=str, nargs='+',
        help='List of folder to scan for .sa files.'
    )
    parser.add_argument(
        '--terms', type=str, nargs='+',
        help='Terms used to calculate total energy. Default: None'
    )
    parser.add_argument(
        '--R_dc_work', type=str, nargs='+', default='',
        help='DC working restraint names to compute R factors (.sa.viols files). Default: None'
    )
    parser.add_argument(
        '--R_dc_free', type=str, nargs='+', default='',
        help='DC free restraint names to compute R factors (.sa.viols files). Default: None'
    )
    parser.add_argument(
        '--R_dc_err', type=float, default=0.5,
        help='Errors used to compute DC R factors. Default: None'
    )
    parser.add_argument(
        '--R_csa_work', type=str, nargs='+', default='',
        help='CSA working restraint names to compute R factors (.sa.viols files). Default: None'
    )
    parser.add_argument(
        '--R_csa_free', type=str, nargs='+', default='',
        help='CSA free restraint names to compute R factors (.sa.viols files). Default: None'
    )
    parser.add_argument(
        '--R_csa_err', type=float, default=5.0,
        help='Errors used to compute CSA R factors. Default: None'
    )
    parser.add_argument(
        '--top', type=int, default=10,
        help='Extract this number of top structures and summarize stats. Default: 10'
    )
    args = parser.parse_args()
    return args


def r_value(viols_file, search_terms, err, rtype):
    fi = open(viols_file, 'r')

    # Initialize lists
    resids = []
    resnames = []
    obss = []
    calcs = []
    
    # Scan for all Dipolar Coupling restraints
    mark = 0
    for line in fi.readlines():
        if mark == 2 and len(line.split()) == 0:
            mark = 0
        
        if mark == 2:
            line_s = line.split()
            if line_s[0] == '*':
                shift = 1
            else:
                shift = 0
                
            if rtype == 'DC':
                #print(line_s)
                resids.append(int(line_s[2+shift]))
                resnames.append(line_s[3+shift])
                obss.append(float(line_s[9+shift]))
                calcs.append(float(line_s[10+shift]))

            if rtype == 'CSA':
                #print(line_s)
                resids.append(int(line_s[2+shift]))
                resnames.append(line_s[3+shift])
                obss.append(float(line_s[5+shift]))
                calcs.append(float(line_s[6+shift]))

        if mark == 1 and re.search('---', line):
            mark = 2

        for t in search_terms:
            if re.search('{} '.format(t), line) and mark == 0:
                mark = 1

    # Compute R-factor
    csum = []
    for e,c in zip(obss,calcs):
        e = abs(e)
        c = abs(c)
        r = (abs(e-c)/err)**2
        csum.append(r)
    r = sum(csum)/len(obss)

    result = []
    result.append(r)
    result.append(zip(resids,resnames,obss,calcs))
    #for a,b,c,d in zip(resids,resnames,obss,calcs):
    #    print('{}\t{}\t{}\t{}'.format(a,b,c,d))
    #print(r)   
    return result


def r_output(indata):
    c = 0
    data = dict()
    mods = []
    rids = []
    for l in indata:
        for resid,resname,obs,calc in l:
            try:
                data[resid].append(str(calc).ljust(10))
            except:
                data[resid] = [resname.ljust(10), str(obs).ljust(10), str(calc).ljust(10)]
                rids.append(resid)
        mods.append(str(c).ljust(10))
        c += 1

    rids = sorted(rids)
    print('RESID'.ljust(10)+' '+'RESNAME'.ljust(10)+' '+'OBS'.ljust(10)+' '+' '.join(mods))
    for i in rids:
        print(str(i).ljust(10)+' '+' '.join(data[i]))



def main():

    # Read arguments
    args = parse_args()
    folders = args.folders
    energy_terms = args.terms
    r_dc_work = args.R_dc_work
    r_dc_free = args.R_dc_free
    r_dc_err = args.R_dc_err
    r_csa_work = args.R_csa_work
    r_csa_free = args.R_csa_free
    r_csa_err = args.R_csa_err
    top = args.top

    # Accumulate list of files from folder list
    files = []
    for folder in folders:
        for f in glob.glob('{}/*.sa'.format(folder)):
            files.append(f)

    # Collect summary data for selected energy terms
    summary = dict()
    results = []
    for f in files:

        # Init summary results
        fi = open(f, 'r')
        tot_energy = 0
        summary[f] = dict()

        # Scan each line in .sa file for energy terms
        for line in fi.readlines():
            info = re.search('REMARK summary', line)
            if info:
                line = line.split()
                if line[2] in energy_terms:
                    summary[f][line[2]] = line[3]
                    tot_energy += float(line[3])
        
        # Record file and total energy for sorting
        results.append((f, tot_energy))
        fi.close()

            
    # Sort files by energy
    print('Detected {} files'.format(len(results)))
    print('Ranking structures according to sum of terms: {}\n'.format(' '.join(energy_terms)))
    sorted_results = sorted(results, key = lambda x: x[1])

    
    # Output sorted results
    # Add fields for R_values
    if r_csa_work:
        energy_terms.append('R_CSA_w')
    if r_csa_free:
        energy_terms.append('R_CSA_f')
    if r_dc_work:
        energy_terms.append('R_DC_w')
    if r_dc_free:
        energy_terms.append('R_DC_f')


    energy_terms_ljust = [ t.ljust(10) for t in energy_terms ]
    print('Filename'.ljust(38)+' '+'TOTAL'.ljust(10)+' '+' '.join(energy_terms_ljust))
    for f,e in sorted_results:
        file_results = []
        ef = '{0:.2f}'.format(e)

        # Calculate R-values for file
        if r_csa_work:
            summary[f]['R_CSA_w'] = '{0:.4f}'.format(r_value(f+'.viols', r_csa_work, r_csa_err, 'CSA')[0])
        if r_csa_free:
            summary[f]['R_CSA_f'] = '{0:.4f}'.format(r_value(f+'.viols', r_csa_free, r_csa_err, 'CSA')[0])
        if r_dc_work:
            summary[f]['R_DC_w'] = '{0:.4f}'.format(r_value(f+'.viols', r_dc_work, r_dc_err, 'DC')[0])
        if r_dc_free:
            summary[f]['R_DC_f'] = '{0:.4f}'.format(r_value(f+'.viols', r_dc_free, r_dc_err, 'DC')[0])
        
        for term in energy_terms:
            file_results.append(str(summary[f][term]).ljust(10))
        print(f.ljust(38)+' '+str(ef).ljust(10)+' '+' '.join(file_results))

    # Print summary statistics
    print()
    c = 0
    r_csa_w_lists = []
    r_csa_f_lists = []
    r_dc_w_lists = []
    r_dc_f_lists = []
    for f,e in sorted_results[:top]:
        print('Copying {} to top.{}.sa'.format(f, c))
        shutil.copyfile(f, 'top.{}.sa'.format(c))
        if r_csa_work:
            r_csa_w_lists.append(r_value(f+'.viols', r_csa_work, r_csa_err, 'CSA')[1])
        if r_csa_free:
            r_csa_f_lists.append(r_value(f+'.viols', r_csa_free, r_csa_err, 'CSA')[1])
        if r_dc_work:
            r_dc_w_lists.append(r_value(f+'.viols', r_dc_work, r_dc_err, 'DC')[1])
        if r_dc_free:
            r_dc_f_lists.append(r_value(f+'.viols', r_dc_free, r_dc_err, 'DC')[1])
        c += 1
    print()

    # Mean and standard deviations
    print('Term\t\tMean\t\tStdev')
    for term in energy_terms:
        vals = []
        tots = []
        for f,e in sorted_results[:top]:
            vals.append(float(summary[f][term]))
            tots.append(e)
        print('{0}\t\t{1:.3f}\t\t{2:.3f}'.format(term, np.mean(vals), np.std(vals)))
    print('TOT\t\t{0:.3f}\t\t{1:.3f}'.format(np.mean(tots), np.std(tots)))
    print()

    # Output OBS vs CALC summary
    print('CSA Working Restraints:')
    if r_csa_work:
        r_output(r_csa_w_lists)

    print('\nCSA Free Restraints:')
    if r_csa_free:
        r_output(r_csa_f_lists)

    print('\nDC Working Restraints:')
    if r_dc_work:
        r_output(r_dc_w_lists)

    print('\nDC Free Restraints:')
    if r_dc_free:
        r_output(r_dc_f_lists)
    
if __name__ == '__main__':
    main()
