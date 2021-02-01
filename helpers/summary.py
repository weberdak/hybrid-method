import glob
import re

files = glob.glob('output_files/*263.sa')

energy_terms = ['DIPL', 'CS', 'CDIH',
                'BOND', 'ANGL', 'IMPR',
                'HBDA', 'HBDB', 'NOE', 'EEFX']

rcalc_dc_terms = ''


def r_value_dc(viols_file, term, err):
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
            #print(line_s)
            resids.append(int(line_s[2]))
            resnames.append(line_s[3])
            obss.append(float(line_s[9]))
            calcs.append(float(line_s[10]))
            
        if mark == 1 and re.search('---', line):
            mark = 2
            
        if re.search('Dipolar Coupling restraints in potential term: {}'.format(term), line) and mark == 0:
            mark = 1

    # Compute R-factor
    csum = []
    for e,c in zip(obss,calcs):
        e = abs(e)
        c = abs(c)
        r = (abs(e-c)/err)**2
        csum.append(r)
    r = sum(csum)/len(obss)
    
    for a,b,c,d in zip(resids,resnames,obss,calcs):
        print('{}\t{}\t{}\t{}'.format(a,b,c,d))

    print(r)
        
    return


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

    # Read R-value
    r_value_dc(f+'.viols', 'amide_NH', 0.5)
    
    # Record file and total energy for sorting
    results.append((f, tot_energy))
    fi.close()

# Sort files by energy
sorted_results = sorted(results, key = lambda x: x[1])

# Output sorted results
energy_terms_ljust = [ t.ljust(10) for t in energy_terms ]
print('#Filename'.ljust(38)+' '+'TOTAL'.ljust(10)+' '+' '.join(energy_terms_ljust))
for f,e in sorted_results:
    file_results = []
    for term in energy_terms:
        file_results.append(str(summary[f][term]).ljust(10))
    print(f.ljust(38)+' '+str(e).ljust(10)+' '+' '.join(file_results))
    
