import numpy as np

order = 0.80
error = 10.0
error_dc = 1.0
flip = 1.0

data = 'ossnmr_restraints.dat'
outfile = 'ossnmr_cs.tbl'
outfile_dc = 'ossnmr_dc.tbl'

f = open(outfile,'w')
f_dc = open(outfile_dc,'w')

resids = np.genfromtxt(data,usecols=(0),dtype=int)
cs_isos = np.genfromtxt(data,usecols=(2),dtype=float)
cs_anis = np.genfromtxt(data,usecols=(3),dtype=float)
dcs = np.genfromtxt(data,usecols=(4),dtype=float)

for r, cs_iso, cs_ani, dc in zip(resids,cs_isos,cs_anis,dcs):
    cs_corr = (1/order*flip)*(cs_ani-cs_iso)
    dc_corr = (1/order*flip)*dc

    f.write('assign\n')
    f.write(' ( resid 9999 and name OO )\n')
    f.write(' ( resid 9999 and name Z )\n')
    f.write(' ( resid 9999 and name X )\n')
    f.write(' ( resid 9999 and name Y )\n')
    f.write(' ( resid {} and name N )\n'.format(r))
    f.write(' ( resid {} and name HN )\n'.format(r))
    f.write(' ( resid {} and name C )\n'.format(r-1))
    f.write(' {0:.3f} {1:.3f} {1:.3f}\n'.format(cs_corr,error))

    f_dc.write('assign\n')
    f_dc.write(' ( resid 9999 and name OO )\n')
    f_dc.write(' ( resid 9999 and name Z )\n')
    f_dc.write(' ( resid 9999 and name X )\n')
    f_dc.write(' ( resid 9999 and name Y )\n')
    f_dc.write(' ( resid {} and name N )\n'.format(r))
    f_dc.write(' ( resid {} and name HN )\n'.format(r))
    f_dc.write(' {0:.3f} {1:.3f} {1:.3f}\n'.format(dc_corr,error_dc))

f.close()
f_dc.close()
