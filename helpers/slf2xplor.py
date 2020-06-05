import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Initial folding for membrane protein structure.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-i', '--infile', type=str,
        help='TSV file of OS-ssNMR chemical shifts and dipolar couplings.'
    )
    parser.add_argument(
        '-o', '--out_prefix', type=str,
        help='Output file prefix for CS and DC restraint tables.'
    )
    parser.add_argument(
        '--order', type=float,
        help='General order parameter for scaling CS and DC restraints. Default: 1.0', default=1.0
    )
    parser.add_argument(
        '--align_order', type=float, choices=[1.0, -0.5],
        help='Alignment order parameter. 1.0 for flipped bicelle or -0.5 for unflipped bicelle. Default: 1.0',
        default=1.0
    )
    parser.add_argument(
        '--pas', type=float, nargs='+',
        help='Backbone amide 15N shift tensor. Default: 57.3 81.2 228.1',
        default=(57.3,81.2,228.1)
    )
    parser.add_argument(
        '--pas_gly', type=float, nargs='+',
        help='Backbone amide 15N shift tensor for glycine. Default: 45.6 66.3 211.6',
        default=(45.6,66.3,211.6)
    )
    parser.add_argument(
        '--error_csa', type=float,
        help='Error associated with chemical shifts in ppm. Default: 5.0', default=5.0
    )
    parser.add_argument(
        '--error_dc', type=float,
        help='Error associated with chemical shifts in ppm. Default: 0.5', default=0.5
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    
    # CSA and DC parameters
    order = args.order
    error_csa = args.error_csa
    error_dc = args.error_dc
    flip = args.align_order
    pas = args.pas
    pas_gly = args.pas_gly
    
    # Filenames
    data = args.infile
    outfile_csa = args.out_prefix+'_cs.tbl'
    outfile_csa_gly = args.out_prefix+'_cs_gly.tbl'
    outfile_dc = args.out_prefix+'_dc.tbl'

    # Read data
    resids = np.genfromtxt(data,usecols=(0),dtype=int)
    resnames = np.genfromtxt(data,usecols=(1),dtype=str)
    css = np.genfromtxt(data,usecols=(2),dtype=float)
    dcs = np.genfromtxt(data,usecols=(3),dtype=float)

    # Initialize output files
    print('\nWriting the following XPLOR-NIH inputs...')
    print('CSA restraints will be written to {}'.format(outfile_csa))
    print('DC restraints will be written to {}'.format(outfile_dc))
    f_csa = open(outfile_csa,'w')
    f_dc = open(outfile_dc,'w')
    if 'G' in resnames:
        print('Glycines detected! CSA restraints will be written to {}'.format(outfile_csa_gly))
        f_csa_gly = open(outfile_csa_gly,'w')
        
    # Calculate isotropic chemical shifts
    iso = np.mean(pas)
    iso_gly = np.mean(pas_gly)

    # Verbose
    print('\nRestraints will be scaled with the following parameters...')
    print('Backbone 15N principal components: {}'.format(pas))
    print('Isotropic shift subtracted: {0:.3f}'.format(iso))
    print('Backbone 15N principal components (glycine): {}'.format(pas_gly))
    print('Isotropic shift subtracted (glycine): {0:.3f}'.format(iso_gly))
    print('Error assigned to all CSA restraints: {} ppm'.format(error_csa))
    print('Error assigned to all DC restraints: {} ppm'.format(error_dc))
    flipmsg = { 1.0: 'Flipped Bicelle/Glass Plate', -0.5: 'Unflipped Bicelle' }
    print('Alignment order: {} ({})'.format(flip,flipmsg[flip]))
    
    print('\nScaling CSA restraints (reduced form)...')
    print('Residue\tCS\tRestraint')
    # Read and scale chemical shifts and couplings
    for r, rn, cs in zip(resids,resnames,css):

        # Start making restraints
        if not np.isnan(cs):
            if rn == 'G':
                cs_corr = ((cs-iso_gly)/(order*flip))
                f = f_csa_gly
                print('{0}{1}*\t{2:.3f}\t{3:.3f}'.format(rn,r,cs,cs_corr))
            else:
                cs_corr = ((cs-iso)/(order*flip))
                f = f_csa
                print('{0}{1}\t{2:.3f}\t{3:.3f}'.format(rn,r,cs,cs_corr))

            # Write to file
            f.write('assign\n')
            f.write(' ( resid 9999 and name OO )\n')
            f.write(' ( resid 9999 and name Z )\n')
            f.write(' ( resid 9999 and name X )\n')
            f.write(' ( resid 9999 and name Y )\n')
            f.write(' ( resid {} and name N )\n'.format(r))
            f.write(' ( resid {} and name HN )\n'.format(r))
            f.write(' ( resid {} and name C )\n'.format(r-1))
            f.write(' {0:.3f} {1:.3f} {1:.3f}\n'.format(cs_corr,error_csa))


    # Write DC restraints
    print('\nScaling DC restraints...')
    for r, rn, dc in zip(resids,resnames,dcs):

        if not np.isnan(dc):
            dc_corr = abs(dc/(order*flip))
            print('{0}{1}\t{2:.3f}\t{3:.3f}'.format(rn,r,dc,dc_corr))
            f_dc.write('assign\n')
            f_dc.write(' ( resid 9999 and name OO )\n')
            f_dc.write(' ( resid 9999 and name Z )\n')
            f_dc.write(' ( resid 9999 and name X )\n')
            f_dc.write(' ( resid 9999 and name Y )\n')
            f_dc.write(' ( resid {} and name N )\n'.format(r))
            f_dc.write(' ( resid {} and name HN )\n'.format(r))
            f_dc.write(' {0:.3f} {1:.3f} {1:.3f}\n'.format(dc_corr,error_dc))

    f_csa.close()
    f_csa_gly.close()
    f_dc.close()

        
if __name__ == '__main__':
    main()
