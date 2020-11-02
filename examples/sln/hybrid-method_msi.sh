#!/bin/bash -l
#PBS -l walltime=12:00:00,nodes=1:ppn=24,mem=2gb 
#PBS -m abe
#PBS -M dweber@umn.edu

cd /home/vegliag/dweber/xplor_data/20201101_SLN_TEST

xplor -py -smp 24 -o logfile.out hybrid-method.py \
      --structure_in      input_xplor/sln_ext.pdb \
      --DC_NH             input_xplor/ossnmr_dc.tbl \
      --CSA_N1            input_xplor/ossnmr_cs.tbl \
      --CSA_N1_gly        input_xplor/ossnmr_cs_gly.tbl \
      --HBDA              input_xplor/sln.hbda.tbl \
      --NOE               input_xplor/sln.hbnoe.tbl \
      --DIHE              input_xplor/iso_shifts.tbl \
      --DC_NH_max         10.735 \
      --nstructures       512 \
      --CSA_N1_tensor     57.3 81.2 228.1 \
      --CSA_N1_tensor_gly 45.6 66.3 211.6 \
      --CSA_N1_beta       -17.0 \
      --CSA_N1_beta_gly   -21.6 \
      --tm_domain         7 26 \
      --immx_thickness    25.8 \
      --immx_nparameter   10 \
      --w_slf             5 \
      --w_r               3

# Mark best structures and move to folder
getBest -num 10 -symlinks
rm -rf output_files
mkdir output_files
mv *.sa* output_files/
mv logfile* output_files/




