#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=24
#SBATCH --mem=16g
#SBATCH --tmp=16g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=dweber@umn.edu

cd /home/vegliag/dweber/xplor_data/20210325_SLN

# Folding
echo "Folding protein."
xplor -py -smp 24 -o logfile.out hybrid.py \
      --structure_in      summary.fold/top.0.sa \
      --DC_NH_work        input/sln_dc.tbl \
      --CSA_N1_work       input/sln_cs.tbl \
      --CSA_N1_gly_work   input/sln_cs_gly.tbl \
      --HBDA              input/sln.hbda.tbl \
      --NOE               input/sln.hbnoe.tbl \
      --DIHE              input/sln.dihe.tbl \
      --DC_NH_max         10.735 \
      --nstructures       1000 \
      --CSA_N1_tensor     57.3 81.2 228.1 \
      --CSA_N1_tensor_gly 45.6 66.3 211.6 \
      --CSA_N1_beta       -17.0 \
      --CSA_N1_beta_gly   -21.6 \
      --tm_domain         12 24 \
      --immx_thickness    25.72 \
      --immx_nparameter   10 \
      --w_slf             5 \
      --w_r               3 \
      --seed              4534 \
      --highTempSteps     1000 \
      --initialTemp       350 \
      --finalTemp         2.5 \
      --stepTemp          1.25 \
      --annealSteps       201 \
      --unfold            no \
      --resetCenter       yes \
      --repelStart        no \
      --ezPot             "resid 0:31"


# Mark best structures and move to folder
echo "Getting best refined structures according to XPLOR."
getBest -num 10 -symlinks
rm -rf out.refine
mkdir out.refine
mv *.sa* out.refine/
mv logfile* out.refine/


# Run summary statistics
echo "Getting best refined structures using summary.py"
rm -rf summary.refine
mkdir summary.refine
cd summary.refine/
python3 ../summary.py \
	--folders ../out.refine* \
	--terms DIPL_w CS_w CDIH BOND ANGL IMPR EEFX \
	--R_dc_work amide_NH_work \
	--R_dc_err 0.5 \
	--R_csa_work amide_N1_work amide_N1_gly_work \
	--R_csa_err 5.0 > summary.out
cd ../

