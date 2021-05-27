#!/bin/bash -l

## CONTROL SCRIPT FOR RUNNING HYBRID METHOD CALCULATION

# Set hybrid method directory
HYBRID_DIR=~/GitHub/hybrid-method

# Switches (set to 0 to turn off)
RUN=1
SUMMARY=1

# Run calculation
if [[ $RUN -gt 0 ]]; then
    xplor -py -smp 4 -o logfile.out $HYBRID_DIR/hybrid.py \
	  --structure_in      input/start_structure/sln_ext.pdb \
	  --DC_NH_work        input/sln_dc.work.tbl \
	  --DC_NH_free        input/sln_dc.free.tbl \
	  --CSA_N1_work       input/sln_cs.work.tbl \
	  --CSA_N1_free       input/sln_cs.free.tbl \
	  --CSA_N1_gly_work   input/sln_cs_gly.tbl \
	  --CSA_N1_gly_free   "" \
	  --HBDA              input/sln.hbda.tbl \
	  --NOE               input/sln.hbnoe.tbl \
	  --DIHE              input/sln.dihe.tbl \
	  --torsionPot        torsionDB \
	  --DC_NH_max         10.735 \
	  --nstructures       128 \
	  --CSA_N1_tensor     57.3 81.2 228.1 \
	  --CSA_N1_tensor_gly 45.6 66.3 211.6 \
	  --CSA_N1_beta       -17.0 \
	  --CSA_N1_beta_gly   -21.6 \
	  --nonbondedPot      repel \
	  --rampeefx          no \
	  --tm_domain         7 26 \
	  --immx_thickness    25.72 \
	  --immx_nparameter   10 \
	  --w_slf             5 \
	  --w_r               3 \
	  --seed              3956 \
	  --highTempSteps     20000 \
	  --initialTemp       3500.0 \
	  --finalTemp         25.0 \
	  --stepTemp          12.5 \
	  --annealSteps       201 \
	  --unfold            yes \
	  --resetCenter       yes \
	  --repelStart        no \
	  --ezPot             "resid 0:31" \
	  --relax             no \
	  --relaxTerms        CDIH NOE CS_w DIPL_w torsionDB BOND IMPR ANGL RAMA \
	  --relaxTemp         25.0 \
	  --relaxSteps        15000

    # Move output files to out folder
    rm -rf out
    mkdir out
    mv *.sa* out/
    mv logfile* out/
fi


# Run summary statistics
if [[ $SUMMARY -gt 0 ]]; then
    rm -rf summary
    mkdir summary
    cd summary/
    python3 $HYBRID_DIR/helpers/summary.py \
	    --folders ../out \
	    --terms DIPL_w CS_w BOND ANGL IMPR repel repel14 CDIH \
	    --R_dc_work amide_NH_work \
	    --R_dc_free amide_NH_free \
	    --R_dc_err 0.5 \
	    --R_csa_work amide_N1_work amide_N1_gly_work \
	    --R_csa_free amide_N1_free amide_N1_gly_free \
	    --R_csa_err 5.0 > summary.out
    cd ../
fi
