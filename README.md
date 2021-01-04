# hybrid-method

Work flow for calculating membrane protein structures using isotropic and anisotropic NMR restraints. 

## Description

For a brief introduction to combining oriented sample solid-state NMR (OS-ssNMR) and isotropic restraints from solution or MAS-ssNMR, please read/cite our new book chapter:

Weber, D. K., Larsen, E. K., Gopinath, T., & Veglia, G. (2020). Chapter 12: Hybridizing isotropic and anisotropic solid-state NMR restraints for membrane protein structure determination. In F. Separovic & M.-A. Sani (Eds.), *Solid-State NMR*. IOP Publishing. https://doi.org/10.1088/978-0-7503-2532-5ch12

This repository details a protocol to determine the structure and topology of membrane proteins by hybridizing isotropic and anisotropic NMR restraints into XPLOR-NIH simulated annealing calculations. The protocol has been updated from our previous versions ([Shi et. al., J. Biol. NMR, 2009](https://doi.org/10.1007/s10858-009-9328-9) and [Mote et al., J. Biol. NMR, 2013](http://link.springer.com/10.1007/s10858-013-9766-2)) top incorporate new functions and force fields that have since been built into XPLOR-NIH. This protocol has also been designed for users who are not expert users of XPLOR-NIH and includes the following features:

* Incorporates the [EEFx force field](http://dx.doi.org/10.1016/j.jmr.2014.03.011), with [updated parameters](https://link.springer.com/article/10.1007/s10858-016-0082-5), and the [IMMx implicit membrane model](http://dx.doi.org/10.1016/j.bpj.2015.06.047) introduced by the Marassi Lab.
* A [hybrid-method.py](hybrid-method.py) script containing a preset protocol and parameter set. This script is not modified and is instead run from a Bash configuration shell script that lists the most critical parameters.
* Supports CSA and DC restraints obtained by OS-ssNMR using bicelles (flipped and unflipped).
* [tsv2talos.py](helpers/tsv2talos.py) helper tool to convert a TSV file of chemical shifts into TALOS-N input format.
* [slf2xplor.py](helpers/slf2xplor.py) helper tool to prepare XPLOR-NIH restraint files for CSs and DCs. This applies appropriate scaling corrections to account for bicelle alignment orientation and residual protein/lipid dynamics.
* [hybridTools.tcl](helpers/hybridTools.tcl) library of VMD functions to analyze results. Including:
	* Alignment of helical segments so topology is unaffected (i.e., not applying rotations)
	* Helix tilt and azimuthal angle measurements
* [fakeHelix.py](helpers/fakeHelix.py) helper tool to generate artificial dihedral and hydrogen bonding restraints for transmembrane segments that can be safely assumed as being helical.
* All restraining potentials are optional. Just leave options blank in the configuration script if data is unavailable.

List of XPLOR-NIH potentials/classes currently applied in the [hybrid-method.py](hybrid-method.py) script:

* [varTensorTools](https://nmr.cit.nih.gov/xplor-nih/doc/current/python/ref/varTensorTools.html) for setting DC and CSA tensors
* [rdcPotTools](https://nmr.cit.nih.gov/xplor-nih/doc/current/python/ref/rdcPotTools.html) for DCs
* [csaPotTools](https://nmr.cit.nih.gov/xplor-nih/doc/current/python/ref/csaPotTools.html) for CSAs
* XplorPot: CDIH, BOND, ANGL and IMPR for dihedrals, bond lengths, bond angles and improper torsion angles.
* [TorsionDBPotTools](https://nmr.cit.nih.gov/xplor-nih/doc/current/python/ref/torsionDBPotTools.html) for database torsion angles
* [XplorPot HBDA and HBDB](https://nmr.cit.nih.gov/xplor-nih/doc/current/python/ref/protocol.html) for hydrogen bonds
* [noePotTools](https://nmr.cit.nih.gov/xplor-nih/doc/current/python/ref/noePotTools.html) for distance restraints
* [eefxPotTools](https://nmr.cit.nih.gov/xplor-nih/doc/current/python/ref/eefxPotTools.html) for non-bonded terms and implicit membrane


## Installation

This method requires a UNIX operating system (i.e., MAC-OS or Linux) to install XPLOR-NIH. All commands shown below are executed from a Bash terminal.

* Python 3 (must include numpy)
* [XPLOR-NIH 3.0](https://nmr.cit.nih.gov/xplor-nih/)
* [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) (for analysis)

To download, click the green "Code" button above and "Download ZIP". Scripts are usually executed in the directory they are copied to.

## Examples

### Sarcolipin Structure

Sarcolipin (SLN) is a single-pass membrane protein and a well-known regulator of the sarco/endoplasmic reticulum Ca<sup>2+</sup>-ATPase (SERCA) pump responsible for the relaxation of skeletal muscle cells post-contraction. SLN is a simple protein and a good starting point for learning how to hybridize isotropic and anisotropic NMR restraints into simulated annealing structure calculations using XPLOR-NIH. While this is a just a simple system, the structure and topology of SLN in lipid bilayers is crucial in order to form functional transmembrane protein-protein interactions with SERCA.



#### Step 1: Prepare dihedral restraints

An adequate set of backbone dihedral restraints are usually necessary to obtain the correct fold during simulated annealing. The hybrid-method.py script will function without these restraints, but calculations will generally run poorly. If chemical shifts are available, copy them into the TSV formatted file like [iso_shifts.dat](examples/sln/input_raw/iso_shifts.dat). The chemical shifts have been taken from MAS-ssNMR of SLN in lipid bilayers ([Mote et al., J Biol NMR, 2013](https://link.springer.com/article/10.1007/s10858-013-9766-2)). These can be prepared using a spreadsheet, then copied and pasted into a text file. Make sure the header line is correct and only standard atom names are used. Convert this file to a TALOS-N input file using the [tsv2talos.py](helpers/tsv2talos.py) helper script as follows:

	python3 tsv2talos.py -i iso_shifts.dat -o iso_shifts.tls -s MGINTRELFLNFTIVLITVILMWLLVRSYQY

This produces the [iso_shifts.tls](examples/sln/input_raw/iso_shifts.tls) that can then be submitted to the [TALOS-N webserver](https://spin.niddk.nih.gov/bax/nmrserver/talosn/) to generate dihedral angle restraints. Note that TALOS inputs can also be generated by Sparky. Take the [pred.tab](examples/sln/input_raw/pred.tab) output file returned by the server and convert it to a restraint file using the [convertTalos](https://nmr.cit.nih.gov/xplor-nih/doc/current/helperPrograms/convertTalos.html) tool shipped with XPLOR-NIH:

	convertTalos -out iso_shifts.tbl -predFile pred.tab

The [iso_shifts.tbl](examples/sln/input_xplor/iso_shifts.tbl) will be used as the XPLOR-NIH input.




#### Step 2: Prepare CSA and DC restraints

Backbone <sup>15</sup>N CSA and <sup>1</sup>H-<sup>15</sup>N DC measurements are usually taken directly from crosspeak peak positions observed in a PISEMA (or other separated local field, SLF) OS-ssNMR spectrum of the membrane protein aligned using either bicelles or glass plates. Note that the hybrid-method.py script can be used if either CSA or DC measurements are missing (i.e., if data was obtained from 1D <sup>15</sup>N CP experiments of selectively labeled samples). CSA and DC restraints must be scaled appropriately to account for the general order parameter associated with residual dynamics of the protein and the alignment medium. To perform this scaling, the [slf2xplor.py](helpers/slf2xplor.py) helper tool is available and will output the appropriate CSA and DC restraint tables in XPLOR-NIH format.

To use this tool, first produce a TSV file of residue numbers, names, chemical shifts and dipolar couplings like [ossnmr.dat](examples/sln/input_raw/ossnmr.dat). If values are missing, enter a non-numeric value ("n", "nan" etc.). Do not leave any values blank as the file will not read correctly. For SLN, a [hcSE-SAMPI4](https://link.springer.com/article/10.1007/s10858-019-00273-1) spectrum was acquired in *unflipped* bicelles (--align_order -0.5) and fitted to a general order parameter of 0.9 (--order 0.9) using the [PISA-SPARKY](https://github.com/weberdak/pisa.py) program we have written specifically to fit SLF spectra of alpha-helical membrane proteins to the PISA model. If the spectrum was acquired using *flipped* bicelles or glass plates, then the --align_order parameter is set to 1.0 (default).  CSA and DC restraint tables are generated from the [ossnmr.dat](examples/sln/input_raw/ossnmr.dat) file by:

	python3 slf2xplor.py -i ossnmr.dat -o ossnmr --order 0.9 --align_order -0.5

This will output three XPLOR-NIH restraint tables: [ossnmr_cs.tbl](examples/sln/input_xplor/ossnmr_cs.tbl), [ossnmr_dc.tbl](examples/sln/input_xplor/ossnmr_dc.tbl) and [ossnmr_cs_gly.tbl](examples/sln/input_xplor/ossnmr_cs_gly.tbl). Note that glycine CSs are treated separately since they require a unique shift tensor. For CS restraints, the isotropic chemical shift, set as the average of the shift tensor components, are subtracted from the oriented CS observed in SLF spectra. Oriented chemical shifts MUST be [externally referenced correctly](http://dx.doi.org/10.1016/j.ssnmr.2014.03.003). The "reduced" shift is then divided by general order parameter and the order of the alignment medium (i.e., 0.9 * -0.5). DCs are scaled the same way. In the above command, default <sup>15</sup>N tensor parameters are used, which for backbone amides of non-glycine residues (--pas 57.3 81.2 228.1) are taken from  [Murray et. al., J. Mag. Res., 2014](https://doi.org/10.1016/j.jmr.2013.12.014) and from [Straus et. al., J. Biol. NMR, 2003](https://doi.org/10.1023/A:1024098123386) for glycines (--pas_gly 45.6 66.3 211.6). CSA and DC positions are assumed to have errors of 5 ppm (--error_csa 5.0) and 0.5 kHz (--error_dc 0.5), respectively. This accounts for linewidths and errors associated with assuming a constant shift tensor for all residues. To modify shift tensors and errors, these options must be explicitly stated in the command line:

	python3 slf2xplor.py \
	   	-i ossnmr.dat \
	   	-o ossnmr \
	   	--order 0.9 \
	   	--align_order -0.5 \
	   	--pas 57.3 81.2 228.1 \
	   	--pas_gly 45.6 66.3 211.6 \
	   	--error_csa 5.0 \
	   	--error_dc 0.5

The [ossnmr_cs.tbl](examples/sln/input_xplor/ossnmr_cs.tbl), [ossnmr_dc.tbl](examples/sln/input_xplor/ossnmr_dc.tbl) and [ossnmr_cs_gly.tbl](examples/sln/input_xplor/ossnmr_cs_gly.tbl) restraint files are now ready to be used as inputs for XPLOR-NIH.



#### Step 3: Prepare additional restraints for helical segments

While the backbone dihedral restraints were obtained from MAS-ssNMR and essentially encode the alpha-helical structure of the TM segment, it is also helpful to reinforce these with synthetic distance (NOE) and hydrogen bonding (HBDA) restraints. The HBDA restraints are also reinforced by the [HBDB knowledge-based potential](https://nmr.cit.nih.gov/xplor-nih/xplorMan/hbdb.html). Since we know that SLN is structured from residues R6 to S28, we can use the [synthHelix.py](synthHelix.py) helper tool:

	python3 synthHelix.py \
	   --start 6 \
	   --stop 28 \
	   --sequence MGINTRELFLNFTIVLITVILMWLLVRSYQY \
	   --start_id 1 \
	   --out_prefix sln

The --sequence and --start_id options do not need to be specified as they are generally required to scan the sequence for prolines, which SLN does not have. Three files will be written:

* [sln.dihe.tbl](examples/sln/input_xplor/sln.dihe.tbl) contains phi (-63 degrees) and psi (-42 degrees) dihedral restraints for the segment. These will be ignored in the calculation since we have experimental ones.
*  [sln.hbnoe.tbl](examples/sln/input_xplor/sln.hbnoe.tbl) includes ideal 2.06 angstrom distance restraints from the amide H to the carbonyl oxygen (i-4) atoms.
*  [sln.hbda.tbl](examples/sln/input_xplor/sln.hbda.tbl) includes the hydrogen bond restraints.



#### Step 4: Generate starting structure

The [pdbutil webserver](https://spin.niddk.nih.gov/bax/nmrserver/pdbutil/ext.html) is a good tool to generate starting structures. The structure may either be in a folded or extended state. Note that regions will be unfolded during high-temperature stages and likely remain unfolded if no restraints are applied. For this example, the [sln_seq.tbl](examples/sln/inputraw/sln_seq.tbl)  file was used as an input and the structure returned was fully helical.



#### Step 5: Running the simulated annealing calculation

Now that all of the restraint tables and input files have been prepared, we are ready to run the first simulated annealing stage. The [hybrid-method.py](hybrid-method.py) is used and is run by from customizable BASH script ([hybrid-method_local.sh](examples/hybrid-method_local.sh)) as follows:

	#!/bin/bash -l
	
	xplor -py -smp 4 -o logfile.out hybrid-method.py \
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

which is run from a Bash terminal by:

	bash hybrid-method_local.sh

The progress of the calculation can be followed using the tail command:

	tail -f logfile.out

Running this shell script allows us to avoid the having to edit the hybrid-method.py program, which requires a good understanding of both XPLOR-NIH and Python scripting. Although, the line "xplor.requireVersion('3.0')" will have to be modified if a different version of XPLOR-NIH is installed. This generic method is very well suited as quick-start/introductory approach to solving structure/topology of simple single-pass helical membrane proteins. However, if more complicated systems are to be solved, such as multi-pass proteins, beta-barrels and oligomers, the hybrid-method.py and Bash scripts will have to be edited/developed further.

It should also be emphasized that all restraints (DC_NH, CSA_N1, CSA_N1_gly, HBDA, NOE and DIHE) are completely optional. If these restraints are not available, then simply put "" in place of the file paths. For example, the following will not apply NOE or hydrogen bonding restraints:

	xplor -py -smp 4 -o logfile.out hybrid-method.py \
	  --structure_in      input_xplor/sln_ext.pdb \
	  --DC_NH             input_xplor/ossnmr_dc.tbl \
	  --CSA_N1            input_xplor/ossnmr_cs.tbl \
	  --CSA_N1_gly        input_xplor/ossnmr_cs_gly.tbl \
	  --HBDA              "" \
	  --NOE               "" \
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

The basic XPLOR-NIH protocol is taken from ([Tian et al., J. Biol. NMR, 2017](https://doi.org/10.1007/s10858-016-0082-5)) and applies EEFx force field and IMMx implicit membrane model ([Tian et al., Biophys. J., 2015](https://doi.org/10.1016/j.bpj.2015.06.047) ) developed by the Marassi Lab.

The IMMx membrane hydrophobic thickness is specified by "--immx_thickness", which for this example was set to the thickness of a DMPC/POPC (3:1 molar ratio) bicelle (25.8 angstroms) as taken from the weighted averages of DMPC and POPC (= 25.4 x 0.75 + 27.0 x 0.25). See p. 379 of the [Handbook of Lipid Bilayers - 2nd Edition](https://books.google.com/books?hl=en&lr=&id=JgnLBQAAQBAJ&oi=fnd&pg=PP1&dq=Marsh,+Handbook+of+lipid+bilayers&ots=t_vN1All4R&sig=yQdbT7LruGwNszFa586dCV1WPX4#v=onepage&q=Marsh%2C%20Handbook%20of%20lipid%20bilayers&f=false) (Marsh, 2013). The IMMx nparameter (--immx_nparameter),  which defines transition at the hydrophobic/hydrophilic interface, was also set to it default value of 10.

The maximum <sup>15</sup>N-<sup>1</sup>H dipolar coupling (--DC_NH) was set to the default 10.735 kHz ([Denny et al., J. Mag. Res, 2001](https://doi.org/10.1006/jmre.2001.2405)).

The principal axis system (PAS) of the non-glycine chemical shift tensor (--CSA_N1_tensor) was set the default  delta11 = 57.3 ppm, delta22 = 81.2 ppm, delta33 = 228.1 ppm ([Murray et al., J. Mag. Res., 2014](http://dx.doi.org/10.1016/j.jmr.2013.12.014)), and the glycine PAS (--CSA_N1_tensor_gly) to delta11 = 45.6 ppm, delta22 = 66.3 ppm, delta33 = 211.6 ppm ([Straus et al., J. Biol. NMR, 2003](https://doi.org/10.1023/A:1024098123386)). Note that the <sup>15</sup>N PAS components will be automatically converted to the XPLOR-NIH input format. The beta values were set to default values of -17.0 and -21.6 degrees as per the above references. Beta values for the <sup>15</sup>N PAS are negative due to the right-hand rule. 

The --w_slf and --w_r are the weighting terms for the DC and CSA restraints in the form typically used the Cross ([Kim et al., J. Am. Chem. Soc., 2001](https://doi.org/10.1021/ja003380x)) and Veglia ([Shi et. al., J. Biol. NMR, 2009](https://doi.org/10.1007/s10858-009-9328-9)) groups. The --w_slf of 5 specifies that the sum of the CSA and DC restraints in 5x the torsion angle (DIHE) term (i.e., 5 x 200 kcal/mol = 1000 kcal/mol) and --w_r of 3 specifies that the DC term (i.e., 750 kcal/mol/kHz) is weighted 3x more heavily than the CSA term (i.e., 250 kcal/mol/ppm). These restraints are applied as flat-well potentials in our protocols.



A basic overview of this simulated annealing protocol follows:

0. Read input PDB structure
1. Initial torsion angle minimization (100 steps)
2. Center protein to membrane (--tm_domain) then high temperature torsion dynamics with REPEL force field (3500 K for 3 ps, 3000 steps)
3. Center protein then high temperature torsions dynamics phasing in EEFx parameters (3500 K for 3 ps, 3000 steps)
4. Center protein then high temperature torsion dynamics with only EEFx parameters (3500 K for 26 ps, 26000 steps)
5. Center protein again then simulated annealing (3500 K to 25 K in 12.5 K steps, 0.2 ps/200 steps per increment).
6. Low temperature torsion dynamics (25 K for 15 ps, 15000 steps)
7. Powell torsion angle minimization (500 steps)
8. Powell Cartesian minimization (500 steps)



#### Step 6: Refinement

For the next stage, 

	#!/bin/bash -l
	
	xplor -py -smp 4 -o logfile.out hybrid-method-refine.py \
	  --structure_in      input_xplor/hybrid-method_326.sa \
	  --DC_NH             input_xplor/ossnmr_dc.tbl \
	  --CSA_N1            input_xplor/ossnmr_cs.tbl \
	  --CSA_N1_gly        input_xplor/ossnmr_cs_gly.tbl \
	  --HBDA              input_xplor/sln.hbda.tbl \
	  --NOE               input_xplor/sln.hbnoe.tbl \
	  --DIHE              input_xplor/iso_shifts.tbl \
	  --DC_NH_max         10.735 \
	  --nstructures       100 \
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

