# hybrid-method

Workflow for calculating membrane protein structures using isotropic and anisotropic NMR restraints. 

## Description

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

This method requires a UNIX operating system (i.e., MAC-OS or Linux) to install XPLOR-NIH. All commands shown below are executed from a BASH terminal.

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
*  [sln.hbda.tbl](examples/sln/input_xplor/sln.hbda.tbl) include the hydrogen bond restraints.



#### Step 4: Generate extended starting structure

The [pdbutil webserver](https://spin.niddk.nih.gov/bax/nmrserver/pdbutil/ext.html) is a good tool to generate a starting structure. Note that the structure calculation will automatically force the starting into an extended state. For this example, the [sln_seq.tbl](examples/sln/inputraw/sln_seq.tbl)  file was used as an input.



#### Step 5: Running the simulated annealing calculation

