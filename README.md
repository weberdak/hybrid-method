# hybrid-method

Calculation of membrane protein structures using isotropic and anisotropic NMR restraints. 

## Description

This repository details a protocol to determine the STRUCTURE AND TOPOLOGY of membrane proteins by hybridizing isotropic and anisotropic NMR restraints into XPLOR-NIH simulated annealing calculations. The protocol has been significantly updated from our previous versions ([Shi et. al., J. Biol. NMR, 2009](https://doi.org/10.1007/s10858-009-9328-9)) in incorportate new functions and force fields that have since been built into XPLOR-NIH. This protocol has also been designed for users who are not expert users of XPLOR-NIH and includes the following features:

* Incorporates the new [EEFx force field]() and [IMMx implicit membrane model]() built into XPLOR-NIH by the Marassi Lab.
* A [hybrid-method.py](hybrid-method.py) script containing a preset protocol and parameter set. This script is never modified and is instead run using a BASH configuration file. This allows users to perform calculations without an advanced knowledge of XPLOR-NIH. Although, users should eventually take the time to review the hybrid-method.py file for a better understanding of what it does.
* Supports CSA and DC restraints obtained by OS-ssNMR using bicelles (flipped and unflipped).
* [tsv2talos.py](helpers/tsv2talos.py) helper tool to convert a TSV file of chemical shifts into TALOS-N input format.
* [slf2xplor.py](helpers/slf2xplor.py) helper tool to prepare XPLOR-NIH restraint files for CSs and DCs. This file applies appropriate scaling corrections to account for bicelle alignment orientation and protein dynamics.
* [hybridTools.tcl](helpers/hybridTools.tcl) library of VMD functions to analyze results. Including:
	* Alignment of helical segments so topology is unaffected (i.e., not applying rotations)
	* Helix tilt and azimuthal angle measurements
* All restraining potentials are optional. Just leave option blank in the BASH configuration script if data is unavailable.

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

This method requires a UNIX operating system (i.e., MAC-OS or Linux) to install XPLOR-NIH. All commands shown below are excuted from a BASH terminal.

* Python 3 (must include numpy)
* [XPLOR-NIH 3.0](https://nmr.cit.nih.gov/xplor-nih/)
* [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) (for analysis)

To download, click the "Clone or Download" button above and "Download ZIP". Scripts are usually copied into the directory they are executed in.

## Examples

### Sarcolipin Structure

Sarcolipin (SLN) is a single-pass membrane protein and a well-known regulator of the sarco/endoplasmic reticulum Ca2+-ATPase (SERCA) pump responsible for the relaxation of skeletal muscle cells post-contraction. SLN is a simple protein and a good starting point for learning how to hybridize isotropic and anisotropic NMR restraints into simulated annealing structure calculations using XPLOR-NIH.

#### Step 1: Prepare dihedral restraints

An adequate set of backbone dihedral restraints are usually neccesary to obtain the correct fold during simulated annealing. The hybrid-method.py script will function without these restraints, but calculations will generally run poorly in our experience. If chemical shifts are available, copy them into the TSV formatted file like [iso_shifts.dat](examples/sln/input_raw/iso_shifts.dat). These can be prepared using a spreadsheet, then copyied and pasted into a text file. Make sure the header line is correct and only standard atom names are used. Convert this file to a TALOS-N input file using the [tsv2talos.py](helpers/tsv2talos.py) helper script as follows:

	python tsv2talos.py -i iso_shifts.dat -o iso_shifts.tls -s MGINTRELFLNFTIVLITVILMWLLVRSYQY

This produces the [iso_shifts.tls](examples/sln/input_raw/iso_shifts.tls) that can then be submitted to the [TALOS-N webserver](https://spin.niddk.nih.gov/bax/nmrserver/talosn/) to generate dihedral angle restraints. Note that TALOS inputs can also be generated by Sparky. Take the [pred.tab](examples/sln/input_raw/pred.tab) output file returned by the server and convert it to a restraint file using the [convertTalos](https://nmr.cit.nih.gov/xplor-nih/doc/current/helperPrograms/convertTalos.html) tool shipped with XPLOR-NIH:
	
	convertTalos -out iso_shifts.tbl -predFile pred.tab
	
The [iso_shifts.tbl](examples/sln/input_xplor/iso_shifts.tbl) will be used as the XPLOR-NIH input.


#### Step 2: Prepare CSA and DC restraints

15N CSA and 1H-15N DC measurements are usually taken directly from crosspeak peak positions observed in a PISEMA (or other seperated local field, SLF) OS-ssNMR spectrum of the membrane protein aligned using either bicelles or glass plates. Note that the hybrid-method.py script can be used if either CSA or DC measurements are missing (i.e., if data is only available from 1D 15N CP experiments of selectively labelled samples). CSA and DC restraints must be scaled appropriately to account for residual dynamics (general order parameter) and the orientation of alignment (bicelle flip angle). To perform this scaling, and then output restraints in XPLOR-NIH format, the [slf2xplor.py](helpers/slf2xplor.py) helper tool is available.

First, produce a TSV file of residue numbers, name, chemical shifts and dipolar coupling like [ossnmr.dat](examples/sln/input_raw/ossnmr.dat). If values are unavailable, specify a non-numeric value ("n", "nan" etc.) - do not leave any values blank as the file will not read correctly. For SLN, a [hcSE-SAMPI4](https://link.springer.com/article/10.1007/s10858-019-00273-1) spectrum was acquired in unflipped bicelles (--align_order -0.5) and fitted to a general order parameter of 0.9 (--order 0.9) using the [PISA-SPARKY](https://github.com/weberdak/pisa.py) program we have written to analyse SLF spectra of alpha-helical membrane proteins. CSA and DC restraint tables are generated from [ossnmr.dat](examples/sln/input_raw/ossnmr.dat) by:

	python slf2xplor.py -i ossnmr.dat -o ossnmr --order 0.9 --align_order -0.5

This will output three XPLOR-NIH restraint tables: [ossnmr_cs.tbl](examples/sln/input_xplor/ossnmr_cs.tbl), [ossnmr_dc.tbl](examples/sln/input_xplor/ossnmr_dc.tbl) and [ossnmr_cs_gly.tbl](examples/sln/input_xplor/ossnmr_cs_gly.tbl). Note that glycines CS are treated seperately since they require a unique shft tensor. For CS restraints, the isotropic chemical shift, determined from the average of the shift tensor compents, are subtratracted from the oriented CS observed in SLF spectra. Oriented chemical shifts MUST be [externally referenced correctly](http://dx.doi.org/10.1016/j.ssnmr.2014.03.003). The "reduced" shift is then divided by dynamic and alignment order parameters (i.e., 0.9 * -0.5). DC value a scaled the same way. In the above command, default 15N tensor parameters are used for backbone amides [Murray et. al., J. Mag. Res., 2014](https://doi.org/10.1016/j.jmr.2013.12.014) and non-glycine residues (--pas 57.3 81.2 228.1) and [Straus et. al., J. Biol. NMR, 2003](https://doi.org/10.1023/A:1024098123386) for glycines (--pas_gly 45.6 66.3 211.6). CSA and DC positions are assumed to have errors of 3 ppm (--error_csa 3.0) and 0.3 kHz (--error_dc 0.3), respectively. This accounts for linewidths and errors accociated with assuming a constant shift tensore for all residues. To modify the shift tensors and errrors, the options must be explicitly stated in the command line:

	python slf2xplor.py \
       	-i ossnmr.dat \
       	-o ossnmr \
       	--order 0.9 \
       	--align_order -0.5 \
       	--pas 57.3 81.2 228.1 \
       	--pas_gly 45.6 66.3 211.6 \
       	--error_csa 3.0 \
       	--error_dc 0.3
  
Now that we have [ossnmr_cs.tbl](examples/sln/input_xplor/ossnmr_cs.tbl), [ossnmr_dc.tbl](examples/sln/input_xplor/ossnmr_dc.tbl) and [ossnmr_cs_gly.tbl](examples/sln/input_xplor/ossnmr_cs_gly.tbl) restraint files, we are ready to do the structure calculation.
  
#### Step 3: Running the simulated annealing calculation



  