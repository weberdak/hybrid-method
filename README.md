# hybrid-method

Calculation of membrane protein structures using isotropic and anisotropic NMR restraints. 

## Description

This repository details a protocol to determine the STRUCTURE AND TOPOLOGY of membrane proteins by hybridizing isotropic and anisotropic NMR restraints into XPLOR-NIH simulated annealing calculations. The protocol has been significantly updated from our previous versions ([Shi et. al., J. Biol. NMR, 2009](https://doi.org/10.1007/s10858-009-9328-9)) in incorportate new functions and force fields that have since been built into XPLOR-NIH. This protocol has also been designed for users who are not expert users of XPLOR-NIH and includes the following features:

* Incorporates the new [EEFx forcefield]() and [IMMx implicit membrane model]() built into XPLOR-NIH by the Marassi Lab.
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

An adequate set of backbone dihedral restraints are usually neccesary to obtain the correct fold during simulated annealing. The hybrid-method.py script will function without these restraints, but calculations will gnerally run poorly in our experience. If chemical shifts are available, copy them into the TSV formatted file like [iso_shifts.dat](examples/sln/input_raw/iso_shifts.dat). These can be prepared using a spreadsheet then copying and pasting into a text file. Make sure the header line is correct and only standard atom names are used. Convert this file to a TALOS-N input file using the [tsv2talos.py](helpers/tsv2talos.py) script as follows:

	python tsv2talos.py -i iso_shifts.dat -o iso_shifts.tls -s MGINTRELFLNFTIVLITVILMWLLVRSYQY

This produces the [iso_shifts.tls](examples/sln/input_raw/iso_shifts.tls) that can then be submitted to the [TALOS-N webserver](https://spin.niddk.nih.gov/bax/nmrserver/talosn/) to generate dihedral angle restraints. Note that TALOS input can also be generated by Sparky. 
	
	convertTalos -out iso_shifts.tbl -predFile pred.tab
	
The [iso_shifts.tbl](examples/sln/input_xplor/iso_shifts.tbl) will be used as the XPLOR-NIH input.


#### Step 2: Prepare CSA and DC restraints

Build a TSV file [ossnmr.dat](examples/sln/input_raw/ossnmr.dat).

	python slf2xplor.py -i ossnmr.dat -o ossnmr --order 0.9 --align_order -0.5

This will output three XPLOR-NIH files: [ossnmr_cs.tbl](examples/sln/input_xplor/ossnmr_cs.tbl), [ossnmr_dc.tbl](examples/sln/input_xplor/ossnmr_dc.tbl) and [ossnmr_cs_gly.tbl](examples/sln/input_xplor/ossnmr_cs_gly.tbl).

	python slf2xplor.py \
       	-i ossnmr.dat \
       	-o ossnmr \
       	--order 0.9 \
       	--align_order -0.5 \
       	--pas 57.3 81.2 228.1 \
       	--pas_gly 45.6 66.3 211.6 \
       	--error_csa 3.0 \
       	--error_dc 0.3
  
  Default parameters for backbone amides are taken from [Murray et. al., J. Mag. Res., 2014](https://doi.org/10.1016/j.jmr.2013.12.014) for non-glycine residues and [Straus et. al., J. Biol. NMR, 2003](https://doi.org/10.1023/A:1024098123386) for glycines
  


  