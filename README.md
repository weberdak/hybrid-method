# hybrid-method

Calculation of membrane protein structures using isotropic and anisotropic NMR restraints.

## Description

This repository details a protocol to determine the STRUCTURE AND TOPOLOGY of membrane proteins by hybridizing isotropic and anisotropic NMR restraints into XPLOR-NIH simulated annealing calculations.

* Incorporates the new [EEFx forcefield]() and [IMMx implicit membrane model]() built into XPLOR-NIH by the Marassi Lab.
* A [hybrid-method.py](hybrid-method.py) script containing a preset protocol and parameter set. This script is never modified and is instead run using a BASH configuration file. This allows users to perform calculations without an advanced knowledge of XPLOR-NIH. Although, users should eventually take the time to review the hybrid-method.py file for a better understanding of what it does.
* Supports CSA and DC restraints obtained by OS-ssNMR using bicelles (flipped and unflipped).
* [tsv2talos.py](helpers/tsv2talos.py) helper tool to convert a TSV file of chemical shifts into TALOS-N input format.
* [slf2xplor.py](helpers/slf2xplor.py) helper tool to prepare XPLOR-NIH restraint files for CSs and DCs. This file applies appropriate scaling corrections to account for bicelle alignment orientation and protein dynamics.
* [hybridTools.tcl](helpers/hybridTools.tcl) library of VMD functions to analyze results. Including:
	* Alignment of helical segments so topology is unaffected (i.e., not applying rotations)
	* Helix tilt and azimuthal angle measurements



## Installation

This method requires a UNIX operating system (i.e., MAC-OS or Linux) to install XPLOR-NIH. All commands shown below are excuted from a BASH terminal.

* Python 3 (must include numpy)
* [XPLOR-NIH 3.0](https://nmr.cit.nih.gov/xplor-nih/)
* [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) (for analysis)

To download 

## Examples

### Sarcolipin Structure

Sarcolipin (SLN) is a single-pass membrane protein and a well-known regulator of the sarco/endoplasmic reticulum Ca2+-ATPase (SERCA) pump responsible for the relaxation of skeletal muscle cells post-contraction. SLN is a simple protein and a good starting point for learning how to hybridize isotropic and anisotropic NMR restraints into simulated annealing structure calculations using XPLOR-NIH.

#### Step 1: Prepare dihedral restraints

An adequate set of dihedral restraints are neccesary to obtain the correct fold 

Prepare [iso_shifts.dat](examples/sln/input_raw/iso_shifts.dat) TSV file

	python tsv2talos.py -i iso_shifts.dat -o iso_shifts.tls -s MGINTRELFLNFTIVLITVILMWLLVRSYQY

Produces TALOS-N input file [iso_shifts.tls](examples/sln/input_raw/iso_shifts.tls)

Upload the iso_shifts.tls file into the [TALOS-N webserver](https://spin.niddk.nih.gov/bax/nmrserver/talosn/) to convert to dihedral angle backbone restraints. 
	
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
  


  