# hybrid-method

Calculation of membrane protein structures using isotropic and anisotropic NMR restraints.

## Example 1: Sarcolipin Structure

#### Step 1: Prepare Dihedral restraints

Prepare [iso_shifts.dat](example/sln/input_raw/iso_shifts.dat) TSV file

	python tsv2talos.py -i iso_shifts.dat -o iso_shifts.tls -s MGINTRELFLNFTIVLITVILMWLLVRSYQY

Produces TALOS-N input file [iso_shifts.tls](example/sln/input_raw/iso_shifts.tls)

Upload the iso_shifts.tls file into the [TALOS-N webserver](https://spin.niddk.nih.gov/bax/nmrserver/talosn/) to convert to dihedral angle backbone restraints. 
	
	convertTalos -out iso_shifts.tbl -predFile pred.tab
	
The [iso_shifts.tbl](example/sln/input_raw/iso_shifts.tbl) will be used as the XPLOR-NIH input.



