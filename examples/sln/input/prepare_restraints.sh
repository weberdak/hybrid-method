# Run script using: bash prepare_restraints.sh

# Set DIR of helper directory
HELPERS=~/GitHub/hybrid-method/helpers


# Generate helical restraints
python3 $HELPERS/synthHelix.py \
	--start 6 \
	--stop 28 \
	--sequence AMGINTRELFLNFTIVLITVILMWLLVRSYQY \
	--start_id 0 \
	--out_prefix sln


# Process oriented restraints
python3 $HELPERS/slf2xplor.py \
	-i raw/ossnmr.dat \
	-o sln \
	--order 0.9 \
	--align_order -0.5 \
	--pas 57.3 81.2 228.1 \
	--pas_gly 45.6 66.3 211.6 \
	--error_csa 5.0 \
	--error_dc 0.5


# Create working and free CSA restraints
python3 $HELPERS/rfree.py \
	-i sln_cs.tbl sln_cs_gly.tbl \
	-r 20.0

# Create working and free DC restraints
python3 $HELPERS/rfree.py \
	-i sln_dc.tbl \
	-r 20.0


