# Make sure the hybridTools.tcl is loaded from the correct directory
source /home/daniel/GitHub/hybrid-method/helpers/hybridTools.tcl


# Select TM residues and CA of residue (i-1) to calculation azimuthal angle
set tm [atomselect top "resid 6 to 28"]
set azi [atomselect top "resid 15 and name CA"]


# Compute tilt and azimuthal (Leu16) angles. Align models (only Z-rotations and X-Y translations).
# Note: the azi is directed at the i-1 Ca of the respective residue.
set r [helix_topology $tm $azi -invert -inside C -out results_tilt6-28_azi16.dat -no_z]


# write aligned PDBs
set s [atomselect top "all protein"]
$s set chain A
set num_frames [molinfo top get numframes]
for {set i 0} {$i < $num_frames} {incr i} {
    $s frame $i
    $s writepdb sln.model$i.pdb
}


## Draw membrane
draw delete all
draw color black
set thickness 25.72
draw material Transparent
draw cylinder [list 0 0 [expr ($thickness/2)+0.1]] [list 0 0 [expr ($thickness/2)-0.1]] radius 20 resolution 100
draw cylinder [list 0 0 [expr (-$thickness/2)+0.1]] [list 0 0 [expr (-$thickness/2)-0.1]] radius 20 resolution 100


## Find clashed (if any)
set s [atomselect top "all protein"]
clashes $s -cutoff 2.2 -allowance 0.7 -o clashes.dat
