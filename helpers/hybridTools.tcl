# hybridTools.tcl
# ---------------
# A VMD function to align and analyze helix topology.
# Written by D. K. Weber (last revised June 7 2020)
#
# Features
# --------
# - Align helix by center of mass without affecting tilt or rotation angles
# - Option to not apply Z translation (if depth of insertion is important)
# - Align helical segment along X axis (good for figures)
# - Measure tilt and azimuthal angles

proc helix_topology { h_selection a_selection { args } } {
    # Centers helix selection to origin. Aligns along Z. Measures tilt and azimuthal angle.
    #
    # Parameters
    # ----------
    # h_selection: obj
    #    VMD selection of helical segment. The code below will automatically extract relevant info.
    # a_selection: obj
    #    VMD selection of residue for which azimuthal angle is calculated (i.e., a CA atom usually).
    #    If more than one atoms is selected, the angle will be determined from the center of mass.
    #
    # Optional flags
    # --------------
    # -fix_tilt <angle>
    #    <angle>: float
    #       Rotate helix in each frame so that tilt angle is fixed to this value. Default: null.
    # -invert
    #    Helix reads from C to N terminus. Default: off.
    # -inside <term>
    #    <term>: N or C
    #       The helix is flipped so that <term> is on lower leaflet. Default: N.
    # -no_z
    #    Do not apply a Z translation when aligning helices. Default: off.
    # -out: <file name>
    #    <file name>: Name of file to output tilt and azimuthal angles. Default: null.

    # Handle optional arguments
    set fix_tilt [argparse $args "-fix_tilt" 1 "null"]
    set invert [ switch_on $args "-invert" off ]
    set inside [argparse $args "-inside" 1 "N"]
    set no_z [ switch_on $args "-no_z" off ]
    set out [argparse $args "-out" 1 "null"]
    
    # Initialize sub-selections
    set indices [lsort -integer [$h_selection get index]]
    set half [expr int([llength $indices]/2)]
    set i 0
    set firsthalf []
    set secondhalf []
    foreach index $indices {
	if {$i <= $half} { lappend firsthalf $index }
	if {$i >= $half} { lappend secondhalf $index }
	incr i
    }
    set indices [$h_selection get index]
    set s1 [atomselect top "index $firsthalf and name N C CA"]
    set s2 [atomselect top "index $secondhalf and name N C CA"]
    set sa [atomselect top "index $indices and name N C CA"]
    set sel_all [ atomselect top all ]

    # Select region surrounding azi residue for local centroid
    set frag [$a_selection get fragment]
    set resid [$a_selection get resid]
    set azi_region [atomselect top "fragment $frag and resid [expr $resid - 2] to [expr $resid +2] and name N C CA"]
    
    # Statistics
    set tilt_angles []
    set azi_angles []
    # Loop through each frame and record helical angle
    set num_frames [molinfo top get numframes]
    for {set i 0} {$i < $num_frames} {incr i} {
	$s1 frame $i
	$s2 frame $i
	$sa frame $i
	$sel_all frame $i
	
	# Measure center of masses for N-term (s1) and C-term (s2) helical segments
	set c1 [measure center $s1]
	set c2 [measure center $s2]

	# Flip helix so that <term> in on lower leaflet (if not already)
	if { $inside == "N" } {
	    if { [lindex $c1 2] > [lindex $c2 2] } {
		set rot [trans axis x 180 deg]
		$sel_all move $rot
	    }
	}
	if { $inside == "C" } {
	    if { [lindex $c1 2] < [lindex $c2 2] } {
		set rot [trans axis x 180 deg]
		$sel_all move $rot
	    }
	}
	
	# Re-measure centers (apply inversion if reading tilt in C-N direction)
	if { $invert == on } {
	    #puts "Helix tilt angle will be measured in C to N direction"
	    set c1 [measure center $s2]
	    set c2 [measure center $s1]
	} else {
	    #puts "Helix tilt angle will be measured in N to C direction"
	    set c1 [measure center $s1]
	    set c2 [measure center $s2]	    
	}
	
	# Measure tilt angle between center of masses
	set vt [vecnorm [vecsub $c2 $c1]]
	set t_ang [expr acos([lindex $vt 2]) * 57.2957795130823]
	#puts $t_ang
	lappend tilt_angles $t_ang

	# Align helix along X axis
	set ang [expr atan([lindex $vt 1]/[lindex $vt 0])]
	
	# Rotate so X component is always positive
	if { [lindex $vt 0] < 0 } {
	    set rot [trans axis z [expr -$ang + 3.1415926535897931] rad]
	} else {
	    set rot [trans axis z -$ang rad]
	}
	$sel_all move $rot

	# Rotate to specified fixed angle if desired
	if { $fix_tilt != "null" } {
	    set rot [trans axis y [expr $fix_tilt - $t_ang] deg]
	    $sel_all move $rot
	}
	
	# Translate center to origin
	set ca [measure center $sa]
	set trans [trans origin $ca]
	if { $no_z == on } {
	    #puts "test"
	    set xrow [lindex $trans 0]
	    set yrow [lindex $trans 1]
	    set zrow [lindex $trans 2]
	    set zrow [list [lindex $zrow 0] [lindex $zrow 1] [lindex $zrow 2] 0.0 ]
	    set irow [lindex $trans 3]
	    set trans [list $xrow $yrow $zrow $irow]
	}
	#puts $trans
	$sel_all move $trans
	#draw sphere [measure center $sa]
	
	# Calculate azimuthal angle of residue
	# Rotate tilt to zero and determine angle from X,Y displacement of azi point from centroid.
	$a_selection frame $i
	$azi_region frame $i
	set azi_point [measure center $a_selection]
	set c [measure center $azi_region]
	set rot [trans axis y [expr 0 - $t_ang] deg]
	# Rotate points along Z axis and meaure X and Y deviations
	set c [vectrans $rot $c]
	set azi_point [vectrans $rot $azi_point]
	#draw color blue
	#draw sphere $c
	#draw color red
	#draw sphere $azi_point
	set diff [vecsub $azi_point $c]
	set diffx [lindex $diff 0]
	set diffy [lindex $diff 1]
	set aziy [expr atan($diffx/$diffy) * 57.2957795130823]
	set azix [expr atan($diffy/$diffx) * 57.2957795130823]
	
	
	# Determine phase (Y = 180, -Y = 0, X = 90, -X = 270)
	# +X (0 to 180 degrees)
	if { $diffx > 0 } {
	    # +X and -Y (0 to 90 degrees)
	    if { $diffy < 0 } { set azi [expr 0 - $aziy]; #puts "+X-Y" }
	    # +X and +Y (90 to 180 degrees)
	    if { $diffy > 0 } { set azi [expr 180 - $aziy]; #puts "+X+Y"}
	}
	# -X (180 to 360 degrees)
	if { $diffx < 0 } {
	    # -X and +Y (180 to 270 degrees)
	    if { $diffy > 0 } { set azi [expr 180 - $aziy]; #puts "-X+Y"}
	    # -X and -Y (270 to 360 degrees)
	    if { $diffy < 0 } { set azi [expr 360 - $aziy]; #puts "-X-Y" }
	}   
	lappend azi_angles $azi
	#puts "$aziy\t$azix\t$azi\t$diffx\t$diffy"	
    }
    
    # Output to terminal
    set t_ang_mean [format "%.2f" [lindex [stats $tilt_angles] 0]]
    set t_ang_std  [format "%.2f" [lindex [stats $tilt_angles] 1]]
    set a_ang_mean [format "%.2f" [lindex [stats $azi_angles] 0]]
    set a_ang_std  [format "%.2f" [lindex [stats $azi_angles] 1]]
    puts "Average tilt angle of segment: $t_ang_mean +/- $t_ang_std"
    puts "Average azi angle residue [lindex $resid 0]: $a_ang_mean +/- $a_ang_std"

    # Output results to files
    if { $out != "null" } {
	puts "Writing results to $out..."
	set outf [open $out w]
	puts $outf "#Frame\tTilt\tazi_[lindex $resid 0]"
	for {set i 0} {$i < $num_frames} {incr i} {
	    set tilt [format "%.2f" [lindex $tilt_angles $i]]
	    set azi [format "%.2f" [lindex $azi_angles $i]]
	    puts $outf "$i\t$tilt\t$azi"
	}
	close $outf
    }
    
    # Clean up
    $s1 delete
    $s2 delete
    $sa delete
    $sel_all delete
    $azi_region delete
    unset s1 s2 sa sel_all azi_region
}


proc argparse { args flag input_index default } {
    # Crude argument parser that does the trick.
    #
    # Parameters
    # ----------
    # args: list of str and values
    #    List of arguments and values
    # flag: str
    #    <args> are searched for the presence of this flag
    # input_index: int
    #    Position that value occurs after <flag>. I.e., if "1", the
    #    value immediately following <flag> will be returned.
    # default: anything
    #    If <flag> isn't found in <args>, then this value is returned.
    #
    # Returns
    # -------
    # value: anything
    #    The value parsed in from <args> of default.
   
    set value $default
    if { [ lsearch $args $flag ] != -1 } { 
	set value [lindex $args [expr ([lsearch $args $flag ] + $input_index)]]
    } 
    return $value
}


# Proc to handle switched
proc switch_on { args flag default } {
    set value $default
    if { [ lsearch $args $flag ] != -1 } { 
	set value on
    }
    return $value
}


proc stats { data } {
    # Compute mean and standard deviation of list
    #
    # Parameters
    # ----------
    # data: list of floats
    #    List of numbers.
    #
    # Returns
    # -------
    # [mean, std]: list of two floats
    #    Mean and standard deviation.
    
    # Comute mean
    set csum 0
    set csum2 0
    set num [llength $data]
    foreach i $data {
	set csum [expr $csum + $i]
	set csum2 [expr $csum2 + pow($i,2)] 
    }
    set mean [expr $csum / $num]
    set std [expr sqrt((($num*$csum2)-pow($csum,2))/($num*($num-1)))]
    return [list $mean $std]
}

