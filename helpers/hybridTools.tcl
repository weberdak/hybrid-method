

proc align_helix { selection { args } } {
    # Centers helix selection to origin  
    #
    # Parameters
    # ----------
    # selection: obj
    #    VMD selection. The code below will automatically extract relevant info.
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
    #
    # Outputs
    # -------
    # Aligns helical segment 

    # Handle optional arguments
    set fix_tilt [argparse $args "-fix_tilt" 1 "null"]
    set invert [ switch_on $args "-invert" off ]
    set inside [argparse $args "-inside" 1 "N"]
    
    # Initialize sub-selections
    set indices [lsort -integer [$selection get index]]
    set half [expr int([llength $indices]/2)]
    set i 0
    set firsthalf []
    set secondhalf []
    foreach index $indices {
	if {$i <= $half} { lappend firsthalf $index }
	if {$i >= $half} { lappend secondhalf $index }
	incr i
    }
    set indices [$selection get index]
    set s1 [atomselect top "index $firsthalf and name N C CA"]
    set s2 [atomselect top "index $secondhalf and name N C CA"]
    set sa [atomselect top "index $indices and name N C CA"]
    set sel_all [ atomselect top all ]

    # Statistics
    set tilt_angles []
    
    # Loop through each frame and record helical angle
    set num_frames [molinfo top get numframes]
    for {set i 0} {$i < $num_frames} {incr i} {
	$s1 frame $i
	$s2 frame $i
	$sa frame $i
	$sel_all frame $i
	
	# Measure center
	set c1 [measure center $s1]
	set c2 [measure center $s2]

	# Flip helix so that <term> in on lower leaflet
	if { $inside == "N" } {
	    if { [lindex $c1 2] > [lindex $c2 2] } {
		puts "bamn"
		set rot [trans axis x 180 deg]
		$sel_all move $rot
	    }
	}
	if { $inside == "C" } {
	    if { [lindex $c1 2] < [lindex $c2 2] } {
		puts "bamc"
		set rot [trans axis x 180 deg]
		$sel_all move $rot
	    }
	}
	
	# Re-measure centers
	set c1 [measure center $s1]
	set c2 [measure center $s2]	    
	
	# Measure angle 
        set vs [vecnorm [vecsub $c2 $c1]]
	set ang [expr atan([lindex $vs 0]/[lindex $vs 1])]
	
	# Rotate so Y component is always positive
	if { [lindex $vs 1] < 0 } {
	    set rot [trans axis z [expr $ang + 3.1415926535897931] rad]
	} else {
	    set rot [trans axis z $ang rad]
	}
	$sel_all move $rot

	# Measure helix angle
	if { $invert == on } {
	    set vt [vecnorm [vecsub $c1 $c2]]
	} else {
	    set vt [vecnorm [vecsub $c2 $c1]]
	}
	set t_ang [expr acos([lindex $vt 2]) * 57.2957795130823]
	puts $t_ang
	lappend tilt_angles $t_ang
	
	# Rotate to specified angle if desired
	if { $fix_tilt != "null" } {
	    set rot [trans axis x [expr $fix_tilt - $t_ang] deg]
	    $sel_all move $rot
	}
	
	# Translate center to origin
	set ca [measure center $sa]
	set trans [trans origin $ca]
	$sel_all move $trans
	
    }
    # Output
    set t_ang_mean [average $tilt_angles] 
    puts "Average tilt angle of segment: $t_ang_mean"
    
    # Clean up
    $s1 delete
    $s2 delete
    $sa delete
    $sel_all delete
    unset s1 s2 sa sel_all
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


proc average { data } {
    set csum 0
    foreach i $data {
	set csum [expr $csum + $i]
    }
    return [expr $csum / [llength $data]]
}
