parray env HOME
package require pbctools
mol load xyz "/home/glensk/Downloads/tmp2/vacancy_pos_angstrom.xyz"
pbc set {11.474929 11.474929 11.474929 60.000000 60.000000 60.000000} -all
set id [molinfo top get id ]
#mol waitfor all molid $id;   # when at some point not all points are read in
set allmyframes [molinfo $id get numframes]
puts "--> frames $allmyframes"
pbc unwrap
#animate write xyz unwrapped.xyz
#set sel1 [atomselect top "index 63" frame 0];
set sel1 [atomselect top "index 0" frame 0];
set x [$sel1 get x ];
set y [$sel1 get y ];
set z [$sel1 get z ];
set r02 [expr $x*$x+$y*$y+$z*$z];

set fo [open "vacancy_radius_squared.dat" a]
for {set i 0} {$i < $allmyframes} {incr i} {
    #set sel1 [atomselect top "index 63" frame $i];
    set sel1 [atomselect top "index 0" frame $i];
    set xyz [$sel1 get {x y z} ];
    set x [$sel1 get x ];
    set y [$sel1 get y ];
    set z [$sel1 get z ];
    #set xx [expr $x*$x];
    #set yy [expr $y*$y];
    #set zz [expr $z*$z];
    set rr [expr $x*$x+$y*$y+$z*$z-$r02];

    #puts $fo "$x\t$y\t$z"
    #puts $fo "$xx\t$y\t$z"
    puts $fo "$rr"
}
close $fo
exit
# atomselect top "index 1" frame all
#
# works:
# atomselect top "index 1" frame all
# atomselect311 get {x y z}
# atomselect top "index 2" frame all
# atomselect312 get {x y z}
#
# for {set i 1} {$i < 6} {incr i} {puts $i }  # loop in tcs
#
# set sel1 [atomselect top "index 1" frame all]
# sel1 get {x y z}
# set sel1 [atomselect top "index 2" frame all]
# sel1 get {x y z}
#
# set sel1 [atomselect top "index 1" frame all];$sel1 get {x y z}
# set sel1 [atomselect top "index 2" frame all];$sel1 get {x y z}
# 
# WORKS!!!
# for {set i 1} {$i < 6} {incr i} {set sel1 [atomselect top "index $i" frame all];set out [$sel1 get {x y z}];puts $out}
# 
# THIS IS WHAT I NEED:
# for {set i 0} {$i < 6} {incr i} {set sel1 [atomselect top "index 63" frame
# $i];set out [$sel1 get {x y z}];puts $out}

