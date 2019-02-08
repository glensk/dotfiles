parray env HOME
package require pbctools
mol load xyz "/home/glensk/Downloads/tmp2/simulation.pos_0.xyz"
pbc set {11.474929 11.474929 11.474929 60.000000 60.000000 60.000000} -all
pbc unwrap
#animate write xyz unwrapped.xyz
set fo [open "myout.dat" a]
for {set i 0} {$i < 10} {incr i} {set sel1 [atomselect top "index 63" frame $i];set out [$sel1 get {x y z}];puts $fo $out}
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

