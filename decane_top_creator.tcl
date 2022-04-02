#!/usr/bin/tclsh

package require topotools

package require pbctools

mol new decane_box_void.xyz

#H20

set selh [atomselect top {name H}]
$selh set type H
$selh set mass 1.00794
$selh set charge 0.045

set selc1 [atomselect top {name C}]
$selc1 set type C
$selc1 set mass 12.0107
$selc1 set charge -0.09

set selc2 [atomselect top {name C2}]
$selc2 set type C
$selc2 set mass 12.0107
$selc2 set charge -0.135

#mol bondsrecalc top
topo retypebonds
topo guessangles
topo guessdihedrals

mol reanalyze top

topo writelammpsdata decane_pre.top

quit
