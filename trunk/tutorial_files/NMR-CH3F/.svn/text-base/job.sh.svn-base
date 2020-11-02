#!/bin/sh

pw=../../PW/pw.x

molecule=CH3F
natoms=5
atoms=(null C H H H F)
cart=(null x y z)

for i in `seq 1 $natoms`
do
  for dir in 1 2 3
  do
    base=$molecule-converse-${atoms[$i]}$i${cart[$dir]}
    sed "s/@ATOM@/$i/;s/@DIR@/$dir/" <$molecule-converse.in.tmpl >$base.in
    $pw <$base.in >$base.out
  done
done

