#!/bin/sh

pw=$HOME/Codes/converse-nmr/trunk/PW/pw.x

molecule=quartz
natoms=9
atoms=(null Si Si Si O O O O O O)
cart=(null x y z)

for i in `seq 1 $natoms`
do
  for dir in 1 2 3
  do
    base=$molecule-converse-${atoms[$i]}$i${cart[$dir]}
    sed "s/@ATOM@/$i/;s/@DIR@/$dir/" <$molecule-converse.in.tmpl >$base.in
    mpirun $pw <$base.in >$base.out
  done
done

