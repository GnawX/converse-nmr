#!/bin/sh

pw=../../PW/pw.x

molecule=H2CN
cart=(dummy x y z)

for dir in 1 2 3
do
  base=$molecule-converse-${cart[$dir]}
  sed "s/@DIR@/$dir/" <$molecule-converse.in.tmpl >$base.in
  $pw <$base.in >$base.out
done

