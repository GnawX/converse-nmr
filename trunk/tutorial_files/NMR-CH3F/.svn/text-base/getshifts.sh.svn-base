#!/bin/sh

echo "H shielding:"
grep Chemical *-H*out | awk ' { tr=$5; getline; tr+=$6; getline; tr+=$7; print tr/3.0}'

echo

echo "C shielding:"
grep Chemical *-C*out | awk ' { tr=$5; getline; tr+=$6; getline; tr+=$7; print 200.5108+tr/3.0}'

echo

echo "F shielding:"
grep Chemical *-F*out | awk ' { tr=$5; getline; tr+=$6; getline; tr+=$7; print 306.4172+tr/3.0}'
