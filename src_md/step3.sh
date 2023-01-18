#!/bin/bash

for run in 8DFN_full 8DFN_half 7SI9_full 7SI9_half
do
    cd ../$run

    mkdir 5_POST
    cp 3_NPT/npt.gro 5_POST/md.gro

    ############################################################################ POST-PROCESSING
    echo 1 0 | gmx trjconv -s 4_MD/md_plain.tpr -f 4_MD/md_plain.xtc -o 5_POST/md-nojump.xtc -pbc nojump -center

    echo 1 0 | gmx trjconv -s 4_MD/md_plain.tpr -f 5_POST/md-nojump.xtc -o 5_POST/md-rottrans.xtc -fit rot+trans
    rm 5_POST/md-nojump.xtc

    echo 1 0 | gmx trjconv -s 4_MD/md_plain.tpr -f 5_POST/md-rottrans.xtc -o 5_POST/md-center.xtc -pbc mol -center
    rm 5_POST/md-rottrans.xtc

done
