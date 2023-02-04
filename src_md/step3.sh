#!/bin/bash

for run in 8DFN_full 8DFN_half 7SI9_full 7SI9_half
do
    cd ../$run

    mkdir 5_POST wetness_analysis pyinteraph
    cp 3_NPT/npt.gro 5_POST/md.gro

    ############################################################################ POST-PROCESSING
    echo 1 0 | gmx trjconv -s 4_MD/md_plain.tpr -f 4_MD/md_plain.xtc -o 5_POST/md-0_nojump.xtc -pbc nojump -center

    echo 3 0 | gmx trjconv -s 4_MD/md_plain.tpr -f 5_POST/md-0_nojump.xtc -o 5_POST/md-1_center.xtc -pbc mol -center
    rm 5_POST/md-0_nojump.xtc

    echo 1 0 | gmx trjconv -s 4_MD/md_plain.tpr -f 5_POST/md-1_center.xtc -o 5_POST/2_md-rottrans.xtc -fit rot+trans
    rm 5_POST/md-1_center.xtc

done
