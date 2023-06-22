cd ../../data_analysis
gmx sasa -f trajectories/mt2_rep0.xtc -s trajectories/mt2_rep0.tpr -o sasa/area.xvg -odg sasa/sfe.xvg -or  sasa/or.xvg -oa sasa/oa.xvg -tv sasa/volume.xvg -surface protein
