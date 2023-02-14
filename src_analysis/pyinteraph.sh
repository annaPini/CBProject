cd ../data_analysis/pyinteraph

PATH_GRO="../trajectories/mt2_rep0.gro"
PATH_XTC="../trajectories/mt2_rep0.xtc"
PATH_TPR="../trajectories/mt2_rep0.tpr"

###################
### salt bridges
pyinteraph -s $PATH_TPR -t $PATH_XTC -r $PATH_GRO --sb-co 5 -b --sb-graph mt2_rep0-sb_graph.dat --ff-masses charmm27 -v --sb-cg-file charged_groups.ini

### H bonds
# pyinteraph -s $PATH_TPR -t $PATH_XTC -r $PATH_GRO -y --hb-graph mt2_rep0-hb_graph.dat --ff-masses charmm27 -v --hb-ad-file hydrogen_bonds.ini

### Hydrophobic interaction
pyinteraph -s $PATH_TPR -t $PATH_XTC -r $PATH_GRO -f --hc-co 5 -f --hc-graph mt2_rep0-hc_graph.dat --ff-masses charmm27 -v --hc-residues ALA,VAL,LEU,ILE,PHE,PRO,MET,TRP

###################
### cutoff for occurrence percentage
filter_graph -d mt2_rep0-sb_graph_all.dat -c mt2_rep0-clusters_size_sb.dat -p mt2_rep0-clusters_plot_sb.pdf
# filter_graph -d mt2_rep0-hb_graph_all.dat -c mt2_rep0-clusters_size_hb.dat -p CLUSTER_PLOT_HB
filter_graph -d mt2_rep0-hc_graph_all.dat -c mt2_rep0-clusters_size_hc.dat -p mt2_rep0-clusters_plot_hc.pdf

### filtering the graphs using the identified cutoff of the occurrence percentage
filter_graph -d mt2_rep0-sb_graph_all.dat -o mt2_rep0-sb_graph_filtered.dat -t 20.0
# filter_graph -d mt2_rep0-hb_graph_all.dat -o mt2_rep0-hb_graph_filtered.dat -t 20.0
filter_graph -d mt2_rep0-hc_graph_all.dat -o mt2_rep0-hc_graph_filtered.dat -t 20.0

###################
