from parameters import *
import matplotlib.pyplot as plt
import os
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj

if __name__ == "__main__":
    
    python3 -m virtualenv pyinteraph_venv

    for run in RUNS:

        PATH_GRO    = DIR_DA_TRAJECTORIES / f"{run}.gro"
        PATH_XTC    = DIR_DA_TRAJECTORIES / f"{run}.xtc"
        PATH_TPR    = DIR_DA_TRAJECTORIES / f"{run}.tpr"

        SB_GRAPH    =  DIR_DA_PYINTERAPH / f"{run}-sb_graph.dat"
        HB_GRAPH    =  DIR_DA_PYINTERAPH / f"{run}-hb_graph.dat"
        HC_GRAPH    =  DIR_DA_PYINTERAPH / f"{run}-HC_graph.dat"

        SB_GRAPH_ALL    =  DIR_DA_PYINTERAPH / f"{run}-sb_graph_all.dat"
        HB_GRAPH_ALL    =  DIR_DA_PYINTERAPH / f"{run}-hb_graph_all.dat"
        HC_GRAPH_ALL    =  DIR_DA_PYINTERAPH / f"{run}-HC_graph_all.dat"

        SB_GRAPH_FILTERED    =  DIR_DA_PYINTERAPH / f"{run}-sb_graph_filtered.dat"
        HB_GRAPH_FILTERED    =  DIR_DA_PYINTERAPH / f"{run}-hb_graph_filtered.dat"
        HC_GRAPH_FILTERED    =  DIR_DA_PYINTERAPH / f"{run}-HC_graph_filtered.dat"

        KBP_GRAPH    =  DIR_DA_PYINTERAPH / f"{run}-kbp-graph.dat"

        CH_GROUPS   =  DIR_DA_PYINTERAPH / f"{run}-ch_groups.ini"
        H_BONDS     =  DIR_DA_PYINTERAPH / f"{run}-H_bonds.ini"

        CLUSTER_PLOT_SB = DIR_DA_PYINTERAPH / f"{run}-clusters_plot_sb.pdf"
        CLUSTER_PLOT_HB = DIR_DA_PYINTERAPH / f"{run}-clusters_plot_hb.pdf"
        CLUSTER_PLOT_HC = DIR_DA_PYINTERAPH / f"{run}-clusters_plot_hc.pdf"

        CLUSTER_SIZE_SB = DIR_DA_PYINTERAPH / f"{run}-clusters_size_sb.dat"
        CLUSTER_SIZE_HB = DIR_DA_PYINTERAPH / f"{run}-clusters_size_hb.dat"
        CLUSTER_SIZE_HC = DIR_DA_PYINTERAPH / f"{run}-clusters_size_hc.dat"
        
        MACRO_IIN_UNWEIGHTED = DIR_DA_PYINTERAPH / f"{run}-macro_IIN_unweighted.dat"

        # salt bridges
        pyinteraph -s PATH_TPR -t PATH_XTC -r PATH_GRO --sb-co 5 -b --sb-graph SB_GRAPH --ff-masses charmm27 -v --sb-cg-file CH_GROUPS

        # H bonds
        pyinteraph -s PATH_TPR -t PATH_XTC -r PATH_GRO -y --hb-graph HB_GRAPH --ff-masses charmm27 -v --hb-ad-file H_BONDS

        # Hydrophobic interaction
        pyinteraph -s PATH_TPR -t PATH_XTC -r PATH_GRO -f --hc-co 5 -f --hc-graph HC_GRAPH --ff-masses charmm27 -v --hc-residues ALA,VAL,LEU,ILE,PHE,PRO,MET,TRP

        # cutoff for occurrence percentage 
        filter_graph -d SB_GRAPH_ALL -c CLUSTER_SIZE_SB -p CLUSTER_PLOT_SB
        filter_graph -d HB_GRAPH_ALL -c CLUSTER_SIZE_HB -p CLUSTER_PLOT_HB
        filter_graph -d HC_GRAPH_ALL -c CLUSTER_SIZE_HC -p CLUSTER_PLOT_HC

        # filtering the graphs using the identified cutoff of the occurrence percentage
        filter_graph -d SB_GRAPH_ALL -o SB_GRAPH_FILTERED -t 20.0
        filter_graph -d HB_GRAPH_ALL -o HB_GRAPH_FILTERED -t 20.0
        filter_graph -d HC_GRAPH_ALL -o HC_GRAPH_FILTERED -t 20.0

        # generating a macroIIN
        filter_graph -d SB_GRAPH_FILTERED -d HC_GRAPH_FILTERED -d SB_GRAPH_FILTERED -o MACRO_IIN_UNWEIGHTED

        # generating a macroIIN weighted
            # generate the energy interaction map
        pyinteraph -s PATH_TPR -t PATH_XTC -r PATH_GRO -p --ff-masses charmm27 -v --kbp-graph KBP_GRAPH
            # generating the macroIIN starting from the KBP_GRAPH
        filter_graph -d hb-graph_filtered.dat -d hc-graph_filtered.dat -d sb-graph_filtered.dat -o macro_IIN_weighted.dat -w kbp-graph.dat