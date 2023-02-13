from parameters import *

dir_pyinteraph = Path("pyinteraph")
dir_trajectories = Path("trajectories")

CH_GROUPS =  dir_pyinteraph / "charged_groups.ini"
H_BONDS   =  dir_pyinteraph / "hydrogen_bonds.ini"

PATH_GRO = dir_trajectories / f"{CURRENT_RUN}.gro"
PATH_XTC = dir_trajectories / f"{CURRENT_RUN}.xtc"
PATH_TPR = dir_trajectories / f"{CURRENT_RUN}.tpr"

SB_GRAPH =  dir_pyinteraph / f"{CURRENT_RUN}-sb_graph.dat"
HB_GRAPH =  dir_pyinteraph / f"{CURRENT_RUN}-hb_graph.dat"
HC_GRAPH =  dir_pyinteraph / f"{CURRENT_RUN}-HC_graph.dat"

SB_GRAPH_ALL =  dir_pyinteraph / f"{CURRENT_RUN}-sb_graph_all.dat"
HB_GRAPH_ALL =  dir_pyinteraph / f"{CURRENT_RUN}-hb_graph_all.dat"
HC_GRAPH_ALL =  dir_pyinteraph / f"{CURRENT_RUN}-HC_graph_all.dat"

SB_GRAPH_FILTERED =  dir_pyinteraph / f"{CURRENT_RUN}-sb_graph_filtered.dat"
HB_GRAPH_FILTERED =  dir_pyinteraph / f"{CURRENT_RUN}-hb_graph_filtered.dat"
HC_GRAPH_FILTERED =  dir_pyinteraph / f"{CURRENT_RUN}-HC_graph_filtered.dat"

KBP_GRAPH =  dir_pyinteraph / f"{CURRENT_RUN}-kbp-graph.dat"


CLUSTER_PLOT_SB = dir_pyinteraph / f"{CURRENT_RUN}-clusters_plot_sb.pdf"
CLUSTER_PLOT_HB = dir_pyinteraph / f"{CURRENT_RUN}-clusters_plot_hb.pdf"
CLUSTER_PLOT_HC = dir_pyinteraph / f"{CURRENT_RUN}-clusters_plot_hc.pdf"

CLUSTER_SIZE_SB = dir_pyinteraph / f"{CURRENT_RUN}-clusters_size_sb.dat"
CLUSTER_SIZE_HB = dir_pyinteraph / f"{CURRENT_RUN}-clusters_size_hb.dat"
CLUSTER_SIZE_HC = dir_pyinteraph / f"{CURRENT_RUN}-clusters_size_hc.dat"

MACRO_IIN_UNWEIGHTED = dir_pyinteraph / f"{CURRENT_RUN}-macro_IIN_unweighted.dat"

PATH_SH = "pyinteraph.sh"

command = '\n'.join([
    "cd ../data_analysis",

    #################################
    # salt bridges
    # f"pyinteraph -s {PATH_TPR} -t {PATH_XTC} -r {PATH_GRO} --sb-co 5 -b --sb-graph {SB_GRAPH} --ff-masses charmm27 -v --sb-cg-file {CH_GROUPS}"
    #
    # # H bonds
    # f"pyinteraph -s {PATH_TPR} -t {PATH_XTC} -r {PATH_GRO} -y --hb-graph {HB_GRAPH} --ff-masses charmm27 -v --hb-ad-file {H_BONDS}"
    #
    # # Hydrophobic interaction
    # f"pyinteraph -s {PATH_TPR} -t {PATH_XTC} -r {PATH_GRO} -f --hc-co 5 -f --hc-graph {HC_GRAPH} --ff-masses charmm27 -v --hc-residues ALA,VAL,LEU,ILE,PHE,PRO,MET,TRP"

    # # cutoff for occurrence percentage
    f"filter_graph -d {SB_GRAPH_ALL} -c {CLUSTER_SIZE_SB} -p {CLUSTER_PLOT_SB}",
    # filter_graph -d HB_GRAPH_ALL -c CLUSTER_SIZE_HB -p CLUSTER_PLOT_HB
    f"filter_graph -d {HC_GRAPH_ALL} -c {CLUSTER_SIZE_HC} -p {CLUSTER_PLOT_HC}",


    # # filtering the graphs using the identified cutoff of the occurrence percentage
    f"filter_graph -d {SB_GRAPH_ALL} -o {SB_GRAPH_FILTERED} -t 20.0",
    # filter_graph -d HB_GRAPH_ALL -o HB_GRAPH_FILTERED -t 20.0
    f"filter_graph -d {HC_GRAPH_ALL} -o {HC_GRAPH_FILTERED} -t 20.0",
    #
    # # generating a macroIIN
    # filter_graph -d HB_GRAPH_FILTERED -d HC_GRAPH_FILTERED -d SB_GRAPH_FILTERED -o MACRO_IIN_UNWEIGHTED
    #
    # # generating a macroIIN weighted
    #     # generate the energy interaction map
    # pyinteraph -s PATH_TPR -t PATH_XTC -r PATH_GRO -p --ff-masses charmm27 -v --kbp-graph KBP_GRAPH
    #     # generating the macroIIN starting from the KBP_GRAPH
    # filter_graph -d hb-graph_filtered.dat -d hc-graph_filtered.dat -d sb-graph_filtered.dat -o macro_IIN_weighted.dat -w kbp-graph.dat
])


with open(PATH_SH, 'w') as file: file.write(command.replace("\\", "/"))
