import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from parameters import *

def get_full_matrix(path_full_matrix):
    return np.genfromtxt(path_full_matrix)

def get_reduced_matrix(path_reduced_matrix):
    header = ["res0_chain", "res0_num", "res0_name", "res0_inter_group", "res1_chain", "res1_num", "res1_name", "res1_inter_group", "occurrence"]
    df = pd.read_csv(path_reduced_matrix, header = None, usecols = range(len(header)), names = header)

    atoms = sorted(set(df.res0_num) | set(df.res1_num))
    matrix = np.zeros((len(atoms), len(atoms)))
    atom_to_index = {atom:id for id,atom in enumerate(atoms)}

    for _,row in df.iterrows():
        atom_0_id = atom_to_index[row["res0_num"]]
        atom_1_id = atom_to_index[row["res1_num"]]
        value = row["occurrence"]

        matrix[atom_0_id, atom_1_id] = value
        matrix[atom_1_id, atom_0_id] = value
    return matrix, atoms

def plot_heatmap(full_matrix, title = "", cmap = "hot"):
    _, ax = plt.subplots()

    kwargs = {
        # "vmin" : 0, "vmax" : 100
        }

    im = ax.imshow(full_matrix, cmap = cmap, interpolation = "nearest", **kwargs)
    ax.figure.colorbar(im, ax = ax)
    ax.set_title(title)
    return ax

def plot_heatmap_labels(reduced_matrix, atoms, title = ""):
    ax = plot_heatmap(reduced_matrix, title)
    ax.set_xticks(np.arange(len(atoms)), labels = atoms)
    ax.set_yticks(np.arange(len(atoms)), labels = atoms)



####################################################################################
if __name__ == "__main__":
    ####################
    # sb_mat = get_full_matrix("sb-graph_all.dat")
    # hc_mat = get_full_matrix("hc-graph_all.dat")
    # hb_mat = get_full_matrix("hb-graph_all.dat")
    #
    # all_mat = sb_mat + hc_mat + hb_mat
    # all_mat[all_mat > 100] = 100
    #
    # all_mat_3_colors = np.zeros(sb_mat.shape)
    # all_mat_3_colors[sb_mat > 0] = 25
    # all_mat_3_colors[hc_mat > 0] = 50
    # all_mat_3_colors[hb_mat > 0] = 100
    #
    # # plot_heatmap(sb_mat, "Salt Bridges")
    # # plot_heatmap(hc_mat, "Hydrophobic Clusters")
    # # plot_heatmap(hb_mat, "Hydrogen Bonds")
    #
    # plot_heatmap(all_mat, "MERGED")
    # plot_heatmap(all_mat_3_colors, "ALL", cmap = "brg")
    #
    # macro_IIN_unweighted = get_full_matrix("macro_IIN_unweighted.dat")
    # plot_heatmap(macro_IIN_unweighted, "macro_IIN_unweighted")

    ####################
    # mat0 = get_full_matrix("hc-graph_all.dat")
    # plot_heatmap(mat0, "HB all")
    #
    # mat1 = get_full_matrix("hc-graph_filtered.dat")
    # plot_heatmap(mat1, "HB filtered")
    #
    # mat_diff = mat0 - mat1
    # plot_heatmap(mat_diff, "Difference",
    #     # cmap = "Greys"
    #     )

    SB_GRAPH_ALL =  DIR_DA_PYINTERAPH / f"{RUN_DETAILED_ANALYSIS}-salt-bridges_all.csv"

    ####################
    sb_mat, sb_atoms = get_reduced_matrix(SB_GRAPH_ALL)
    # hc_mat, hc_atoms = get_reduced_matrix("hydrophobic-clusters_all.csv")
    # hb_mat, hb_atoms = get_reduced_matrix("hydrogen-bonds_all.csv")
    #
    plot_heatmap_labels(sb_mat, sb_atoms, "Salt Bridges")
    # plot_heatmap_labels(hc_mat, hc_atoms, "Hydrophobic Clusters")
    # plot_heatmap_labels(hb_mat, hb_atoms, "Hydrogen Bonds")

    ####################

    plt.show()
