import argparse
from _params import *
from calculations import *
from visualizations import *
from documentation import docs_main
from argparse import RawDescriptionHelpFormatter

# //////////////////////////////////////////////////////////////////////////////
def run_calculations():
    for run in RUNS:
        print(f"***** CURRENT RUN: {run} *****")
        PATH_GRO       = DIR_DA_TRAJECTORIES / f"{run}.gro"
        PATH_XTC       = DIR_DA_TRAJECTORIES / f"{run}.xtc"
        PATH_COORDS    = DIR_DA_COORDS / f"{run}.npy"
        PATH_COORDS_AS = DIR_DA_COORDS / f"{run}-AS0.npy"

        PATH_BSE_NAIVE  = DIR_DA_BSE  / f"{run}.naive.npy"
        PATH_BSE_ALTERN = DIR_DA_BSE  / f"{run}.altern.npy"
        PATH_CLUSTERING = DIR_DA_CLUSTERING / f"{run}.link.npy"
        PATH_PCA_SPACE  = DIR_DA_PCA  / f"{run}.space.npy"
        PATH_PCA_CUMVAR = DIR_DA_PCA  / f"{run}.cumvar.npy"
        PATH_PCA_PCOMPS = DIR_DA_PCA  / f"{run}.pcomponents.npy"
        PATH_RAMA       = DIR_DA_RAMA / f"{run}.npy"
        PATH_RGYR       = DIR_DA_RGYR / f"{run}.npy"
        PATH_RMSD       = DIR_DA_RMSD / f"{run}.npy"
        PATH_RMSF       = DIR_DA_RMSF / f"{run}.npy"


        traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))
        coords = np.load(PATH_COORDS) if PATH_COORDS.exists() else extract_ca_coords(traj, PATH_COORDS)
        if not PATH_COORDS_AS.exists(): extract_AS_coords(traj, PATH_COORDS_AS, single_subunit = "t1" in run)

        if not PATH_RAMA.exists(): calc_rama(traj, PATH_RAMA)
        if not PATH_RGYR.exists(): calc_rgyr(traj, PATH_RGYR)
        if not PATH_RMSD.exists(): calc_rmsd(coords, PATH_RMSD)
        if not PATH_RMSF.exists(): calc_rmsf(traj, PATH_RMSF)

        if not PATH_PCA_SPACE.exists(): calc_pca(traj, PATH_PCA_SPACE, PATH_PCA_CUMVAR, PATH_PCA_PCOMPS)
        if not PATH_CLUSTERING.exists(): calc_link(PATH_RMSD, PATH_CLUSTERING, CLUSTERING_LINK_METHOD)

        if not PATH_BSE_NAIVE.exists(): calc_bse_naive(PATH_RMSD, PATH_BSE_NAIVE)
        if not PATH_BSE_ALTERN.exists(): calc_bse_alternate(PATH_RMSD, PATH_BSE_ALTERN)


# ------------------------------------------------------------------------------
def vis_bse():
    plotter = Plotter_BSE()
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        plotter.vis_BSE(
            bse_naive_0 = np.load(DIR_DA_BSE / f"{run0}.naive.npy"),
            bse_altern_0 = np.load(DIR_DA_BSE / f"{run0}.altern.npy"),
            bse_naive_1 = np.load(DIR_DA_BSE / f"{run1}.naive.npy"),
            bse_altern_1 = np.load(DIR_DA_BSE / f"{run1}.altern.npy"),
            title = run_preffix,
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_cluster(plotter_objs):
    for run_preffix in RUN_PREFFIXES:
        plotter_objs.append(Plotter_RMSD1D_Clustering())
        run0 = run_preffix + "_rep0"

        plotter_objs[-1].vis_rmsd1d_clustering(
            rmsd_mat0 = np.load(DIR_DA_RMSD / f"{run0}.npy"),
            Z = np.load(DIR_DA_CLUSTERING / f"{run0}.link.npy"),
            color_map = COLOR_MAP_RMSD,
            title = run0
        )


# ------------------------------------------------------------------------------
def vis_cmap0(): print("placeholder cmap0")


# ------------------------------------------------------------------------------
def vis_cmap1(): print("placeholder cmap1")


# ------------------------------------------------------------------------------
def vis_cmap2(): print("placeholder cmap2")


# ------------------------------------------------------------------------------
def vis_pca(plotter_objs):
    for run_preffix in RUN_PREFFIXES:
        plotter_objs.append(Plotter_PCA())
        run0 = run_preffix + "_rep0"

        cumvar = np.load(DIR_DA_PCA / f"{run0}.cumvar.npy")
        n_pcs = np.where(cumvar > 0.95)[0][0]
        print("...>>> NPCS:", n_pcs)

        ### Alternatively, can use the 'vis_1pca_plotly' method
        ### With it, interaction with 3D scatter plot is faster
        ### However, it opens several windows, so it's advised to use it for only 1 run

        plotter_objs[-1].vis_1pca(
            cumvar = cumvar,
            space = np.load(DIR_DA_PCA / f"{run0}.space.npy"),
            pcomps = np.load(DIR_DA_PCA / f"{run0}.pcomponents.npy"),
            title = run_preffix
        )


# ------------------------------------------------------------------------------
def vis_pyinteraph(): print("placeholder pyinteraph")


# ------------------------------------------------------------------------------
def vis_rama(plotter_objs):
    for run_preffix in RUN_PREFFIXES:
        plotter_objs.append(Plotter_RAMA())
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        plotter_objs[-1].vis_2rama(
            rama0 = np.load(DIR_DA_RAMA / f"{run0}.npy"),
            rama1 = np.load(DIR_DA_RAMA / f"{run1}.npy"),
            title = run_preffix,
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_rgyr():
    plotter = Plotter_RGYR()
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"
        plotter.vis_2rgyr(
            rgyr0 = np.load(DIR_DA_RGYR / f"{run0}.npy"),
            rgyr1 = np.load(DIR_DA_RGYR / f"{run1}.npy"),
            title = run_preffix,
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_rmsd0(plotter_objs):
    """RMSD 2D (memory intensive)"""

    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        rmsd_mat0 = np.load(DIR_DA_RMSD / f"{run0}.npy")
        rmsd_mat1 = np.load(DIR_DA_RMSD / f"{run1}.npy")

        plotter_rmsd2d = Plotter_RMSD2D()
        plotter_rmsd2d.vis_rmsd2d(rmsd_mat0, title = run0)
        plotter_rmsd2d.vis_rmsd2d(rmsd_mat1, title = run1)


# ------------------------------------------------------------------------------
def vis_rmsd1(plotter_objs):
    """RMSD 1D compare (intra comparisons)"""

    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        rmsd_mat0 = np.load(DIR_DA_RMSD / f"{run0}.npy")
        rmsd_mat1 = np.load(DIR_DA_RMSD / f"{run1}.npy")

        plotter_objs.append(Plotter_RMSD1D_Compare())
        plotter_objs[-1].vis_rmsd1d(
            rmsd_mat0, rmsd_mat1, run_preffix, run0, run1
        )


# ------------------------------------------------------------------------------
def vis_rmsd2(plotter_objs):
    """RMSD 1D compare (inter comparisons)"""

    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        rmsd_mat0 = np.load(DIR_DA_RMSD / f"{run0}.npy")
        rmsd_mat1 = np.load(DIR_DA_RMSD / f"{run1}.npy")

        # WIP... (placeholder)
        plotter_objs.append(Plotter_RMSD1D())
        plotter_objs[-1].vis_rmsd1d(rmsd_mat0, title = run0)


# ------------------------------------------------------------------------------
def vis_rmsf0():
    """RMSF 1D compare (intra comparisons)"""

    plotter = Plotter_RMSF()
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        plotter.vis_2rmsf(
            rmsf0 = np.load(DIR_DA_RMSF / f"{run0}.npy"),
            rmsf1 = np.load(DIR_DA_RMSF / f"{run1}.npy"),
            title = run0[:3] if (run0[:3] == run1[:3]) else f"{run0} vs {run1}",
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_rmsf1():
    """RMSF 1D compare (inter comparisons)"""

    comparisons = [
        ("mt2_rep0", "wt2_rep0"),
        ("mt1_rep0", "mt2_rep0"),
        ("wt1_rep0", "wt2_rep0"),
        ("mt2_rep0", "wt1_rep0"),
    ]

    plotter = Plotter_RMSF()
    for run0, run1 in comparisons:
        plotter.vis_2rmsf(
            rmsf0 = np.load(DIR_DA_RMSF / f"{run0}.npy"),
            rmsf1 = np.load(DIR_DA_RMSF / f"{run1}.npy"),
            title = run0[:3] if (run0[:3] == run1[:3]) else f"{run0} vs {run1}",
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_sasa(): print("placeholder sasa")


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = docs_main["main"], formatter_class = RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--calc", action = "store_true", help = docs_main["calc"])
    parser.add_argument("-b", "--bse", action = "store_true", help = docs_main["bse"])
    parser.add_argument("-l", "--cluster", action = "store_true", help = docs_main["cluster"])
    parser.add_argument("-m0", "--cmap0", action = "store_true", help = docs_main["cmap"]) # TODO GENERAL ANALYSIS
    parser.add_argument("-m1", "--cmap1", action = "store_true", help = docs_main["cmap"]) # TODO SPECIFIC ANALYSIS
    parser.add_argument("-m2", "--cmap2", action = "store_true", help = docs_main["cmap"]) # TODO SPECIFIC ANALYSIS ACTIVE SITE
    parser.add_argument("-p", "--pca", action = "store_true", help = docs_main["pca"])
    parser.add_argument("-y", "--pyinteraph", action = "store_true", help = docs_main["pyinteraph"]) # TODO
    parser.add_argument("-a", "--rama", action = "store_true", help = docs_main["rama"])
    parser.add_argument("-g", "--rgyr", action = "store_true", help = docs_main["rgyr"])
    parser.add_argument("-d0", "--rmsd0", action = "store_true", help = docs_main["rmsd"]) # 2D
    parser.add_argument("-d1", "--rmsd1", action = "store_true", help = docs_main["rmsd"]) # 1D COMPARE INTRA
    parser.add_argument("-d2", "--rmsd2", action = "store_true", help = docs_main["rmsd"]) # 1D COMPARE INTER
    parser.add_argument("-f0", "--rmsf0", action = "store_true", help = docs_main["rmsf"]) # COMPARE INTRA
    parser.add_argument("-f1", "--rmsf1", action = "store_true", help = docs_main["rmsf"]) # COMPARE INTER
    parser.add_argument("-s", "--sasa", action = "store_true", help = docs_main["sasa"]) # TODO
    args = parser.parse_args()

    plot_flags = {**vars(args)}
    plot_flags.pop("calc")
    plot_requested = any(plot_flags.values())

    if not (args.calc | plot_requested):
        print(">>> Refer to the documentation (main.py -h) for description of usage. Closing...")
        exit()

    if args.calc: run_calculations()

    if plot_requested:
        if not any(DIR_DA_COORDS.iterdir()):
            print("xxx Visualization(s) requested, but it seems calculations haven't been perform yet (_coords folder empty). Aborting...")
            exit()
    else:
        print("/// No visualizations requested, closing program...")
        exit()

    ### must keep an active pointer of the plotter objects to avoid
    ### garbage collector to dispose of them, which would break interactivity
    ### i.e. functionality of widgets (sliders)
    plotter_objs = []
    if args.bse: vis_bse()
    if args.cluster: vis_cluster(plotter_objs)
    if args.cmap0: vis_cmap0()
    if args.cmap1: vis_cmap1()
    if args.cmap2: vis_cmap2()
    if args.pca: vis_pca(plotter_objs)
    if args.pyinteraph: vis_pyinteraph()
    if args.rama: vis_rama(plotter_objs)
    if args.rgyr: vis_rgyr()
    if args.rmsd0: vis_rmsd0(plotter_objs)
    if args.rmsd1: vis_rmsd1(plotter_objs)
    if args.rmsd2: vis_rmsd2(plotter_objs)
    if args.rmsf0: vis_rmsf0()
    if args.rmsf1: vis_rmsf1()
    if args.sasa: vis_sasa()

    plt.show()


# //////////////////////////////////////////////////////////////////////////////
