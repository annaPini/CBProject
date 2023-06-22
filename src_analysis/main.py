from _params import *
from calculations import *

# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run in RUNS:
        print(f"***** CURRENT RUN: {run} *****")
        PATH_GRO       = DIR_DA_TRAJECTORIES / f"{run}.gro"
        PATH_XTC       = DIR_DA_TRAJECTORIES / f"{run}.xtc"
        PATH_COORDS    = DIR_DA_TRAJECTORIES / f"{run}.coords.npy"
        PATH_COORDS_AS = DIR_DA_TRAJECTORIES / f"{run}-AS.coords.npy"

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

        ########################################################################
        traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))

        coords = np.load(PATH_COORDS) if PATH_COORDS.exists() else extract_ca_coords(traj, PATH_COORDS)
        if not PATH_COORDS_AS.exists(): extract_AS_coords(traj, PATH_COORDS_AS, single_subunit = "t1" in run)

        if not PATH_RAMA.exists(): calc_rama(traj, PATH_RAMA)
        if not PATH_RGYR.exists(): calc_rgyr(traj, PATH_RGYR)
        if not PATH_RMSD.exists(): calc_rmsd(coords, PATH_RMSD)
        if not PATH_RMSF.exists(): calc_rmsf(traj, PATH_RMSF)

        if not PATH_PCA_SPACE.exists(): calc_pca(traj, PATH_PCA_SPACE, PATH_PCA_CUMVAR, PATH_PCA_PCOMPS)
        if not PATH_CLUSTERING.exists(): calc_link(PATH_RMSD, PATH_CLUSTERING, CLUSTERING_LINK_METHOD)

        rmsd_1d = np.load(PATH_RMSD)[0]
        if not PATH_BSE_NAIVE.exists(): calc_bse_naive(rmsd_1d, PATH_BSE_NAIVE)
        if not PATH_BSE_ALTERN.exists(): calc_bse_alternate(rmsd_1d, PATH_BSE_ALTERN)


# //////////////////////////////////////////////////////////////////////////////
