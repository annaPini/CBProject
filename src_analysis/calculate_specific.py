from parameters import *
import MDAnalysis as mda
from MDAnalysis.lib import distances

# //////////////////////////////////////////////////////////////////////////////
def calc_pca(traj, pca_space_dir, pca_cumvar_dir, path_pca_pcomps):
    print(f">>> Preparing '{pca_space_dir.name}'...")

    print(">>> Running PCA...")
    rec_pca = PCA(traj, select = "name CA")
    rec_pca.run()
    print(">>> Done.")

    space = rec_pca.transform(traj.select_atoms("name CA"), 3)
    cumvar = rec_pca.cumulated_variance
    pcomps = rec_pca.p_components

    print(f"...>>> Saving '{pca_space_dir.name}'...")
    np.save(pca_space_dir, space, allow_pickle = False)
    np.save(pca_cumvar_dir, cumvar, allow_pickle = False)
    np.save(path_pca_pcomps, pcomps, allow_pickle = False)


def extract_specific_coords(traj, path_coords):
    pass
    # print(f">>> Preparing '{coords_dir.name}'...")
    # coords = np.array([traj.select_atoms("protein and name CA").positions for _ in traj.trajectory[:]])
    # np.save(coords_dir, coords, allow_pickle = False)
    # return coords

# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    PATH_GRO = DIR_DA_TRAJECTORIES / f"{CURRENT_RUN}.gro"
    PATH_XTC = DIR_DA_TRAJECTORIES / f"{CURRENT_RUN}.xtc"

    PATH_PCA_SPACE  = DIR_DA_SPECIFIC / f"{CURRENT_RUN}-pca_space.npy"
    PATH_PCA_CUMVAR = DIR_DA_SPECIFIC / f"{CURRENT_RUN}-pca_cumvar.npy"
    PATH_PCA_PCOMPS = DIR_DA_SPECIFIC / f"{CURRENT_RUN}-pca_pcomponents.npy"

    traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))

    if not PATH_PCA_CUMVAR.exists(): calc_pca(traj, PATH_PCA_SPACE, PATH_PCA_CUMVAR, PATH_PCA_PCOMPS)

# //////////////////////////////////////////////////////////////////////////////
