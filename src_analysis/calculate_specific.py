from parameters import *
import MDAnalysis as mda
from MDAnalysis.analysis.pca import PCA, cosine_content
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

# //////////////////////////////////////////////////////////////////////////////
def calc_link(path_rmsd, path_link, link_method):
    print(f">>> Preparing '{path_link.name}'...")
    rmsd_mat = np.load(path_rmsd)
    Z = linkage(rmsd_mat, link_method)
    np.save(path_link, Z, allow_pickle = False)
    print("...>>> Done.")

def calc_cluster(Z, t, label_criterion):
    return fcluster(Z, t = t, criterion = label_criterion)

# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
def extract_AS_coords(traj, path_coords, single_subunit = False): # active site coordinates
    print(f">>> Preparing '{path_coords.name}'...")

    selections = [
        "resid 145 and name SG",
        "resid 41 and name NE2",
        "resid 41 and name CE1",
        "resid 164 and name O",
        "resid 41 and name ND1",
        "resid 164 and name ND?",
        "resid 187 and name OD2",
        "resid 40 and name NE",
        "resid 187 and name OD1",
        "resid 40 and name NH2",
    ]

    # create an atom group
    atoms = sum(traj.select_atoms(sel) for sel in selections)

    coords = np.zeros((len(traj.trajectory), len(atoms), 3))
    for i,_ in enumerate(traj.trajectory[:]):
        for j,atom in enumerate(atoms):
            coords[i,j] = atom.position

    if single_subunit:
        np.save(PATH_AS_COORDS, coords, allow_pickle = False)
    else:
        np.save(PATH_AS_COORDS, coords[:, 0::2, :], allow_pickle = False)
        np.save(PATH_AS_COORDS.parent / f"{PATH_AS_COORDS.stem}1.npy", coords[:, 1::2, :], allow_pickle = False)

# ------------------------------------------------------------------------------
def calc_cmap_AS(coords_AS):
    pass

# ------------------------------------------------------------------------------
# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run in RUNS:
        print(f"***** CURRENT RUN: {run} *****")
        PATH_GRO  = DIR_DA_TRAJECTORIES / f"{run}.gro"
        PATH_XTC  = DIR_DA_TRAJECTORIES / f"{run}.xtc"

        PATH_RMSD = DIR_DA_GENERAL / f"{run}-rmsd.npy"
        PATH_LINK = DIR_DA_SPECIFIC / f"{run}-link.npy"

        PATH_PCA_SPACE  = DIR_DA_SPECIFIC / f"{run}-pca_space.npy"
        PATH_PCA_CUMVAR = DIR_DA_SPECIFIC / f"{run}-pca_cumvar.npy"
        PATH_PCA_PCOMPS = DIR_DA_SPECIFIC / f"{run}-pca_pcomponents.npy"

        PATH_AS_COORDS = DIR_DA_SPECIFIC / f"{run}-AS_coords.npy"

        #########################################
        traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))

        if not PATH_LINK.exists(): calc_link(PATH_RMSD, PATH_LINK, link_method = "ward")
        if not PATH_PCA_CUMVAR.exists(): calc_pca(traj, PATH_PCA_SPACE, PATH_PCA_CUMVAR, PATH_PCA_PCOMPS)
        if not PATH_AS_COORDS.exists(): extract_AS_coords(traj, PATH_AS_COORDS, single_subunit = "t1" in run)

# //////////////////////////////////////////////////////////////////////////////
