from parameters import *
import MDAnalysis as mda
from MDAnalysis.lib import distances
from MDAnalysis.analysis.rms import rmsd, RMSF
from MDAnalysis.analysis.dihedrals import Ramachandran

# //////////////////////////////////////////////////////////////////////////////
def extract_ca_coords(traj, path_coords):
    print(f">>> Preparing '{path_coords.name}'...")
    coords = np.array([traj.select_atoms("protein and name CA").positions for _ in traj.trajectory[:]])
    np.save(path_coords, coords, allow_pickle = False)
    return coords

# ------------------------------------------------------------------------------
def calc_rmsd(coords, path_rmsd):
    print(f">>> Preparing '{path_rmsd.name}'...")
    n_frames = coords.shape[0]

    print(f"...>>> Calculating (non-redundant) all-to-all RMSD for {n_frames} frames...")
    result = [[rmsd(coords[i], coords[j]) for j in range(i)] for i in range(1, n_frames)]

    print(f"...>>> Fitting into a {n_frames}x{n_frames} matrix...")
    rmsd_map = np.zeros((n_frames, n_frames))
    for i,res in enumerate(result, 1):
        rmsd_map[i,:i] = rmsd_map[:i,i] = res

    print(f"...>>> Saving '{path_rmsd.name}'...")
    np.save(path_rmsd, rmsd_map, allow_pickle = False)

# ------------------------------------------------------------------------------
def calc_rmsf(traj, path_rmsf):
    print(f">>> Preparing '{path_rmsf.name}'...")

    ca_atoms = traj.select_atoms("protein and name CA")
    print("...>>> Calculating RMSF...")
    rmsf_CA = RMSF(ca_atoms).run()
    np.save(path_rmsf, rmsf_CA.rmsf, allow_pickle = False)
    print("...>>> Done.")

# ------------------------------------------------------------------------------
def calc_rgyr(traj, path_rgyr):
    print(f">>> Preparing '{path_rgyr.name}'...")

    rgyr = [traj.select_atoms("protein and name CA").radius_of_gyration() for _ in traj.trajectory[:]]
    np.save(path_rgyr, rgyr, allow_pickle = False)
    print("...>>> Done.")

# ------------------------------------------------------------------------------
def calc_rama(traj, path_rama):
    print(f">>> Preparing '{path_rama.name}'...")

    prot = traj.select_atoms("protein")
    rama = Ramachandran(prot).run()
    np.save(path_rama, rama.angles, allow_pickle = False)
    print("...>>> Done.")

# ------------------------------------------------------------------------------
def calc_cmap(coords_frame):
    return distances.distance_array(coords_frame, coords_frame)

def calc_cmap_AS(coords0, coords1):
    return distances.distance_array(coords0, coords1)

# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run in RUNS:
        print(f"***** CURRENT RUN: {run} *****")
        PATH_GRO    = DIR_DA_TRAJECTORIES / f"{run}.gro"
        PATH_XTC    = DIR_DA_TRAJECTORIES / f"{run}.xtc"
        PATH_COORDS = DIR_DA_TRAJECTORIES / f"{run}-coords.npy"

        PATH_RMSD   = DIR_DA_GENERAL      / f"{run}-rmsd.npy"
        PATH_RMSF   = DIR_DA_GENERAL      / f"{run}-rmsf.npy"
        PATH_RGYR   = DIR_DA_GENERAL      / f"{run}-rgyr.npy"
        PATH_RAMA   = DIR_DA_GENERAL      / f"{run}-rama.npy"

        #########################################
        traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))
        coords = np.load(PATH_COORDS) if PATH_COORDS.exists() else extract_ca_coords(traj, PATH_COORDS)

        if not PATH_RMSD.exists(): calc_rmsd(coords, PATH_RMSD)
        if not PATH_RMSF.exists(): calc_rmsf(traj, PATH_RMSF)
        if not PATH_RGYR.exists(): calc_rgyr(traj, PATH_RGYR)
        if not PATH_RAMA.exists(): calc_rama(traj, PATH_RAMA)

# //////////////////////////////////////////////////////////////////////////////
