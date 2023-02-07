from parameters import *
import MDAnalysis as mda
from MDAnalysis.lib import distances
from MDAnalysis.analysis.rms import rmsd, RMSF

# //////////////////////////////////////////////////////////////////////////////
def extract_ca_coords(traj, coords_dir):
    print(f">>> Preparing '{coords_dir.name}'...")
    coords = np.array([traj.select_atoms("protein and name CA").positions for _ in traj.trajectory[:]])
    np.save(coords_dir, coords, allow_pickle = False)
    return coords

# ------------------------------------------------------------------------------
def calc_rmsd(coords, rmsd_dir):
    print(f">>> Preparing '{rmsd_dir.name}'...")
    n_frames = coords.shape[0]

    print(f"...>>> Calculating (non-redundant) all-to-all RMSD for {n_frames} frames...")
    result = [[rmsd(coords[i], coords[j]) for j in range(i)] for i in range(1, n_frames)]

    print(f"...>>> Fitting into a {n_frames}x{n_frames} matrix...")
    rmsd_map = np.zeros((n_frames, n_frames))
    for i,res in enumerate(result, 1):
        rmsd_map[i,:i] = rmsd_map[:i,i] = res

    print(f"...>>> Saving rmsd.npy file...")
    np.save(rmsd_dir, rmsd_map, allow_pickle = False)

# ------------------------------------------------------------------------------
def calc_rmsf(traj, rmsf_dir):
    print(f">>> Preparing '{rmsf_dir.name}'...")

    ca_atoms = traj.select_atoms("protein and name CA")
    print("calculating RMSF...")
    rmsf_CA = RMSF(ca_atoms).run()
    np.save(rmsf_dir, rmsf_CA.rmsf, allow_pickle = False)
    print("done")

# ------------------------------------------------------------------------------
def calc_rgyr(traj, rgyr_dir):
    print(f">>> Preparing '{rgyr_dir.name}'...")

    rgyr = [traj.select_atoms("protein and name CA").radius_of_gyration() for _ in traj.trajectory[:]]
    np.save(rgyr_dir, rgyr, allow_pickle = False)
    print("done")

# ------------------------------------------------------------------------------
def calc_cmap(coords, cmap_dir):
    print(f">>> Preparing '{cmap_dir.name}'...")

    print("...>>> Calculating contact map...")
    d_CaCa = distances.distance_array(coords, coords)
    np.save(cmap_dir, d_CaCa, allow_pickle = False)
    print("done")

# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run in RUNS[:1]:
        PATH_GRO    = DIR_DA_TRAJECTORIES / f"{run}.gro"
        PATH_XTC    = DIR_DA_TRAJECTORIES / f"{run}.xtc"
        PATH_COORDS = DIR_DA_TRAJECTORIES / f"{run}-coords.npy"

        PATH_RMSD   = DIR_DA_GENERAL      / f"{run}-rmsd.npy"
        PATH_RMSF   = DIR_DA_GENERAL      / f"{run}-rmsf.npy"
        PATH_RGYR   = DIR_DA_GENERAL      / f"{run}-rgyr.npy"
        PATH_CMAP   = DIR_DA_GENERAL      / f"{run}-cmap.npy"

        #########################################
        traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))
        coords = np.load(PATH_COORDS) if PATH_COORDS.exists() else extract_ca_coords(traj, PATH_COORDS)

        if not PATH_RMSD.exists(): calc_rmsd(coords, PATH_RMSD)
        if not PATH_RMSF.exists(): calc_rmsf(traj, PATH_RMSF)
        if not PATH_RGYR.exists(): calc_rgyr(traj, PATH_RGYR)
        if not PATH_CMAP.exists(): calc_cmap(coords[0], PATH_CMAP)

# //////////////////////////////////////////////////////////////////////////////
