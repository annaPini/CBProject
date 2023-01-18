from parameters import *
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd

# //////////////////////////////////////////////////////////////////////////////
def gen_coords_file(gro_dir, xtc_dir, coords_dir):
    print("...>>> Reading trajectory...")
    traj = mda.Universe(str(gro_dir), str(xtc_dir))

    print("...>>> Extracting coordinates...")
    coords = np.array([traj.select_atoms("name CA").positions for _ in traj.trajectory[:]])

    print(f"...>>> Saving coords.npy file...")
    np.save(coords_dir, coords, allow_pickle = False)

def gen_rmsd_file(coords, rmsd_dir):
    n_frames = coords.shape[0]

    print(f"...>>> Calculating (non-redundant) all-to-all RMSD for {n_frames} frames...")
    result = [[rmsd(coords[i], coords[j]) for j in range(i)] for i in range(1,n_frames)]

    print(f"...>>> Fitting into a {n_frames}x{n_frames} matrix...")
    rmsd_map = np.zeros((n_frames, n_frames))
    for i,res in enumerate(result, 1):
        rmsd_map[i,:i] = rmsd_map[:i,i] = res

    print(f"...>>> Saving rmsd.npy file...")
    np.save(rmsd_dir, rmsd_map, allow_pickle = False)


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run in RUNS:
        #######################################
        POSTDIR = get_POSTDIR(run)
        GRO = POSTDIR / f"{SYSTEM_NAME}.gro"
        XTC = POSTDIR / f"{TRAJECTORY_NAME}.xtc"
        COORDS = POSTDIR / f"{TRAJECTORY_NAME}-coords.npy"
        RMSD = POSTDIR / f"{TRAJECTORY_NAME}-rmsd.npy"

        #######################################
        if RMSD.exists():
            print(f">>> rmsd.npy already exists for run '{run}', skipping...")
        else:
            if not COORDS.exists():
                print(f">>> Couldn't find coords.npy file for run '{run}', creating...")
                gen_coords_file(GRO, XTC, COORDS)

            print(f">>> Couldn't find rmsd.npy file for run '{run}', creating...")
            gen_rmsd_file(np.load(COORDS), RMSD)

# //////////////////////////////////////////////////////////////////////////////
