from parameters import *
import matplotlib.pyplot as plt 
import MDAnalysis as mda
from MDAnalysis.lib import distances

def calc_cmap(gro_dir, traj_dir, out_dir):
    traj = mda.Universe(str(gro_dir), str(traj_dir))
    ca_atoms = traj.select_atoms("protein and name CA")
    print("...>>> Calculating contact map...")
    d_CaCa = distances.distance_array(ca_atoms.positions, ca_atoms.positions)
    np.save(out_dir, d_CaCa, allow_pickle=False)
    print("done")

def vis_cmap(cmap_dir,run):
    cmap = np.load(cmap_dir)
    plt.figure(figsize=((10,10)))
    plt.title(run)
    img = plt.imshow(cmap)
    plt.colorbar(img)
    plt.show()

def compare_cmap(gro_dir, traj_dir, f1, f2, out_dir1, out_dir2):

    traj1 = mda.Universe(str(gro_dir), str(traj_dir))[f1]
    ca_atoms1 = traj1.select_atoms("protein and name CA")
    traj2 = mda.Universe(str(gro_dir), str(traj_dir))[f2]
    ca_atoms2 = traj2.select_atoms("protein and name CA")

    print("...>>> Calculating contact map for frame {f1}...")
    d_CaCa1 = distances.distance_array(ca_atoms1.positions, ca_atoms1.positions)
    np.save(out_dir1, d_CaCa1, allow_pickle=False)

    print("...>>> Calculating contact map frame {f2}...")
    d_CaCa2 = distances.distance_array(ca_atoms2.positions, ca_atoms2.positions)
    np.save(out_dir2, d_CaCa2, allow_pickle=False)

    print("done")


if __name__ == "__main__":

    for run in RUNS:
        POSTDIR = get_POSTDIR(run)
        CMAP_OUT = POSTDIR / f"{TRAJECTORY_NAME}-cmap.npy"
        #CMAP_OUT1 = POSTDIR / f"{TRAJECTORY_NAME}-cmap1.npy"
        #CMAP_OUT2 = POSTDIR / f"{TRAJECTORY_NAME}-cmap2.npy"
        GRO = POSTDIR / f"{SYSTEM_NAME}.gro"
        XTC = POSTDIR / f"{TRAJECTORY_NAME}.xtc"

        if not CMAP_OUT.exists():
            calc_cmap(GRO,XTC,CMAP_OUT)

        #plt.title(run)
        vis_cmap(CMAP_OUT,run)
        # compare 2 cmaps
        #compare_cmap(GRO, XTC, 3, 5, CMAP_OUT1, CMAP_OUT2)
        #vis_cmap(CMAP_OUT1, run)
        #vis_cmap(CMAP_OUT2, run)


