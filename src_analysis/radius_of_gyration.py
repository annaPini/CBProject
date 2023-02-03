from parameters import *
import numpy as np
import matplotlib.pyplot as plt

def calc_rgyr(gro_dir, traj_dir, out_dir):
    traj = mda.Universe(str(gro_dir),str(traj_dir))
    ca_atoms = traj.select_atoms("protein and name CA")
    print("...>>> Calculating Radius of Gyration...")
    rgyr_CA = ca_atoms.radius_of_gyration()
    np.save(out_dir, rgyr_CA, allow_pickle=False)
    print("done")

def vis_rgyr(rgyr_dir):
    rgyr = np.load(rgyr_dir)
    plt.plot(rgyr)
    plt.show()

def show_rgyr(gro_dir,traj_dir, run):
    print(f"...>>> Radius of gyration for {run}...")
    traj = mda.Universe(str(gro_dir),str(traj_dir))
    ca_atoms = traj.select_atoms("protein and name CA")
    for ts in traj.trajectory[:50]:
        time = traj.trajectory.time
        rgyr = ca_atoms.radius_of_gyration()
        print(f"Frame: {ts.frame:3d}, Time: {time:4.0f} ps, Rgyr: {rgyr:.4f}")
    print("done")


if __name__ == "__main__":

    for run in RUNS:
        POSTDIR = get_POSTDIR(run)
        RGYR_OUT = POSTDIR / f"{TRAJECTORY_NAME}-rgyr.npy"
        GRO = POSTDIR / f"{SYSTEM_NAME}.gro"
        XTC = POSTDIR / f"{TRAJECTORY_NAME}.xtc"

        if not RGYR_OUT.exists():
           calc_rgyr(GRO,XTC,RGYR_OUT)

    plt.title(run)
    vis_rgyr(RGYR_OUT)
    
    # show_rgyr(GRO,XTC,run)