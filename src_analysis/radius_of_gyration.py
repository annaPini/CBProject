from parameters import *
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt

def calc_rgyr(gro_dir,traj_dir,out_dir):
    traj = mda.Universe(str(gro_dir),str(traj_dir))
    ca_atoms = traj.select_atoms("protein and name CA")
    print("...>>> Calculating Radius of Gyration...")
    Rgyr = []
    for ts in traj.trajectory:
        Rgyr.append(ca_atoms.radius_of_gyration())
    print("...>>> Saving RGYR...")
    np.save(out_dir,Rgyr,allow_pickle=False)
    print("done")

def vis_rgyr(rgyr_dir):
    rgyr=np.load(rgyr_dir)
    plt.plot(rgyr)
    plt.xlabel("frames")
    plt.ylabel("Radius of Gyration")
    plt.show()


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
    