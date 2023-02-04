from parameters import *
import MDAnalysis as mda
import matplotlib.pyplot as plt

def calc_rgyr(gro_dir,traj_dir,out_dir):
    traj = mda.Universe(str(gro_dir),str(traj_dir))

    print("...>>> Calculating Radius of Gyration...")
    Rgyr = [traj.select_atoms("protein and name CA").radius_of_gyration() for _ in traj.trajectory[:]]

    print("...>>> Saving RGYR...")
    np.save(out_dir,Rgyr,allow_pickle=False)
    print("done")

def vis_rgyr(rgyr_dir):
    rgyr = np.load(rgyr_dir)
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
