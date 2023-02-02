from parameters import *
import matplotlib.pyplot as plt 
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF

def calc_rmsf(gro_dir, traj_dir, out_dir):
    traj = mda.Universe(str(gro_dir), str(traj_dir))
    ca_atoms = traj.select_atoms("protein and name CA")
    print("calculating RMSF...")
    rmsf_CA = RMSF(ca_atoms).run()
    np.save(out_dir, rmsf_CA.rmsf, allow_pickle=False)
    print("done")

def vis_rmsf(rmsf_dir):
    rmsf = np.load(rmsf_dir)

    plt.plot(rmsf)
    # plt.show()
    # plt.figure()

if __name__ == "__main__":
    for run in RUNS[1::2]:
        POSTDIR = get_POSTDIR(run)
        RMSF_OUT = POSTDIR / f"{TRAJECTORY_NAME}-rmsf.npy"
        GRO = POSTDIR / f"{SYSTEM_NAME}.gro"
        XTC = POSTDIR / f"{TRAJECTORY_NAME}.xtc"
        
        if not RMSF_OUT.exists():
            calc_rmsf(GRO, XTC, RMSF_OUT)

        # plt.title(run)
        vis_rmsf(RMSF_OUT)

        rmsf=np.load(RMSF_OUT)
        print(f"Max rmsf for run {run} is {np.argmax(rmsf)}")

    plt.show()