from parameters import *
import matplotlib.pyplot as plt

def vis_rmsf(rmsf_dir):
    rmsf = np.load(rmsf_dir)
    plt.plot(rmsf)
    plt.xlabel("atoms")
    plt.ylabel("RMSF")
    plt.show()
    print(f"Max rmsf for run {run} is {np.argmax(rmsf)}")

def vis_2rmsf(rmsf_dir1,rmsf_dir2):
    rmsf1 = np.load(rmsf_dir1)
    rmsf2 = np.load(rmsf_dir2)

    if rmsf1.size != rmsf2.size:
        arrs = sorted([rmsf1, rmsf2], key=lambda arr:arr.size)
        arrs[0] = np.append(arrs[0],arrs[0])
        plt.plot(arrs[0])

    else:

        plt.plot(rmsf1)
        plt.plot(rmsf2)

    plt.show()

if __name__ == "__main__":
    for run in RUNS:
        PATH_RMSF = DIR_DA_GENERAL / f"{run}-rmsf.npy"
        vis_rmsf(PATH_RMSF)

# vis_2rmsf(rmsf_dir1,rmsf_dir2)
