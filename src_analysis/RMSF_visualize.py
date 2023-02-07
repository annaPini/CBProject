from parameters import *
import matplotlib.pyplot as plt

def vis_rmsf(rmsf_dir):
    rmsf = np.load(RMSF_OUT)
    plt.plot(rmsf)
    print(f"Max rmsf for run {run} is {np.argmax(rmsf)}")

if __name__ == "__main__":
    for run in RUNS[1::2]:
        PATH_RMSF = DIR_DA_GENERAL / f"{run}-rmsf.npy"
        vis_rmsf(PATH_RMSF)

    plt.show()
