from parameters import *
import matplotlib.pyplot as plt

def vis_rgyr(rgyr_dir):
    rgyr = np.load(rgyr_dir)
    plt.plot(rgyr)
    plt.xlabel("frames")
    plt.ylabel("Radius of Gyration")
    plt.show()


if __name__ == "__main__":
    for run in RUNS:
        PATH_RGYR = DIR_DA_GENERAL / f"{run}-rgyr.npy"
        plt.title(run)
        vis_rgyr(PATH_RGYR)
