from parameters import *
import matplotlib.pyplot as plt

def vis_rgyr(rgyr_dir):
    rgyr = np.load(rgyr_dir)
    plt.plot(rgyr)
    plt.xlabel("frames")
    plt.ylabel("Radius of Gyration")
    plt.show()

def vis_2rgyr(rgyr_dir1,rgyr_dir2):
    rgyr1 = np.load(rgyr_dir1)
    rgyr2 = np.load(rgyr_dir2)

    if rgyr1.size != rgyr2.size:
        arrs = sorted([rgyr1, rgyr2], key=lambda arr:arr.size)
        arrs[0] = np.append(arrs[0],arrs[0])

        plt.plot(arrs[0])
        plt.plot(arrs[1])

    else:
       
        plt.plot(rgyr1)
        plt.plot(rgyr2) 

    plt.show()

if __name__ == "__main__":
    for run in RUNS:
        PATH_RGYR = DIR_DA_GENERAL / f"{run}-rgyr.npy"
        plt.title(run)
        vis_rgyr(PATH_RGYR)
# vis_2rgyr(rgyr_dir1,rgyr_dir2)
