from parameters import *
from calculate_general import calc_cmap
import matplotlib.pyplot as plt
import MDAnalysis as mda


def vis_cmap(cmap,run):
    plt.figure(figsize=((10,10)))
    plt.title(run)
    img = plt.imshow(cmap)
    plt.colorbar(img)

def vis_2cmap(cmap1, cmap2, run):
    # cmap1 = np.ones((10,10)) + np.arange(10)
    # cmap2 = np.ones((10,10)) + np.arange(10)*10

    fig, axs = plt.subplots(1,2)
    img1 = axs[0].imshow(cmap1)
    img2 = axs[1].imshow(cmap2)
    fig.colorbar(img1)


def diff_cmap(cmap1, cmap2, run):
    diff = abs(cmap2 - cmap1)

    plt.figure(figsize=((10,10)))
    plt.title(run)
    img = plt.imshow(diff)
    plt.colorbar(img)


if __name__ == "__main__":
    for run in RUNS[:1]:
        PATH_COORDS = DIR_DA_TRAJECTORIES / f"{run}-coords.npy"
        PATH_CMAP1 = DIR_DA_GENERAL / f"{run}-cmap1.npy"
        PATH_CMAP2 = DIR_DA_GENERAL / f"{run}-cmap2.npy"

        f1 = 100
        f2 = 5000


        coords = np.load(PATH_COORDS)
        cmap1 = calc_cmap(coords, PATH_CMAP1, frame = f1)
        cmap2 = calc_cmap(coords, PATH_CMAP2, frame = f2)


        # PATH_GRO     = DIR_DA_TRAJECTORIES / f"{run}.gro"
        # PATH_XTC     = DIR_DA_TRAJECTORIES / f"{run}.xtc"
        # traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))
        #
        # traj.trajectory[f1]
        # coords1 = traj.select_atoms("name CA").positions
        #
        # traj.trajectory[f2]
        # coords2 = traj.select_atoms("name CA").positions

        # calc_cmap(coords1, PATH_CMAP1)#, frame = 100)
        # calc_cmap(coords2, PATH_CMAP2)#, frame = 5000)

        vis_cmap(cmap1, run)
        vis_cmap(cmap2, run)

        # vis_2cmap(cmap1, cmap2, "DIFF")
        diff_cmap(cmap1, cmap2, "DIFF")

        plt.show()

        # compare 2 cmaps
        #compare_cmap(GRO, XTC, 3, 5, CMAP_OUT1, CMAP_OUT2)
        #vis_cmap(CMAP_OUT1, run)
        #vis_cmap(CMAP_OUT2, run)
