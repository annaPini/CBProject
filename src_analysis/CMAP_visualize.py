from parameters import *
from calculate_general import calc_cmap
import matplotlib.pyplot as plt

# //////////////////////////////////////////////////////////////////////////////
def vis_cmap(cmap, title = '', color_method = "Reds"):
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom = 0.1, top = 0.9, left = 0.15, right = 0.95)

    ax.set_title(title, fontdict = dict(fontsize = 20))
    ax.set_xlabel("Atom A", fontdict = dict(fontsize = 16))
    ax.set_ylabel("Atom B", fontdict = dict(fontsize = 16))

    colorbar = fig.colorbar( plt.imshow(-cmap, cmap = color_method) )

    ax.tick_params(labelsize = 12)
    colorbar.ax.tick_params(labelsize = 12)

# ------------------------------------------------------------------------------
def diff_cmap(cmap1, cmap2, title = ''):
    vis_cmap(abs(cmap2 - cmap1), title, color_method = "Blues")



# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":

    run0 = "mt1_rep0"
    run1 = "mt2_rep0"
    frame0 = 0
    frame1 = 0

    PATH_COORDS0 = DIR_DA_TRAJECTORIES / f"{run0}-coords.npy"
    PATH_COORDS1 = DIR_DA_TRAJECTORIES / f"{run1}-coords.npy"

    coords0 = np.load(PATH_COORDS0)
    coords1 = np.load(PATH_COORDS1)

    cmap1 = calc_cmap(coords0[frame0])
    cmap2 = calc_cmap(coords1[frame1])

    vis_cmap(cmap1, f"{run0}-frame_{frame0}")
    vis_cmap(cmap2, f"{run1}-frame_{frame1}")
    # diff_cmap(cmap1, cmap2, f"{run0} frame {frame0} vs {run1} frame {frame1}")

    plt.show()



# //////////////////////////////////////////////////////////////////////////////
