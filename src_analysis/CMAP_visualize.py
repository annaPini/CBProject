from parameters import *
from calculate_general import calc_cmap
import matplotlib.pyplot as plt

# //////////////////////////////////////////////////////////////////////////////
def vis_cmap(cmap, title = '', color_method = "Reds"):
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom = 0.1, top = 0.9, left = 0.15, right = 0.95)

    ax.set_title(title, fontdict = dict(fontsize = 20))
    ax.set_xlabel("CA atom", fontdict = dict(fontsize = 16))
    ax.set_ylabel("CA atom", fontdict = dict(fontsize = 16))

    colorbar = fig.colorbar( plt.imshow(-cmap, cmap = color_method) )

    ax.tick_params(labelsize = 12)
    colorbar.ax.tick_params(labelsize = 12)

# ------------------------------------------------------------------------------
def diff_cmap(cmap1, cmap2, title = ''):
    vis_cmap(abs(cmap2 - cmap1), title, color_method = "Blues")



# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    coords_mt1 = np.load(DIR_DA_TRAJECTORIES / "mt1_rep0-coords.npy")
    coords_mt2 = np.load(DIR_DA_TRAJECTORIES / "mt2_rep0-coords.npy")


    # --------------------------------------------------------------------------
    cmap_mt1_f0 = calc_cmap(coords_mt1[0])
    cmap_mt2_f0 = calc_cmap(coords_mt2[0])

    vis_cmap(cmap_mt1_f0, "mt1_rep0-frame_0")
    vis_cmap(cmap_mt2_f0, "mt2_rep0-frame_0")

    # --------------------------------------------------------------------------
    cmap_mt2_f5472 = calc_cmap(coords_mt2[5472]) # before opening
    cmap_mt2_f5510 = calc_cmap(coords_mt2[5510]) # open
    cmap_mt2_f5692 = calc_cmap(coords_mt2[5692]) # after opening

    # vis_cmap(cmap_mt2_f5472, "mt2_rep0-frame_5472")
    # vis_cmap(cmap_mt2_f5510, "mt2_rep0-frame_5510")
    # vis_cmap(cmap_mt2_f5692, "mt2_rep0-frame_5692")
    diff_cmap(cmap_mt2_f5472, cmap_mt2_f5510, "Before and during")
    diff_cmap(cmap_mt2_f5510, cmap_mt2_f5692, "During and After")
    diff_cmap(cmap_mt2_f5472, cmap_mt2_f5692, "Before and After")

    # --------------------------------------------------------------------------
    plt.show()



# //////////////////////////////////////////////////////////////////////////////
