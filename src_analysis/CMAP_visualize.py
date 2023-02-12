from parameters import *
from calculate_general import calc_cmap
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# //////////////////////////////////////////////////////////////////////////////
class CMAPPlotter:
    def __init__(self, coords, title = ''):
        print(f">>> Plotting CMAP for '{title}'...")
        self.coords = coords
        self.fig, self.ax = plt.subplots()

        self.ax.set_title(title, fontdict = dict(fontsize = 20))
        self.ax.set_xlabel("CA atom", fontdict = dict(fontsize = 16))
        self.ax.set_ylabel("CA atom", fontdict = dict(fontsize = 16))
        self.ax.tick_params(labelsize = 12)

    def vis_cmap(self, cmap, color_method):
        self.im = self.ax.imshow(cmap, cmap = color_method)
        self.colorbar = self.fig.colorbar(self.im)
        self.colorbar.ax.tick_params(labelsize = 12)

    def base_cmap(self, frame):
        self.vis_cmap(-calc_cmap(self.coords[frame]), "Reds")

    def diff_cmap(self, frame0, frame1, color_method = "Blues"):
        cmap0 = calc_cmap(self.coords[frame0])
        cmap1 = calc_cmap(self.coords[frame1])
        self.vis_cmap(abs(cmap1 - cmap0), color_method)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CMAP_Basic(CMAPPlotter):
    def __init(self, coords, title = ''):
        super().__init__(coords, title)
        self.fig.subplots_adjust(bottom = 0.1, top = 0.9, left = 0.15, right = 0.95)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CMAP_Interactive(CMAPPlotter):
    def __init__(self, coords, title = '', min_frame = 0, max_frame = 6000):
        self.min_frame = min_frame
        self.max_frame = max_frame

        super().__init__(coords, title)

        self.fig.subplots_adjust(bottom = 0.25, top = 0.9, left = 0.15, right = 0.8)
        self.init_plot()

    # --------------------------------------------------------------------------
    def init_plot(self):
        self.base_cmap(self.min_frame)
        self.init_slider_frame0()

    # --------------------------------------------------------------------------
    def init_slider_frame0(self):
        self.slid_frame0 = Slider(
            ax = plt.axes((.15, .10, .5, .03)),
            label = "frame", color = "orange" ,
            valstep = 1, valinit = self.min_frame,
            valmin = self.min_frame, valmax = self.max_frame,
        )
        self.slid_frame0.on_changed(self.update_cmap)

    def update_cmap(self, val):
        frame = self.slid_frame0.val
        cmap = calc_cmap(self.coords[frame])

        self.im.set_data(-cmap)
        self.im.set_clim(vmin = np.min(-cmap), vmax = 0)
        self.colorbar.update_normal(self.im)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CMAP_Diff_1mol_2frames(CMAP_Interactive):
    def init_plot(self):
        self.diff_cmap(self.min_frame, self.min_frame + 1, "seismic")
        self.init_slider_frame0()
        self.init_slider_frame1()

    # --------------------------------------------------------------------------
    def init_slider_frame1(self):
        self.slid_frame1 = Slider(
            ax = plt.axes((.15, .05, .5, .03)),
            label = "frame", color = "blue",
            valstep = 1, valinit = self.min_frame + 1,
            valmin = self.min_frame, valmax = self.max_frame,
        )
        self.slid_frame1.on_changed(self.update_cmap)

    def update_cmap(self, val):
        frame0 = self.slid_frame0.val
        frame1 = self.slid_frame1.val
        # diff = calc_cmap(self.coords[frame1]) - calc_cmap(self.coords[frame0])
        diff = abs(calc_cmap(self.coords[frame1]) - calc_cmap(self.coords[frame0]))

        self.im.set_data(diff)
        # self.im.set_clim(vmin = np.min(diff), vmax = np.max(diff))
        self.im.set_clim(0, vmax = np.max(diff))
        self.colorbar.update_normal(self.im)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CMAP_Diff_2mol_1frame(CMAP_Interactive):
    def __init__(self, coords0, coords1, title = '', min_frame = 0, max_frame = 6000):
        self.coords1 = coords1
        super().__init__(coords0, title, min_frame, max_frame)

    def init_plot(self):
        cmap0 = calc_cmap(self.coords[self.min_frame])
        cmap1 = calc_cmap(self.coords1[self.min_frame])
        self.vis_cmap(abs(cmap1 - cmap0), "seismic")

        self.init_slider_frame0()

    # --------------------------------------------------------------------------
    def update_cmap(self, val):
        frame = self.slid_frame0.val
        diff = abs(calc_cmap(self.coords1[frame]) - calc_cmap(self.coords[frame]))

        self.im.set_data(diff)
        self.im.set_clim(0, vmax = np.max(diff))
        self.colorbar.update_normal(self.im)


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    coords_mt1 = np.load(DIR_DA_TRAJECTORIES / "mt1_rep0-coords.npy")
    coords_mt2 = np.load(DIR_DA_TRAJECTORIES / "mt2_rep0-coords.npy")
    coords_wt2 = np.load(DIR_DA_TRAJECTORIES / "wt2_rep0-coords.npy")

    # --------------------------------------------------------------------------
    # CMAP_Basic(coords_mt1, "mt1_rep0").base_cmap(frame = 0)
    # CMAP_Basic(coords_mt2, "mt2_rep0").base_cmap(frame = 0)

    # --------------------------------------------------------------------------
    # frame_before = 5472
    # frame_during = 5510
    # frame_after = 5692
    #
    # CMAP_Basic(coords_mt2, "Before and during").diff_cmap(frame_before, frame_during)
    # CMAP_Basic(coords_mt2, "During and After").diff_cmap(frame_during, frame_after)
    # CMAP_Basic(coords_mt2, "Before and After").diff_cmap(frame_before, frame_after)

    # --------------------------------------------------------------------------
    # cmi = CMAP_Interactive(coords_mt2, "mt2_rep0")
    # cmdiff = CMAP_Diff_1mol_2frames(coords_mt2, "mt2_rep0")
    cmdiff_a = CMAP_Diff_1mol_2frames(coords_wt2, "mt2_rep0", min_frame = 4000, max_frame = 6000)
    cmdiff_b = CMAP_Diff_2mol_1frame(coords_mt2, coords_wt2, "mt2 vs wt2", min_frame = 4000, max_frame = 6000)

    # --------------------------------------------------------------------------
    plt.show()


# //////////////////////////////////////////////////////////////////////////////
