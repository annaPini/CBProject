from parameters import *
from calculate_general import calc_cmap, calc_cmap_AS
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


# //////////////////////////////////////////////////////////////////////////////
class CMAP_Plotter:
    do_abs_on_diff = True

    def __init__(self, coords, title = ''):
        print(f">>> Plotting CMAP for '{title}'...")
        self.coords = coords
        self.fig, self.ax = plt.subplots()

        self.ax.set_title(title, fontdict = dict(fontsize = 20))
        self.ax.set_xlabel("CA atom", fontdict = dict(fontsize = 16))
        self.ax.set_ylabel("CA atom", fontdict = dict(fontsize = 16))
        self.ax.tick_params(labelsize = 12)

        self.do_abs = abs if self.do_abs_on_diff else (lambda x:x)

    def vis_cmap(self, cmap, color_method):
        self.im = self.ax.imshow(cmap, cmap = color_method)
        self.colorbar = self.fig.colorbar(self.im)
        self.colorbar.ax.tick_params(labelsize = 12)

    def base_cmap(self, frame):
        self.vis_cmap(-calc_cmap(self.coords[frame]), "Reds")

    def diff_cmap(self, frame0, frame1, color_method = "Blues"):
        cmap0 = calc_cmap(self.coords[frame0])
        cmap1 = calc_cmap(self.coords[frame1])
        self.vis_cmap(self.do_abs(cmap1 - cmap0), color_method)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CMAP_Basic(CMAP_Plotter):
    def __init__(self, coords, title = ''):
        super().__init__(coords, title)
        self.fig.subplots_adjust(bottom = 0.1, top = 0.9, left = 0.15, right = 0.95)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Dynamic_CMAP(CMAP_Plotter):
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
class DCMD_1mol_2frames(Dynamic_CMAP):
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
        diff = self.do_abs(calc_cmap(self.coords[frame1]) - calc_cmap(self.coords[frame0]))

        self.im.set_data(diff)
        self.im.set_clim(vmin = np.min(diff), vmax = np.max(diff))
        self.colorbar.update_normal(self.im)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class DCMD_2mol_1frame(Dynamic_CMAP):
    def __init__(self, coords0, coords1, title = '', min_frame = 0, max_frame = 6000):
        self.coords1 = coords1
        super().__init__(coords0, title, min_frame, max_frame)

    def init_plot(self):
        cmap0 = calc_cmap(self.coords[self.min_frame])
        cmap1 = calc_cmap(self.coords1[self.min_frame])
        self.vis_cmap(self.do_abs(cmap1 - cmap0), "seismic")

        self.init_slider_frame0()

    # --------------------------------------------------------------------------
    def update_cmap(self, val):
        frame = self.slid_frame0.val
        diff = self.do_abs(calc_cmap(self.coords1[frame]) - calc_cmap(self.coords[frame]))

        self.im.set_data(diff)
        self.im.set_clim(0, vmax = np.max(diff))
        self.colorbar.update_normal(self.im)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plot_cmap_AS(coords_AS):

    sns.set_palette("bwr")

    fig, ax = plt.subplots()

    interactions = ["AB", "CD", "EF", "FG", "EG", "GH", "IJ"]

    df = pd.DataFrame(columns = ["system", "interaction", "distance"])
    print (f">>> Calculating CMAPs...")

    for system,coords in coords_AS.items():

        atoms = {atom : coords[:,i,:] for i,atom in enumerate("ABCDEFGHIJ")}

        for atom0, atom1 in interactions:
            new_df = pd.DataFrame({"distance" : np.diag(calc_cmap_AS(atoms[atom0], atoms[atom1]))})
            new_df["interaction"] = atom0 + atom1
            new_df["system"] = system
            df = df.append(new_df)


    print (f">>> Plotting...")
    # # sns.catplot(data = df, x = "interaction", y = "distance", hue = "system", kind = "boxen")
    # sns.boxplot(data = df, x = "interaction", y = "distance", hue = "system")
    # sns.catplot(data = df, x = "interaction", y = "distance", hue = "system", kind = "point")
    # plt.show(block = True)

    # sns.catplot(data = df, x = "interaction", y = "distance", hue = "system", kind = "boxen", ax = ax)
    # sns.catplot(data = df, x = "interaction", y = "distance", hue = "system", kind = "point", ax = ax)
    #

    # sns.boxplot(data = df, x = "interaction", y = "distance", hue = "system", ax = ax)
    sns.pointplot(data = df, x = "interaction", y = "distance", hue = "system", ax = ax)


    plt.show(block = True)




# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":

    rep = 1

    coords_mt1 = np.load(DIR_DA_TRAJECTORIES / f"mt1_rep{rep}-coords.npy")
    coords_mt2 = np.load(DIR_DA_TRAJECTORIES / f"mt2_rep{rep}-coords.npy")
    coords_wt2 = np.load(DIR_DA_TRAJECTORIES / f"wt2_rep{rep}-coords.npy")

    active_sites = {
        "mt1_as0" : np.load(DIR_DA_SPECIFIC / f"mt1_rep{rep}-AS_coords.npy"),
        "mt2_as0" : np.load(DIR_DA_SPECIFIC / f"mt2_rep{rep}-AS_coords.npy"),
        "mt2_as1" : np.load(DIR_DA_SPECIFIC / f"mt2_rep{rep}-AS_coords1.npy"),
        "wt1_as0" : np.load(DIR_DA_SPECIFIC / f"wt1_rep{rep}-AS_coords.npy"),
        "wt2_as0" : np.load(DIR_DA_SPECIFIC / f"wt2_rep{rep}-AS_coords.npy"),
        "wt2_as1" : np.load(DIR_DA_SPECIFIC / f"wt2_rep{rep}-AS_coords1.npy"),
    }

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
    # # cmi = Dynamic_CMAP(coords_mt2, "mt2_rep0")
    # # cmdiff = DCMD_1mol_2frames(coords_mt2, "mt2_rep0")
    # cmdiff_a = DCMD_1mol_2frames(coords_wt2, "mt2_rep0", min_frame = 4000, max_frame = 6000)
    # cmdiff_b = DCMD_2mol_1frame(coords_mt2, coords_wt2, "mt2 vs wt2", min_frame = 4000, max_frame = 6000)

    # --------------------------------------------------------------------------
    # plot_cmap_AS(active_sites)

    # --------------------------------------------------------------------------
    plt.show()

# //////////////////////////////////////////////////////////////////////////////
