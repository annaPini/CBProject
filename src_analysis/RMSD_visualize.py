from parameters import *
from calculate_specific import calc_cluster

from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.pyplot import cm

# //////////////////////////////////////////////////////////////////////////////
class RMSD_Plotter:
    def __init__(self, rmsd_mat, title = ''):
        print(f">>> Plotting RMSD for '{title}'...")
        self.rmsd_mat = rmsd_mat

        self.fig, self.ax = plt.subplots()
        self.ax.set_title(title, fontdict = dict(fontsize = 20))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RMSD_2D(RMSD_Plotter):
    def __init__(self, rmsd_mat, title = ''):
        super().__init__(rmsd_mat, title)

        self.fig.subplots_adjust(bottom = 0.1, top = 0.9, left = 0.15, right = 0.95)
        self.ax.set_xlabel("Frame", fontdict = dict(fontsize = 16))
        self.ax.set_ylabel("Frame", fontdict = dict(fontsize = 16))

        colorbar = self.fig.colorbar( self.ax.imshow(self.rmsd_mat, cmap = "Reds") )

        self.ax.tick_params(labelsize = 12)
        colorbar.ax.tick_params(labelsize = 12)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RMSD_1D(RMSD_Plotter):
    init_ref_frame = 0

    def __init__(self, rmsd_mat, title = ''):
        super().__init__(rmsd_mat, title)

        self.x = np.arange(self.rmsd_mat.shape[0])
        self.colors = np.zeros((self.rmsd_mat.shape[0], 4))

        self.fig.subplots_adjust(bottom = 0.25, top = 0.9)
        self.ax.set_xlabel("Frame", fontdict = dict(fontsize = 16))
        self.ax.set_ylabel("RMSD", fontdict = dict(fontsize = 16))
        self.ax.tick_params(labelsize = 12)

        self.init_gui()

    # --------------------------------------------------------------------------
    def init_gui(self):
        self.set_color(self.init_ref_frame)
        self.init_plot()
        self.init_slider_ref_frame()

    def update_plot(self, val = None):
        self.set_color(self.slid_ref_frame.val)
        self.update_line_rmsd()

    def set_color(self, ref_frame):
        self.rmsd_arr = self.rmsd_mat[ref_frame]
        self.colors[:] = BLUE
        self.colors[ref_frame] = RED

    # --------------------------------------------------------------------------
    def init_plot(self):
        self.line_rmsd = self.ax.scatter(self.x, self.rmsd_arr, c = self.colors, marker = '.')

    def init_slider_ref_frame(self):
        self.slid_ref_frame = Slider(
            ax = plt.axes((.15, .10, .5, .03)),
            label = "ref_frame", color = "orange" ,
            valstep = 1, valinit = self.init_ref_frame,
            valmin = 0, valmax = self.rmsd_mat.shape[0] - 1,
        )
        self.slid_ref_frame.on_changed(self.update_plot)

    def update_line_rmsd(self):
        self.line_rmsd.set_offsets(
            np.append(self.x, self.rmsd_arr).reshape((2, self.x.size)).T
        )
        self.line_rmsd.set_color(self.colors)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RMSD_1D_WAD(RMSD_1D):
    init_rmsd_treshold = 2

    def init_gui(self):
        self.y = np.zeros(self.rmsd_mat.shape[0])
        self.set_color(self.init_ref_frame, self.init_rmsd_treshold)
        self.init_plot()
        self.init_slider_ref_frame()
        self.init_line_treshold()

    def update_plot(self, val = None):
        self.set_color(self.slid_ref_frame.val, self.slid_rmsd_treshold.val)
        self.update_line_rmsd()
        self.update_line_treshold()

    def set_color(self, ref_frame, rmsd_treshold):
        self.rmsd_arr = self.rmsd_mat[ref_frame]

        self.colors[:] = BLUE
        self.colors[self.rmsd_arr <= rmsd_treshold] = GREEN
        self.colors[ref_frame] = RED

        self.y[:] = rmsd_treshold

    # --------------------------------------------------------------------------
    def init_line_treshold(self):
        self.line_treshold = self.ax.plot(self.x, self.y)
        self.slid_rmsd_treshold = Slider(
            ax = plt.axes((.15, .05, .5, .03)),
            label = "rmsd_treshold", color = "blue",
            valstep = .05, valinit = self.init_rmsd_treshold,
            valmin = 0, valmax = np.max(self.rmsd_mat, title = ''),
        )
        self.slid_rmsd_treshold.on_changed(self.update_plot)

        self.button_save = Button(ax = plt.axes((.75, .05, .2, .1)), label = "Select WA frames")
        self.button_save.on_clicked(self.save_selected_frames)

    def update_line_treshold(self):
        self.line_treshold[0].set_data(self.x, self.y)

    # --------------------------------------------------------------------------
    def save_selected_frames(self, val = None):
        ref_frame = self.slid_ref_frame.val
        rmsd_treshold = self.slid_rmsd_treshold.val

        self.rmsd_arr = self.rmsd_mat[ref_frame]
        frames = [int(f) for f in self.x[self.rmsd_arr <= rmsd_treshold]]

        info = Info(DIR_DA_WAD / f"{CURRENT_RUN}-{ref_frame}-info.json")

        info.update(
            rsmd_treshold = rmsd_treshold,
            ref_frame = ref_frame,
            ref_frame_index = frames.index(ref_frame),
            frames = frames,
        )

        self.button_save.label.set_text(f"Saved {len(frames)} frames")




# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RMSD_1D_Clustering(RMSD_1D):
    init_clustering_t = 900
    label_criterion = "distance"

    def __init__(self, rmsd_mat, Z, color_map):
        self.Z = Z
        self.color_map = color_map
        self.highlight_cluster = False

        self.set_cluster_colors(self.init_clustering_t)

        super().__init__(rmsd_mat, title)

    # --------------------------------------------------------------------------
    def init_gui(self):
        super().init_gui()
        self.init_gui_clustering()

    def set_color(self, ref_frame):
        self.rmsd_arr = self.rmsd_mat[ref_frame]
        current_cluster = self.cluster_mat[ref_frame]

        for i,color in enumerate(self.cluster_colors[:,:]):
            r,g,b,a = color
            if (i + 1 != current_cluster) and self.highlight_cluster:
                a = .1
            self.colors[self.cluster_mat == i + 1] = r,g,b,a

    # --------------------------------------------------------------------------
    def init_gui_clustering(self):
        self.slid_clustering_t = Slider(
            ax = plt.axes((.15, .05, .5, .03)),
            label = "t value", color = "blue",
            valstep = 25, valinit = self.init_clustering_t,
            valmin = 25, valmax = 2000,
        )
        self.slid_clustering_t.on_changed(self.update_clustering_t)

        self.button_toggle_alpha = Button(ax = plt.axes((.75, .05, .2, .1)), label = f"{self.n_clusters} clusters")
        self.button_toggle_alpha.on_clicked(self.toggle_alpha)

    def set_cluster_colors(self, t):
        self.cluster_mat = calc_cluster(self.Z, t, label_criterion = self.label_criterion)
        self.n_clusters = np.max(self.cluster_mat)
        self.cluster_colors = cm.__dict__[self.color_map](np.linspace(0, 1, self.n_clusters))


    # --------------------------------------------------------------------------
    def update_clustering_t(self, vale = None):
        self.set_cluster_colors(self.slid_clustering_t.val)
        self.button_toggle_alpha.label.set_text(f"{self.n_clusters} clusters")
        self.update_plot()

    def toggle_alpha(self, val = None):
        self.highlight_cluster = not self.highlight_cluster
        self.ax.set_facecolor((0, 0, 0) if self.highlight_cluster else (1, 1, 1))
        self.update_plot()


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RMSD_1D_Compare(RMSD_1D):
    def __init__(self, rmsd_mat_0, rmsd_mat_1, title, label0, label1):
        self.rmsd_mat1 = rmsd_mat_1
        self.label0 = label0
        self.label1 = label1
        super().__init__(rmsd_mat_0, title)

    # --------------------------------------------------------------------------
    def set_color(self, ref_frame):
        super().set_color(ref_frame)
        self.rmsd_arr1 = self.rmsd_mat1[ref_frame]

    def init_plot(self):
        self.line_rmsd0 = self.ax.scatter(self.x, self.rmsd_arr, color = BLUE, marker = '+', s = 12, label = self.label0)
        self.line_rmsd1 = self.ax.scatter(self.x, self.rmsd_arr1, color = HALF_RED, marker = 'x', s = 12, label = self.label1)
        self.ax.legend()

    def update_line_rmsd(self):
        self.line_rmsd0.set_offsets(
            np.append(self.x, self.rmsd_arr).reshape((2, self.x.size)).T
        )
        self.line_rmsd1.set_offsets(
            np.append(self.x, self.rmsd_arr1).reshape((2, self.x.size)).T
        )


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# //////////////////////////////////////////////////////////////////////////////


if __name__ == "__main__":
    # COLOR_MAP = "rainbow"
    # COLOR_MAP = "prism"
    COLOR_MAP = "hsv"

    rmsds = []

    ############################################################################
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        PATH_RMSD0 = DIR_DA_GENERAL / f"{run0}-rmsd.npy"
        PATH_RMSD1 = DIR_DA_GENERAL / f"{run1}-rmsd.npy"
        rmsd_mat0 = np.load(PATH_RMSD0)
        rmsd_mat1 = np.load(PATH_RMSD1)

        # rmsds.append(RMSD_2D(rmsd_mat0, title = run0))
        # rmsds.append(RMSD_2D(rmsd_mat1, title = run1))
        rmsds.append(RMSD_1D_Compare(rmsd_mat0, rmsd_mat1, run_preffix, run0, run1))


    ############################################################################

    # run0 = "mt2_rep0"
    # run1 = "mt1_rep1"
    #
    # PATH_RMSD0 = DIR_DA_GENERAL / f"{run0}-rmsd.npy"
    # PATH_RMSD1 = DIR_DA_GENERAL / f"{run1}-rmsd.npy"
    #
    # PATH_LINK = DIR_DA_SPECIFIC / f"{run0}-link.npy"
    # Z = np.load(PATH_LINK)
    #
    #
    # rmsd_mat0 = np.load(PATH_RMSD0)
    # rmsd_mat1 = np.load(PATH_RMSD1)
    #
    # # rmsds.append(RMSD_2D(rmsd_mat0))
    # # rmsds.append(RMSD_1D(rmsd_mat0))
    # # rmsds.append(RMSD_1D_WAD(rmsd_mat0))
    # # rmsds.append(RMSD_1D_Clustering(rmsd_mat0, Z, COLOR_MAP))
    # rmsds.append(RMSD_1D_Compare(rmsd_mat0, rmsd_mat1))
    #

    ############################################################################

    plt.show()

# //////////////////////////////////////////////////////////////////////////////

# https://matplotlib.org/stable/tutorials/colors/colormaps.html
# ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'ColormapRegistry', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'Mapping', 'MutableMapping', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'ScalarMappable', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', '_DeprecatedCmapDictWrapper', '_LUTSIZE', '__builtin_cmaps', '__builtins__', '__cached__', '__doc__', '__file__', '__getattr__', '__loader__', '__name__', '__package__', '__spec__', '_api', '_cmap_registry', '_colormaps', '_gen_cmap_registry', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cbook', 'cividis', 'cividis_r', 'cmap_d', 'cmaps_listed', 'colors', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'datad', 'flag', 'flag_r', 'get_cmap', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'ma', 'magma', 'magma_r', 'mpl', 'nipy_spectral', 'nipy_spectral_r', 'np', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'register_cmap', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'unregister_cmap', 'viridis', 'viridis_r', 'winter', 'winter_r']
