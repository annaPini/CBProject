from parameters import *
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.pyplot import cm
# import matplotlib as mpl
# cmap = mpl.colormaps['viridis']

# print(dir(cm))
# print(*mpl.colormaps.keys(), sep = ' ')
# exit()

# //////////////////////////////////////////////////////////////////////////////
class RMSDPlotter:
    def __init__(self, rmsd_mat):
        print(f">>> Loading RMSD data...")
        self.rmsd_mat = rmsd_mat


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RMSD_2D(RMSDPlotter):
    def __init__(self, rmsd_mat):
        super().__init__(rmsd_mat)
        plt.colorbar(plt.imshow(self.rmsd_mat, vmax = 5))


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RMSD_1D(RMSDPlotter):
    init_ref_frame = 0

    def __init__(self, rmsd_mat):
        super().__init__(rmsd_mat)

        self.x = np.arange(self.rmsd_mat.shape[0])
        self.colors = np.zeros((self.rmsd_mat.shape[0], 4))

        self.fig, self.ax = plt.subplots()
        self.fig.subplots_adjust(bottom = 0.25, top = 0.9)
        # self.ax.set_facecolor((0, 0, 0))
        # plt.xlabel("Frame")
        # plt.ylabel("RMSD")

        self.init_plots()

    # --------------------------------------------------------------------------
    def init_plots(self):
        self.set_color(self.init_ref_frame)
        self.init_line_rmsd()

    def update_plot(self, val = None):
        self.set_color(self.slid_ref_frame.val)
        self.update_line_rmsd()

    def set_color(self, ref_frame):
        self.rmsd_arr = self.rmsd_mat[ref_frame]
        self.colors[:] = BLUE
        self.colors[ref_frame] = RED

    # --------------------------------------------------------------------------
    def init_line_rmsd(self):
        self.line_rmsd = self.ax.scatter(self.x, self.rmsd_arr, c = self.colors, marker = '.')
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

    def init_plots(self):
        self.y = np.zeros(self.rmsd_mat.shape[0])
        self.set_color(self.init_ref_frame, self.init_rmsd_treshold)
        self.init_line_rmsd()
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
            valmin = 0, valmax = np.max(self.rmsd_mat),
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
    def __init__(self, rmsd_mat, cluster_mat):
        self.cluster_mat = cluster_mat

        n = np.max(self.cluster_mat)
        print(f">>> Coloring {n} clusters...")

        # self.cluster_colors = cm.rainbow(np.linspace(0, 1, n))
        # self.cluster_colors = cm.prism(np.linspace(0, 1, n))
        self.cluster_colors = cm.hsv(np.linspace(0, 1, n))
        # np.random.shuffle(self.cluster_colors)

        self.highlight_cluster = False

        super().__init__(rmsd_mat)

        self.button_toggle_alpha = Button(ax = plt.axes((.75, .05, .2, .1)), label = "Highlight cluster")
        self.button_toggle_alpha.on_clicked(self.toggle_alpha)


    def set_color(self, ref_frame):
        self.rmsd_arr = self.rmsd_mat[ref_frame]
        current_cluster = self.cluster_mat[ref_frame]

        for i,color in enumerate(self.cluster_colors[:,:]):
            r,g,b,a = color
            if (i + 1 != current_cluster) and self.highlight_cluster:
                a = .1
            self.colors[self.cluster_mat == i + 1] = r,g,b,a

    def toggle_alpha(self, val = None):
        self.highlight_cluster = not self.highlight_cluster
        self.update_plot()

        # self.ax.set_facecolor((0, 0, 0) if self.highlight_cluster else (1,1,1))
        if self.highlight_cluster:
            self.ax.set_facecolor((0, 0, 0))
            self.button_toggle_alpha.label.set_text("Show all")
        else:
            self.ax.set_facecolor((1, 1, 1))
            self.button_toggle_alpha.label.set_text("Highlight cluster")



    # --------------------------------------------------------------------------
    # this time colors don't change when updating ref_frame
    # def update_plot(self, val = None):
    #     ref_frame = self.slid_ref_frame.val
    #     self.rmsd_arr = self.rmsd_mat[ref_frame]
    #
    #     self.line_rmsd.set_offsets(
    #         np.append(self.x, self.rmsd_arr).reshape((2, self.x.size)).T
    #     )


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":

    rmsds = []

    for run in ["wt2_rep0", "wt2_rep1"]:
    # for run in RUNS[1:2]:
        PATH_RMSD = DIR_DA_GENERAL / f"{run}-rmsd.npy"
        PATH_CLUSTER = DIR_DA_GENERAL / f"{run}-cluster.npy"

        rmsd_mat = np.load(PATH_RMSD)
        cluster_mat = np.load(PATH_CLUSTER)

        # print(cluster_mat)
        # print(cluster_mat.shape)

        ### should assign instances to variables, otherwise sliders might become irresponsive! (garbage collecion?)
        # rmsds.append(RMSD_2D(rmsd_mat))
        rmsds.append(RMSD_1D(rmsd_mat))
        # rmsds.append(RMSD_1D_WAD(rmsd_mat))
        # rmsds.append(RMSD_1D_Clustering(rmsd_mat, cluster_mat))

    plt.show()

# //////////////////////////////////////////////////////////////////////////////

# https://matplotlib.org/stable/tutorials/colors/colormaps.html
# ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'ColormapRegistry', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'Mapping', 'MutableMapping', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'ScalarMappable', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', '_DeprecatedCmapDictWrapper', '_LUTSIZE', '__builtin_cmaps', '__builtins__', '__cached__', '__doc__', '__file__', '__getattr__', '__loader__', '__name__', '__package__', '__spec__', '_api', '_cmap_registry', '_colormaps', '_gen_cmap_registry', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cbook', 'cividis', 'cividis_r', 'cmap_d', 'cmaps_listed', 'colors', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'datad', 'flag', 'flag_r', 'get_cmap', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'ma', 'magma', 'magma_r', 'mpl', 'nipy_spectral', 'nipy_spectral_r', 'np', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'register_cmap', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'unregister_cmap', 'viridis', 'viridis_r', 'winter', 'winter_r']
