from _params import *
from calculations import calc_cluster

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.pyplot import cm

# //////////////////////////////////////////////////////////////////////////////
class Plotter:
    def stylize_ax(self, ax, title, xlabel, ylabel):
        ax.set_title(title, fontdict = dict(fontsize = 20))
        ax.set_xlabel(xlabel, fontdict = dict(fontsize = 16))
        ax.set_ylabel(ylabel, fontdict = dict(fontsize = 16))
        ax.tick_params(labelsize = 12)

    def plot_scatter(self, ax, x, y, colors, marker = '.'):
        line = ax.scatter(x, y, c = colors, marker = marker)
        return line


    def plot_heatmap(self, fig, ax, data, cmap):
        im = ax.imshow(data, cmap)
        colorbar = fig.colorbar(im)
        colorbar.ax.tick_params(labelsize = 12)
        return im, colorbar


    def update_line_data(self, line, x, data):
        line.set_offsets(
            np.append(x, data).reshape((2, x.size)).T
        )


# ////////////////////////////////////////////////////////////////////////////// BSE
class Plotter_BSE(Plotter):
    def vis_BSE(self, bse_naive_0, bse_altern_0, bse_naive_1, bse_altern_1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting BSE for '{title}'...")

        _, ax = plt.subplots()
        self.stylize_ax(ax, title, "Chunk Size", "Standard Deviation")
        ax.plot(bse_naive_0, color = BLUE, linestyle = ':', label = label0 + "_naive")
        ax.plot(bse_altern_0, color = BLUE, linestyle = '-', label = label0 + "_altern")
        ax.plot(bse_naive_1, color = RED, linestyle = ':', label = label1 + "_naive")
        ax.plot(bse_altern_1, color = RED, linestyle = '-', label = label1 + "_altern")
        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// CMAP
class Plotter_CMAP(Plotter): pass


# ////////////////////////////////////////////////////////////////////////////// PCA
class Plotter_PCA(Plotter):
    def update_plot(self, idx):
        self.update_line_data(self.line_2d, self.x, self.pcomps[idx])

    def vis_1pca(self, cumvar, space, pcomps, title = ''):
        print(f">>> Plotting PCA for '{title}'...")

        ##### DATA
        self.x = np.arange(pcomps.shape[0])
        self.pcomps = np.copy(pcomps)

        # print("*****")
        # print("cumvar:", cumvar.shape)
        # print("space:", space.shape)
        # print("pcomps:", pcomps.shape)
        # print("*****")

        ##### FIGURES CREATION
        _, ax_cumvar = plt.subplots()
        self.stylize_ax(ax_cumvar, title, "X", "Cumulative Variance")
        ax_cumvar.set_xlim([0, 100])

        fig = plt.figure()
        fleft, self.fright = fig.subfigures(ncols = 2)

        ax_3d = fleft.add_subplot(projection = "3d")
        ax_3d.set_title(title, fontdict = dict(fontsize = 20))
        self.ax_dict = self.fright.subplot_mosaic("a;a;a;b")

        ##### WIDGETS
        self.slid_pc = Slider(
            ax = self.ax_dict['b'],
            label = "pc", color = "orange",
            valstep = 1, valinit = 0,
            valmin = 0, valmax = pcomps.shape[0] - 1
        )
        self.slid_pc.on_changed(self.update_plot)

        ##### PLOTTING
        ax_cumvar.plot(cumvar)
        self.line_2d = self.ax_dict['a'].scatter(self.x, pcomps[0], marker = '.')

        sc = ax_3d.scatter(space[:,0], space[:,1], space[:,2], c = np.arange(space.shape[0]), s = 20, alpha = 0.5)
        ax_3d.legend(*sc.legend_elements(), loc = "lower left")


    def vis_1pca_plotly(self, cumvar, space, pcomps, title = ''):
        import pandas as pd
        import plotly.express as px

        ##### CUMVAR PLOT
        fg = px.line(
            x = np.arange(cumvar.shape[0]),
            y = cumvar,
            labels = {"x" : "components", "y" : "cumulated variance"},
            range_x = [0, 100]
        )
        fg.show()

        ##### 3D SCATTER
        pca_data = pd.DataFrame(space, columns = ["first_comp", "second_comp", "third_comp"])
        fig = px.scatter_3d(
            pca_data, x = "first_comp", y = "second_comp", z = "third_comp",
            color = pca_data.index.values, width = 900, height = 800
        )
        fig.show()

        ##### 2D PLOTS
        fg0 = px.line(x = np.arange(pcomps.shape[0]), y = pcomps[0])
        fg0.show()

        fg1 = px.line(x = np.arange(pcomps.shape[0]), y = pcomps[1])
        fg1.show()

        fg2 = px.line(x = np.arange(pcomps.shape[0]), y = pcomps[2])
        fg2.show()


# ////////////////////////////////////////////////////////////////////////////// PYINTERAPH
class Plotter_Pyinteraph(Plotter): pass


# ////////////////////////////////////////////////////////////////////////////// RAMA
class Plotter_RAMA(Plotter): pass


# ////////////////////////////////////////////////////////////////////////////// RGYR
class Plotter_RGYR(Plotter):
    def vis_2rgyr(self, rgyr0, rgyr1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting RGYR for '{title}'...")

        frames = np.arange(rgyr0.size)
        _, ax = plt.subplots()
        self.stylize_ax(ax, title, "Frame", "Radius of Gyration")
        ax.scatter(frames, rgyr0, color = BLUE, marker = '+', s = 20, label = label0)
        ax.scatter(frames, rgyr1, color = HALF_RED, marker = 'x', s = 20, label = label1)
        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// RMSD
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_RMSD2D(Plotter):
    def vis_rmsd2d(self, rmsd_mat0, title = ''):
        print(f">>> Plotting RMSD 2D for '{title}'...")

        fig, ax = plt.subplots()
        self.stylize_ax(ax, title, "Frame", "Frame")
        self.plot_heatmap(fig, ax, rmsd_mat0, cmap = "Reds")


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_RMSD1D(Plotter):
    init_ref_frame = 0

    def vis_rmsd1d(self, rmsd_mat0, title = ''):
        print(f">>> Plotting RMSD 1D for '{title}'...")

        ##### DATA
        self.x = np.arange(rmsd_mat0.shape[0])
        self.rmsd_mat0 = np.copy(rmsd_mat0)
        self.colors = np.zeros((rmsd_mat0.shape[0], 4))
        self.update_values(self.init_ref_frame)

        ##### FIGURES CREATION
        self.init_axes()
        self.stylize_ax(self.ax_dict['a'], title, "Frame", "RMSD")

        ##### PLOTTING
        self.init_plot()

        ##### WIDGETS
        self.init_widgets()

    # --------------------------------------------------------------------------
    def update_values(self, ref_frame):
        self.rmsd_arr0 = self.rmsd_mat0[ref_frame]
        self.colors[:] = BLUE
        self.colors[ref_frame] = RED

    def init_axes(self):
        fig = plt.figure(layout = "constrained")
        self.ax_dict = fig.subplot_mosaic("a;a;a;a;a;a;b")

    def init_plot(self):
        self.line_rmsd0 = self.ax_dict['a'].scatter(self.x, self.rmsd_arr0, color = self.colors, marker = '.')

    def init_widgets(self):
        self.init_slider_ref_frame()

    def init_slider_ref_frame(self):
        self.slid_ref_frame = Slider(
            ax = self.ax_dict['b'],
            label = "ref_frame", color = "orange",
            valstep = 1, valinit = self.init_ref_frame,
            valmin = 0, valmax = self.rmsd_mat0.shape[0] - 1,
            valfmt = "%04i"
        )
        self.slid_ref_frame.on_changed(self.update_plot)

    # --------------------------------------------------------------------------
    def update_plot(self, ref_frame):
        self.update_values(ref_frame)
        self.update_line_rmsd()

    def update_line_rmsd(self):
        self.update_line_data(self.line_rmsd0, self.x, self.rmsd_arr0)
        self.line_rmsd0.set_color(self.colors)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_RMSD1D_Compare(Plotter_RMSD1D):
    def vis_rmsd1d(self, rmsd_mat0, rmsd_mat1, title = '', label0 = '', label1 = ''):
        self.rmsd_mat1 = rmsd_mat1
        self.label0 = label0
        self.label1 = label1
        super().vis_rmsd1d(rmsd_mat0, title)
        print(f"...>>> Comparison mode.")

    # --------------------------------------------------------------------------
    def update_values(self, ref_frame):
        super().update_values(ref_frame)
        self.rmsd_arr1 = self.rmsd_mat1[ref_frame]

    def init_plot(self):
        self.line_rmsd0 = self.ax_dict['a'].scatter(self.x, self.rmsd_arr0, color = BLUE, marker = '+', s = 12, label = self.label0)
        self.line_rmsd1 = self.ax_dict['a'].scatter(self.x, self.rmsd_arr1, color = HALF_RED, marker = 'x', s = 12, label = self.label1)
        self.ax_dict['a'].legend()

    def update_line_rmsd(self):
        self.update_line_data(self.line_rmsd0, self.x, self.rmsd_arr0)
        self.update_line_data(self.line_rmsd1, self.x, self.rmsd_arr1)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ CLUSTERING
class Plotter_RMSD1D_Clustering(Plotter_RMSD1D):
    init_clustering_t = 900
    label_criterion = "distance"

    def vis_rmsd1d_clustering(self, rmsd_mat0, Z, color_map, title = ''):
        self.Z = Z
        self.color_map = color_map
        self.highlight_cluster = False
        self.update_values_cluster(self.init_clustering_t)

        super().vis_rmsd1d(rmsd_mat0, title)
        print(f"...>>> Clustering mode.")


    def update_values_cluster(self, t):
        self.cluster_mat = calc_cluster(self.Z, t, label_criterion = self.label_criterion)
        self.n_clusters = np.max(self.cluster_mat)
        self.cluster_colors = cm.__dict__[self.color_map](np.linspace(0, 1, self.n_clusters))

    # --------------------------------------------------------------------------
    def init_axes(self):
        self.fig = plt.figure(layout = "constrained")
        self.ax_dict = self.fig.subplot_mosaic("aaaa;aaaa;aaaa;aaaa;aaaa;aaaa;bbbd;cccd")

    def init_widgets(self):
        super().init_widgets()

        self.slid_clustering_t = Slider(
            ax = self.ax_dict['c'],
            label = "t value", color = "blue",
            valstep = 25, valinit = self.init_clustering_t,
            valmin = 25, valmax = 2000,
            valfmt = "%04i"
        )
        self.slid_clustering_t.on_changed(self.update_clustering_t)

        self.button_toggle_alpha = Button(
            ax = self.ax_dict['d'],
            label = f"{self.n_clusters} clusters"
        )
        self.button_toggle_alpha.on_clicked(self.toggle_alpha)

    # --------------------------------------------------------------------------
    def update_clustering_t(self, clustering_t):
        self.update_values_cluster(clustering_t)
        self.button_toggle_alpha.label.set_text(f"{self.n_clusters} clusters")
        self.update_plot(self.slid_ref_frame.val)

    def toggle_alpha(self, event):
        self.highlight_cluster = not self.highlight_cluster
        self.ax_dict['a'].set_facecolor((0, 0, 0) if self.highlight_cluster else (1, 1, 1))
        self.fig.canvas.draw_idle()
        self.update_plot(self.slid_ref_frame.val)

    def update_values(self, ref_frame):
        self.rmsd_arr0 = self.rmsd_mat0[ref_frame]
        current_cluster = self.cluster_mat[ref_frame]

        for i,color in enumerate(self.cluster_colors[:,:]):
            r,g,b,a = color
            if (i + 1 != current_cluster) and self.highlight_cluster:
                a = .1
            self.colors[self.cluster_mat == i + 1] = r,g,b,a

    # --------------------------------------------------------------------------



# ////////////////////////////////////////////////////////////////////////////// RMSF
class Plotter_RMSF(Plotter):
    def vis_2rmsf(self, rmsf0, rmsf1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting RMSF for '{title}'...")

        frames = np.arange(max(rmsf0.size, rmsf1.size))
        _, ax = plt.subplots()
        self.stylize_ax(ax, title, "CA atom", "RMSF")

        if rmsf0.size == rmsf1.size:
            ax.plot(frames, rmsf0, color = BLUE, linestyle = ':', linewidth = 2, label = label0)
            ax.plot(frames, rmsf1, color = HALF_RED, linestyle = '-', linewidth = 2, label = label1)
        else:
            rmsf_half, rmsf_full = sorted([rmsf0, rmsf1], key = lambda arr: arr.size)
            label_half, label_full = (label0, label1) if rmsf0 is rmsf_half else (label1, label0)

            ax.plot(frames, rmsf_full, color = BLUE, linestyle = ':', linewidth = 2, label = label_full)
            ax.plot(frames[:rmsf_half.size], rmsf_half, color = HALF_RED, linestyle = '-', linewidth = 2, label = label_half)
            ax.plot(frames[rmsf_half.size:], rmsf_half, color = HALF_GREEN, linestyle = '-', linewidth = 2, label = label_half)

        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// SASA
class Plotter_SASA(Plotter): pass




# //////////////////////////////////////////////////////////////////////////////
