from _params import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# //////////////////////////////////////////////////////////////////////////////
class Plotter:
    def __init__(self):
        self.figs = []
        self.axs = []

    def add_fig(self, title, xlabel, ylabel):
        fig, ax = plt.subplots()

        ax.set_title(title, fontdict = dict(fontsize = 20))
        ax.set_xlabel(xlabel, fontdict = dict(fontsize = 16))
        ax.set_ylabel(ylabel, fontdict = dict(fontsize = 16))
        ax.tick_params(labelsize = 12)

        self.figs.append(fig)
        self.axs.append(ax)
        return fig, ax

    def update_line_data(self, line, x, data):
        line.set_offsets(
            np.append(x, data).reshape((2, x.size)).T
        )


# ////////////////////////////////////////////////////////////////////////////// BSE
class Plotter_BSE(Plotter):
    def vis_BSE(self, bse_naive_0, bse_altern_0, bse_naive_1, bse_altern_1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting BSE for '{title}'...")

        _, ax = self.add_fig(title, "Chunk Size", "Standard Deviation")
        ax.plot(bse_naive_0, color = BLUE, linestyle = ':', label = label0 + "_naive")
        ax.plot(bse_altern_0, color = BLUE, linestyle = '-', label = label0 + "_altern")
        ax.plot(bse_naive_1, color = RED, linestyle = ':', label = label1 + "_naive")
        ax.plot(bse_altern_1, color = RED, linestyle = '-', label = label1 + "_altern")
        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// CLUSTERING
class Plotter_Clustering(Plotter): pass


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
        _, ax_cumvar = self.add_fig(title, "X", "Cumulative Variance")
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
        _, ax = self.add_fig(title, "Frame", "Radius of Gyration")
        ax.scatter(frames, rgyr0, color = BLUE, marker = '+', s = 20, label = label0)
        ax.scatter(frames, rgyr1, color = HALF_RED, marker = 'x', s = 20, label = label1)
        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// RMSD
class Plotter_RMSD(Plotter): pass


# ////////////////////////////////////////////////////////////////////////////// RMSF
class Plotter_RMSF(Plotter):
    def vis_2rmsf(self, rmsf0, rmsf1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting RMSF for '{title}'...")

        frames = np.arange(max(rmsf0.size, rmsf1.size))
        _, ax = self.add_fig(title, "CA atom", "RMSF")

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
