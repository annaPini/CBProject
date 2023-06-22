from _params import *
import numpy as np
import matplotlib.pyplot as plt

# //////////////////////////////////////////////////////////////////////////////
class Plotter:
    def __init__(self):
        self.figs, self.axs = [], []

    def add_fig(self, title, xlabel, ylabel):
        fig, ax = plt.subplots()

        ax.set_title(title, fontdict = dict(fontsize = 20))
        ax.set_xlabel(xlabel, fontdict = dict(fontsize = 16))
        ax.set_ylabel(ylabel, fontdict = dict(fontsize = 16))
        ax.tick_params(labelsize = 12)

        self.figs.append(fig)
        self.axs.append(ax)
        return fig, ax


# ////////////////////////////////////////////////////////////////////////////// BSE
class Plotter_BSE(Plotter):
    def vis_BSE(self, bse_naive, bse_altern, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting BSE for '{title}'...")

        _, ax = self.add_fig(title, "Block", "SE")
        ax.plot(bse_naive, label = label0 + "_naive")
        ax.plot(bse_altern, label = label1 + "_altern")
        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// CLUSTERING
class Plotter_Clustering(Plotter):
    pass


# ////////////////////////////////////////////////////////////////////////////// CMAP
class Plotter_CMAP(Plotter):
    pass


# ////////////////////////////////////////////////////////////////////////////// PCA
class Plotter_PCA(Plotter):
    pass


# ////////////////////////////////////////////////////////////////////////////// PYINTERAPH
class Plotter_Pyinteraph(Plotter):
    pass


# ////////////////////////////////////////////////////////////////////////////// RAMA
class Plotter_RAMA(Plotter):
    pass


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
class Plotter_RMSD(Plotter):
    pass


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
class Plotter_SASA(Plotter):
    pass


# //////////////////////////////////////////////////////////////////////////////
