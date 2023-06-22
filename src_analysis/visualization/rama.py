from parameters import *
import matplotlib.pyplot as plt
from MDAnalysis.analysis.data.filenames import Rama_ref

# # //////////////////////////////////////////////////////////////////////////////
def plot_rama_ref(**kwargs):
    """Adapted from MDAnalysis/analysis/dihedrals.Ramachandran.plot"""

    fig, ax = plt.subplots()

    ax.set_facecolor((0,0,0,1))
    ax.axis([-180, 180, -180, 180])
    ax.axhline(0, color = 'w', lw = 1)
    ax.axvline(0, color = 'w', lw = 1)
    ax.set(xticks = range(-180, 181, 60), yticks = range(-180, 181, 60))
    degree_formatter = plt.matplotlib.ticker.StrMethodFormatter("{x:g}°")
    ax.xaxis.set_major_formatter(degree_formatter)
    ax.yaxis.set_major_formatter(degree_formatter)

    X, Y = np.meshgrid(
        np.arange(-180, 180, 4),
        np.arange(-180, 180, 4)
    )
    levels = [1, 17, 15000]
    colors = ['#888888', '#AAAAAA']
    ax.contourf(X, Y, np.load(Rama_ref), levels = levels, colors = colors)
    return ax


def vis_2rama(rama0, rama1, title = '', label0 = '', label1 = ''):
    print(f">>> Plotting RAMA for '{title}'...")
    ax = plot_rama_ref()

    ax.set_title(title, fontdict = dict(fontsize = 20))
    ax.set_xlabel(r"$\phi$", fontdict = dict(fontsize = 16))
    ax.set_ylabel(r"$\psi$", fontdict = dict(fontsize = 16))
    ax.tick_params(labelsize = 12)


    a = rama0[0].reshape(rama0.shape[1], 2)
    ax.scatter(a[:,0], a[:,1], color = RED,  marker = 'x', s = 20, label = label0)

    b = rama1[0].reshape(rama1.shape[1], 2)
    ax.scatter(b[:,0], b[:,1], color = BLUE, marker = '+', s = 20, label = label1)

    ax.legend()


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        PATH_RAMA0 = DIR_DA_GENERAL / f"{run0}-rama.npy"
        PATH_RAMA1 = DIR_DA_GENERAL / f"{run1}-rama.npy"

        vis_2rama(np.load(PATH_RAMA0), np.load(PATH_RAMA1), run_preffix, run0, run1)

    plt.show()

# //////////////////////////////////////////////////////////////////////////////
