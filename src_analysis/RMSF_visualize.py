from parameters import *
import matplotlib.pyplot as plt

# //////////////////////////////////////////////////////////////////////////////
def vis_2rmsf(rmsf0, rmsf1, title = '', label0 = '', label1 = ''):
    print(f">>> Plotting RMSF for '{title}'...")

    frames = np.arange(max(rmsf0.size, rmsf1.size))
    fig, ax = plt.subplots()

    # ax.set_facecolor((0, 0, 0))
    ax.set_title(title, fontdict = dict(fontsize = 20))
    ax.set_xlabel("CA atom", fontdict = dict(fontsize = 16))
    ax.set_ylabel("RMSF", fontdict = dict(fontsize = 16))
    ax.tick_params(labelsize = 12)

    if rmsf0.size != rmsf1.size:
        rmsf_half, rmsf_full = sorted([rmsf0, rmsf1], key = lambda arr: arr.size)
        label_half, label_full = (label0, label1) if rmsf0 is rmsf_half else (label1, label0)

        ax.plot(frames, rmsf_full, color = BLUE, linestyle = ':', linewidth = 2, label = label_full)
        ax.plot(frames[:rmsf_half.size], rmsf_half, color = HALF_RED, linestyle = '-', linewidth = 2, label = label_half)
        ax.plot(frames[rmsf_half.size:], rmsf_half, color = HALF_GREEN, linestyle = '-', linewidth = 2, label = label_half)

    else:
        ax.plot(frames, rmsf0, color = BLUE, linestyle = ':', linewidth = 2, label = label0)
        ax.plot(frames, rmsf1, color = HALF_RED, linestyle = '-', linewidth = 2, label = label1)

    ax.legend()

# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        PATH_RMSF0 = DIR_DA_GENERAL / f"{run0}-rmsf.npy"
        PATH_RMSF1 = DIR_DA_GENERAL / f"{run1}-rmsf.npy"

        rmsf0 = np.load(PATH_RMSF0)
        rmsf1 = np.load(PATH_RMSF1)

        vis_2rmsf(rmsf0, rmsf1, run_preffix, run0, run1)

    comparisons = [
        ("mt2_rep0", "wt2_rep0"),
        ("mt1_rep0", "mt2_rep0"),
        ("wt1_rep0", "wt2_rep0"),
    ]

    for run0, run1 in comparisons:
        PATH_RMSF0 = DIR_DA_GENERAL / f"{run0}-rmsf.npy"
        PATH_RMSF1 = DIR_DA_GENERAL / f"{run1}-rmsf.npy"

        rmsf0 = np.load(PATH_RMSF0)
        rmsf1 = np.load(PATH_RMSF1)

        vis_2rmsf(rmsf0, rmsf1, f"{run0} vs {run1}", run0, run1)

    plt.show()

# //////////////////////////////////////////////////////////////////////////////
