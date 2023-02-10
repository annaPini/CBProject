from parameters import *
import matplotlib.pyplot as plt

# //////////////////////////////////////////////////////////////////////////////
def vis_2rgyr(rgyr0, rgyr1, title = '', label0 = '', label1 = ''):
    print(f">>> Plotting RGYR for '{title}'...")

    frames = np.arange(rgyr0.size)
    fig, ax = plt.subplots()

    ax.set_title(title, fontdict = dict(fontsize = 20))
    ax.set_xlabel("Frame", fontdict = dict(fontsize = 16))
    ax.set_ylabel("Radius of Gyration", fontdict = dict(fontsize = 16))
    ax.tick_params(labelsize = 12)

    ax.scatter(frames, rgyr0, color = BLUE, marker = '+', s = 20, label = label0)
    ax.scatter(frames, rgyr1, color = HALF_RED, marker = 'x', s = 20, label = label1)

    ax.legend()

# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        PATH_RGYR0 = DIR_DA_GENERAL / f"{run0}-rgyr.npy"
        PATH_RGYR1 = DIR_DA_GENERAL / f"{run1}-rgyr.npy"

        rgyr0 = np.load(PATH_RGYR0)
        rgyr1 = np.load(PATH_RGYR1)

        vis_2rgyr(rgyr0, rgyr1, run_preffix, run0, run1)

    plt.show()

# //////////////////////////////////////////////////////////////////////////////
