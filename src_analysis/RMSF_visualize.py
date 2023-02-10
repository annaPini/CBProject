from parameters import *
import matplotlib.pyplot as plt

def vis_rmsf(rmsf_dir):
    rmsf = np.load(rmsf_dir)

    plt.figure()

    plt.plot(rmsf)
    plt.xlabel("atoms")
    plt.ylabel("RMSF")
    plt.show()
    print(f"Max rmsf for run {run} is {np.argmax(rmsf)}")

def vis_2rmsf(rmsf0, rmsf1, title = ''):
    print(f">>> Plotting RMSF for '{title}'...")

    frames = np.arange(rmsf0.size)
    fig, ax = plt.subplots()

    ax.set_title(title, fontdict = dict(fontsize = 20))
    ax.set_xlabel("CA atom", fontdict = dict(fontsize = 16))
    ax.set_ylabel("RMSF", fontdict = dict(fontsize = 16))
    ax.tick_params(labelsize = 12)

    if rmsf0.size != rmsf1.size:
        return

        # arrs = sorted([rmsf0, rmsf1], key=lambda arr:arr.size)
        # arrs[0] = np.append(arrs[0],arrs[0])
        # plt.plot(arrs[0])

    else:
        ax.plot(frames, rmsf0, color = BLUE, linestyle = ':', linewidth = 2)
        ax.plot(frames, rmsf1, color = HALF_RED, linestyle = '-', linewidth = 2)

        # ax.scatter(frames, rmsf0, color = BLUE, marker = '+', s = 12)
        # ax.scatter(frames, rmsf1, color = HALF_RED, marker = 'x', s = 12)


if __name__ == "__main__":
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        PATH_RMSF0 = DIR_DA_GENERAL / f"{run0}-rmsf.npy"
        PATH_RMSF1 = DIR_DA_GENERAL / f"{run1}-rmsf.npy"

        rmsf0 = np.load(PATH_RMSF0)
        rmsf1 = np.load(PATH_RMSF1)

        vis_2rmsf(rmsf0, rmsf1, title = run_preffix)

    plt.show()
