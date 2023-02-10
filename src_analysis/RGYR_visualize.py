from parameters import *
import matplotlib.pyplot as plt

def vis_2rgyr(rgyr0, rgyr1, title = ''):
    print(f">>> Plotting RGYR for '{title}'...")

    frames = np.arange(rgyr0.size)
    fig, ax = plt.subplots()

    ax.set_title(title, fontdict = dict(fontsize = 20))
    ax.set_xlabel("Frame", fontdict = dict(fontsize = 16))
    ax.set_ylabel("Radius of Gyration", fontdict = dict(fontsize = 16))
    ax.tick_params(labelsize = 12)

    ax.scatter(frames, rgyr0, color = BLUE, marker = '+', s = 12)
    ax.scatter(frames, rgyr1, color = HALF_RED, marker = 'x', s = 12)


##### RGYR arrays are always the same size (n = number of frames), so this function from RMSF (x = number of atoms) doesn't apply here
# def vis_2rgyr(rgyr_dir1,rgyr_dir2):
#     rgyr1 = np.load(rgyr_dir1)
#     rgyr2 = np.load(rgyr_dir2)
#
#     if rgyr1.size != rgyr2.size:
#         arrs = sorted([rgyr1, rgyr2], key=lambda arr:arr.size)
#         arrs[0] = np.append(arrs[0],arrs[0])
#
#         plt.plot(arrs[0])
#         plt.plot(arrs[1])
#
#     else:
#
#         plt.plot(rgyr1)
#         plt.plot(rgyr2)
#
#     plt.show()

if __name__ == "__main__":
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        PATH_RGYR0 = DIR_DA_GENERAL / f"{run0}-rgyr.npy"
        PATH_RGYR1 = DIR_DA_GENERAL / f"{run1}-rgyr.npy"

        rgyr0 = np.load(PATH_RGYR0)
        rgyr1 = np.load(PATH_RGYR1)

        vis_2rgyr(rgyr0, rgyr1, title = run_preffix)

    plt.show()
