from parameters import *
import matplotlib.pyplot as plt

def vis_cmap(cmap_dir,run):
    cmap = np.load(cmap_dir)
    plt.figure(figsize=((10,10)))
    plt.title(run)
    img = plt.imshow(cmap)
    plt.colorbar(img)
    plt.show()

def compare_cmap(gro_dir, traj_dir, f1, f2, out_dir1, out_dir2):
    ### so my original idea was that we could calculate the cmap for all frames
    ### in calculate_general and then just load the cmap of specific frames when needed,
    ### but I just tried it and the output is large (4GB), also the computation of cmap is
    ### relatively light, so indeed it makes sense to calculate it each time

    traj1 = mda.Universe(str(gro_dir), str(traj_dir))[f1]
    ca_atoms1 = traj1.select_atoms("protein and name CA")
    traj2 = mda.Universe(str(gro_dir), str(traj_dir))[f2]
    ca_atoms2 = traj2.select_atoms("protein and name CA")

    print("...>>> Calculating contact map for frame {f1}...")
    d_CaCa1 = distances.distance_array(ca_atoms1.positions, ca_atoms1.positions)
    np.save(out_dir1, d_CaCa1, allow_pickle=False)

    print("...>>> Calculating contact map frame {f2}...")
    d_CaCa2 = distances.distance_array(ca_atoms2.positions, ca_atoms2.positions)
    np.save(out_dir2, d_CaCa2, allow_pickle=False)

    print("done")


if __name__ == "__main__":
    for run in RUNS:
        PATH_CMAP = DIR_DA_GENERAL / f"{run}-cmap.npy"

        #plt.title(run)
        vis_cmap(CMAP_OUT,run)
        # compare 2 cmaps
        #compare_cmap(GRO, XTC, 3, 5, CMAP_OUT1, CMAP_OUT2)
        #vis_cmap(CMAP_OUT1, run)
        #vis_cmap(CMAP_OUT2, run)
