from parameters import *
import MDAnalysis as mda

# //////////////////////////////////////////////////////////////////////////////
def get_box_dimensions(traj, info):
    traj.trajectory[info.ref_frame]

    xbox, ybox, zbox, _, _, _ = traj.dimensions
    info.update(
        xbox = float(xbox),
        ybox = float(ybox),
        zbox = float(zbox),
    )

def calc_water_density(traj, info):
    # initialization of values
    xstep = info.xbox / WA_DIVISIONS
    ystep = info.ybox / WA_DIVISIONS
    zstep = info.zbox / WA_DIVISIONS
    steps = np.array([xstep, ystep, zstep])
    frames_inverse = 1 / len(info.frames)

    # initialize the density matrix with very high values
    high_concentration = (DENSITY_TRESHOLD_UPPER + 1) * len(info.frames)
    density = np.zeros((WA_DIVISIONS, WA_DIVISIONS, WA_DIVISIONS)) + high_concentration

    # finding the smallest sub-cube that spans the entire protein
    traj.trajectory[info.ref_frame]
    prot_x, prot_y, prot_z = traj.select_atoms("protein").positions.T

    min_prot_x = prot_x.min()
    min_prot_y = prot_y.min()
    min_prot_z = prot_z.min()

    max_prot_x = prot_x.max()
    max_prot_y = prot_y.max()
    max_prot_z = prot_z.max()

    # set the inside of the sub-cube to 0 in the density matrix
    density[
        int(min_prot_x // xstep) : int(max_prot_x // xstep) + 1,
        int(min_prot_y // ystep) : int(max_prot_y // ystep) + 1,
        int(min_prot_z // zstep) : int(max_prot_z // zstep) + 1,
        ] = 0


    print(">>> Calculating average water density, progress:", end = ' ')
    for i,ts in enumerate(traj.trajectory[info.frames]):
        if not i % 25: print(f"{i+1}/{len(info.frames)}", end = ' ')

        water = traj.select_atoms("name OH2").positions
        water_in_box = water[
            (water[:,0] > 0) & (water[:,0] < info.xbox) &
            (water[:,1] > 0) & (water[:,1] < info.ybox) &
            (water[:,2] > 0) & (water[:,2] < info.zbox)
        ]

        for i,j,k in (water_in_box // steps).astype(int):
            density[i,j,k] += frames_inverse

    ############################################################################
    np.save(WA_DIR / (f"{WA_NAME}-wet_density.npy"), density, allow_pickle = False)


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    info = Info(WA_INFO_PATH)

    #######################################
    GRO = WA_POSTDIR / f"{SYSTEM_NAME}.gro"
    XTC = WA_POSTDIR / f"{TRAJECTORY_NAME}.xtc"

    traj = mda.Universe(str(GRO), str(XTC))

    get_box_dimensions(traj, info)
    calc_water_density(traj, info)
    print()

# //////////////////////////////////////////////////////////////////////////////
