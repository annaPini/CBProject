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

def calc_water_density(traj, info, frames):
    xstep = info.xbox / WA_DIVISIONS
    ystep = info.ybox / WA_DIVISIONS
    zstep = info.zbox / WA_DIVISIONS

    density = np.zeros((WA_DIVISIONS, WA_DIVISIONS, WA_DIVISIONS))

    for i,ts in enumerate(traj.trajectory[frames]):
        print(f"{i+1}/{len(frames)}", end = ' ')

        for x,y,z in traj.select_atoms("name OH2").positions:
            if not (0 < x < info.xbox): continue
            if not (0 < y < info.ybox): continue
            if not (0 < z < info.zbox): continue

            i = int(x // xstep)
            j = int(y // ystep)
            k = int(z // zstep)

            density[i,j,k] += 1

    density /= len(frames)

    # ################################################################################ Output
    np.save(WA_DIR / (f"{WA_NAME}-wet_density.npy"), density, allow_pickle = False)


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    info = Info(WA_INFO_PATH)

    #######################################
    GRO = WA_POSTDIR / f"{SYSTEM_NAME}.gro"
    XTC = WA_POSTDIR / f"{TRAJECTORY_NAME}.xtc"

    traj = mda.Universe(str(GRO), str(XTC))

    # frames = list(range(20))
    # frames = list(range(50))
    # frames = info.frames[:20]
    # frames = info.frames[:50]
    frames = info.frames[:20]

    get_box_dimensions(traj, info)
    calc_water_density(traj, info, frames)
    print()

# //////////////////////////////////////////////////////////////////////////////
