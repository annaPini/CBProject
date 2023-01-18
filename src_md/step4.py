from header import *
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj

for run in RUNS:
    POSTDIR = get_POSTDIR(run)
    GRO = POSTDIR / f"{SYSTEM_NAME}.gro"
    XTC = POSTDIR / f"{SYSTEM_NAME}-rottrans.xtc"
    OUT = POSTDIR / f"{SYSTEM_NAME}-MDA.xtc"

    print(f">>> Reading trajectory for '{run}'...")
    traj = mda.Universe(str(GRO), str(XTC))

    print("...>>> Centering trajectory...")
    AlignTraj(traj, traj, select = "name CA", filename = str(OUT)).run()

    print("...>>> Done.")

    wa_dir = SIMDIR / run / "wetness_analysis"
    if not WA_DIR.exists(): WA_DIR.mkdir()
