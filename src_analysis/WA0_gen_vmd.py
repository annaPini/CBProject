from parameters import *
import MDAnalysis as mda

# //////////////////////////////////////////////////////////////////////////////
def gen_cut_xtc(gro_path, xtc_path, wa_gro_path, wa_xtc_path, info):
    traj = mda.Universe(str(gro_path), str(xtc_path))
    protein = traj.select_atoms("protein")

    with mda.Writer(str(wa_gro_path), protein.n_atoms) as W:
        W.write(protein)

    with mda.Writer(str(wa_xtc_path), protein.n_atoms) as W:
        for ts in traj.trajectory[info.frames]:
            W.write(protein)


def gen_vmd(gro_path, xtc_path, vmd_path, info):
    gro_path = str(gro_path).replace('\\', '/')
    xtc_path = str(xtc_path).replace('\\', '/')

    out = """color change rgb red 0 0 1
color change rgb blue 0 0 1
color change rgb green 0 0 1
color change rgb white 0 1 1
color change rgb white 0 1 1
color change rgb cyan 0 1 0\n\n"""

    out += f"mol new {gro_path} type gro first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"
    out += f"mol addfile {xtc_path} type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"
    out += """
mol delrep 0 top
mol representation QuickSurf
mol color ResType
mol selection {protein and not resid 41 145 164}
mol material Opaque
mol addrep top
mol smoothrep top 0 4

mol representation QuickSurf
mol color Element
mol selection {resid 41 145 164}
mol material Opaque
mol addrep top
mol smoothrep top 0 4

pbc box\n"""
    out += f"puts \">>> Reference frame at index {info.ref_frame_index}\""

    with open(vmd_path, 'w') as file: file.write(out)


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    info = Info(WA_INFO_PATH)

    #######################################
    GRO = WA_POSTDIR / f"{SYSTEM_NAME}.gro"
    XTC = WA_POSTDIR / f"{TRAJECTORY_NAME}.xtc"

    WA_VMD = WA_DIR / f"{WA_NAME}.vmd"
    WA_GRO = WA_DIR / f"{WA_NAME}-prot.gro"
    WA_XTC = WA_DIR / f"{WA_NAME}-prot-wa_frames.xtc"


    gen_cut_xtc(GRO, XTC, WA_GRO, WA_XTC, info)
    gen_vmd(WA_GRO, WA_XTC, WA_VMD, info)

# //////////////////////////////////////////////////////////////////////////////
