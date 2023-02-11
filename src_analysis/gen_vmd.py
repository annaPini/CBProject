from parameters import *

# //////////////////////////////////////////////////////////////////////////////
def vmd_start(gro_path, xtc_path, f0, f1):
    gro_path = str(gro_path).replace('\\', '/')
    xtc_path = str(xtc_path).replace('\\', '/')

    return f"""
mol new {gro_path} type gro first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile {xtc_path} type xtc first {f0} last {f1} step 1 filebonds 1 autobonds 1 waitfor all\n
mol delrep 0 top\n"""[1:]

def vmd_end(pbc):
    return """\n
lappend viewplist [molinfo top]
set topmol [molinfo top]
# done with molecule 0
foreach v $viewplist {
  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)
}
set fixedlist {}
foreach v $fixedlist {
  molinfo $v set fixed 1
}
unset viewplist
unset fixedlist
mol top $topmol
unset topmol""" + ("\n\npbc box" if pbc else '')

def vmd_writer(vmd_middle):
    def vmd_generation(vmd_name, gro_path, xtc_path, f0 = 0, f1 = 2000, pbc = False, *args, **kwargs):
        out = vmd_start(gro_path, xtc_path, f0, f1) + vmd_middle(*args, **kwargs) + vmd_end(pbc)
        vmd_path = DIR_DA_VMD / f"{vmd_name}.vmd"
        with open(vmd_path, 'w') as file: file.write(out)
    return vmd_generation

def mol_commands(*commands):
    return '\n'.join(("mol " + command for command in commands))



# ------------------------------------------------------------------------------
@vmd_writer
def gen_vmd_basic():
    return "\n\n".join((
        mol_commands("representation QuickSurf", "color Structure", "selection {protein}", "material Opaque", "addrep top", "smoothrep top 0 4"),
        mol_commands("representation VDW 0.3", "color Name", "selection {not protein}", "material Glass1", "addrep top"),
        "set viewpoints([molinfo top]) {{{1 0 0 -50.9224} {0 1 0 -51.1731} {0 0 1 -50.9419} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} {{0.0144732 0 0 0} {0 0.0144732 0 0} {0 0 0.0144732 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}",
    ))


# @vmd_writer
# def gen_vmd_focus_active_site():
#     return "\n\n".join((
#         mol_commands("representation NewCartoon 0.300000 10.000000 4.100000 0", "color ColorID 6", "selection {protein}", "material BlownGlass", "addrep top", "smoothrep top 0 4"),
#         mol_commands("representation NewCartoon 0.300000 10.000000 4.100000 0", "color ColorID 1", "selection {resid 163 164 165}", "material Opaque", "addrep top", "smoothrep top 1 4"),
#         mol_commands("representation NewCartoon 0.300000 10.000000 4.100000 0", "color ColorID 7", "selection {resid 40 41 42 144 145 146}", "material Opaque", "addrep top", "smoothrep top 2 4"),
#         mol_commands("representation Licorice 0.300000 12.000000 12.000000", "color Name", "selection {resid 41 145 164}", "material Opaque", "addrep top", "smoothrep top 3 4"),
#         mol_commands("representation VDW 1.000000 12.000000", "color Name", "selection {resid 41 145 164}", "material Glass3", "addrep top", "smoothrep top 4 4"),
#         "set viewpoints([molinfo top]) {{{1 0 0 -72.9} {0 1 0 -54.06} {0 0 1 -50.94} {0 0 0 1}} {{0.380485 -0.721862 0.577897 0} {-0.74395 -0.610101 -0.272297 0} {0.549197 -0.326324 -0.769205 0} {0 0 0 1}} {{0.156988 0 0 0} {0 0.156988 0 0} {0 0 0.156988 0} {0 0 0 1}} {{1 0 0 0.152241} {0 1 0 -0.82334} {0 0 1 -0.311366} {0 0 0 1}}}",
#     ))


@vmd_writer
def gen_vmd_focus_active_site():
    return "\n\n".join((
        mol_commands("representation Ribbons 0.300000 12.000000 2.000000", "color ColorID 6", "selection {protein}", "material BlownGlass", "addrep top", "smoothrep top 0 4"),
        mol_commands("representation Ribbons 0.300000 12.000000 2.000000", "color ColorID 1", "selection {resid 163 164 165}", "material Opaque", "addrep top", "smoothrep top 1 4"),
        mol_commands("representation Ribbons 0.300000 12.000000 2.000000", "color ColorID 7", "selection {resid 40 41 42 144 145 146}", "material Opaque", "addrep top", "smoothrep top 2 4"),
        mol_commands("representation Licorice 0.300000 12.000000 12.000000", "color ColorID 2", "selection {resid 40 187}", "material Opaque", "addrep top", "smoothrep top 4 4"),
        mol_commands("representation Licorice 0.300000 12.000000 12.000000", "color Name", "selection {resid 41 145 164}", "material Opaque", "addrep top", "smoothrep top 3 4"),
        mol_commands("representation VDW 1.000000 12.000000", "color Name", "selection {resid 41 145 164}", "material BlownGlass", "addrep top", "smoothrep top 5 4"),
        "set viewpoints([molinfo top]) {{{1 0 0 -72.9} {0 1 0 -54.06} {0 0 1 -50.94} {0 0 0 1}} {{-0.864712 -0.499059 0.0551708 0} {-0.397742 0.61376 -0.681844 0} {0.306443 -0.611571 -0.729263 0} {0 0 0 1}} {{0.188235 0 0 0} {0 0.188235 0 0} {0 0 0.188235 0} {0 0 0 1}} {{1 0 0 0.152241} {0 1 0 -0.82334} {0 0 1 -0.311366} {0 0 0 1}}}",
    ))


@vmd_writer
def gen_vmd_focus_opening():
    return "\n\n".join((
        mol_commands("representation NewCartoon 0.300000 10.000000 4.100000 0", "color ColorID 11", "selection {residue <= 306}", "material Opaque", "addrep top", "smoothrep top 0 4"),
        mol_commands("representation NewCartoon 0.300000 10.000000 4.100000 0", "color ColorID 12", "selection {residue > 306}", "material Opaque", "addrep top", "smoothrep top 1 4"),
        "set viewpoints([molinfo top]) {{{1 0 0 -72.9} {0 1 0 -54.06} {0 0 1 -50.94} {0 0 0 1}} {{0.080682 -0.0753226 -0.993764 0} {0.406577 0.912779 -0.0361878 0} {0.909954 -0.401147 0.104286 0} {0 0 0 1}} {{0.0434286 0 0 0} {0 0.0434286 0 0} {0 0 0.0434286 0} {0 0 0 1}} {{1 0 0 -0.397759} {0 1 0 0.29666} {0 0 1 -0.311366} {0 0 0 1}}}",
    ))




# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":

    ############### Progress of the MD post-processing workflow
    GRO_MT2R0   = DIR_DA_TRAJECTORIES / "mt2_rep0.gro"
    XTC_MT2R0_A = DIR_DEXTRA / "md_plain.xtc"
    XTC_MT2R0_B = DIR_DEXTRA / "md-nojump.xtc"
    XTC_MT2R0_C = DIR_DEXTRA / "md-center.xtc"
    XTC_MT2R0   = DIR_DA_TRAJECTORIES / "mt2_rep0.xtc"

    GRO_MT2R1 = DIR_DA_TRAJECTORIES / "mt2_rep1.gro"
    XTC_MT2R1 = DIR_DA_TRAJECTORIES / "mt2_rep1.xtc"

    GRO_WT2R0 = DIR_DA_TRAJECTORIES / "wt2_rep0.gro"
    XTC_WT2R0 = DIR_DA_TRAJECTORIES / "wt2_rep0.xtc"

    GRO_WT2R1 = DIR_DA_TRAJECTORIES / "wt2_rep1.gro"
    XTC_WT2R1 = DIR_DA_TRAJECTORIES / "wt2_rep1.xtc"

    # gen_vmd_basic("postproc_0", GRO_MT2R0, XTC_MT2R0_A, pbc = True)
    # gen_vmd_basic("postproc_1", GRO_MT2R0, XTC_MT2R0_B, pbc = True)
    # gen_vmd_basic("postproc_2", GRO_MT2R0, XTC_MT2R0_C, pbc = True)
    # gen_vmd_basic("postproc_3", GRO_MT2R0, XTC_MT2R0  , pbc = True)

    # gen_vmd_focus_active_site("focus_AS_mtr0", GRO_MT2R0, XTC_MT2R0, f1 = -1)
    # gen_vmd_focus_active_site("focus_AS_mtr1", GRO_MT2R1, XTC_MT2R1, f1 = -1)
    # gen_vmd_focus_active_site("focus_AS_wtr0", GRO_WT2R0, XTC_WT2R0, f1 = -1)
    # gen_vmd_focus_active_site("focus_AS_wtr1", GRO_WT2R1, XTC_WT2R1, f1 = -1)

    gen_vmd_focus_opening("focus_open_mtr0", GRO_MT2R0, XTC_MT2R0, f0 = 5000, f1 = 6000)
    gen_vmd_focus_opening("focus_open_mtr1", GRO_MT2R1, XTC_MT2R1, f0 = 4000, f1 = 5000)
    gen_vmd_focus_opening("focus_open_wtr0", GRO_WT2R0, XTC_WT2R0, f0 = 5000, f1 = 6000)
    gen_vmd_focus_opening("focus_open_wtr1", GRO_WT2R1, XTC_WT2R1, f0 = 5000, f1 = 6000)

    #### residues with high RMSF on mutated, compared to wt:
    # (residue >= 45 and residue <= 50) or (residue >= 289 and residue <= 312) or (residue >= 350 and residue <= 356) or (residue >= 593)
