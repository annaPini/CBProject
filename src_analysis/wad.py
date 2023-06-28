from wad_pipeline._params import *
from wad_pipeline.info import Info
from wad_pipeline.step0_rmsd import Plotter_RMSD1D_WAD
from wad_pipeline.step1_vmd import gen_cut_xtc, gen_vmd
from wad_pipeline.step2_calc import get_box_dimensions, calc_water_density
from wad_pipeline.step3_meshes import extract_colors, extract_vertices, prepare_faces, extract_faces
from wad_pipeline.step4_plot import load_vertices_faces, plot_layer, plot_density, plot_protein

from _params import *
from documentation import docs_wad

import re, argparse
import numpy as np
import MDAnalysis as mda
from argparse import RawDescriptionHelpFormatter

import matplotlib.pyplot as plt
import plotly.graph_objects as go


# //////////////////////////////////////////////////////////////////////////////
def choose_frames():
    plotter = Plotter_RMSD1D_WAD()
    plotter.vis_rmsd1d(
        rmsd_mat0 = np.load(DIR_DA_RMSD / f"{RUN_WAD}.npy"),
        title = RUN_WAD
    )
    plt.show()


# ------------------------------------------------------------------------------
def prepare_vmd():
    info = Info(PATH_WAD_INFO)

    PATH_GRO = DIR_DA_TRAJECTORIES / f"{RUN_DETAILED_ANALYSIS}.gro"
    PATH_XTC = DIR_DA_TRAJECTORIES / f"{RUN_DETAILED_ANALYSIS}.xtc"

    PATH_WAD_GRO = DIR_DA_WAD / f"{WAD_NAME}.gro"
    PATH_WAD_XTC = DIR_DA_WAD / f"{WAD_NAME}.xtc"
    PATH_WAD_VMD = DIR_DA_WAD / f"{WAD_NAME}.vmd"

    if not PATH_WAD_GRO.exists():
        gen_cut_xtc(PATH_GRO, PATH_XTC, PATH_WAD_GRO, PATH_WAD_XTC, info)

    gen_vmd(PATH_WAD_GRO, PATH_WAD_XTC, PATH_WAD_VMD, info)


# ------------------------------------------------------------------------------
def calc_density():
    info = Info(PATH_WAD_INFO)

    PATH_GRO = DIR_DA_TRAJECTORIES / f"{RUN_WAD}.gro"
    PATH_XTC = DIR_DA_TRAJECTORIES / f"{RUN_WAD}.xtc"

    traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))

    get_box_dimensions(traj, info)
    calc_water_density(traj, info)
    print()


# ------------------------------------------------------------------------------
def obj_converter():
    info = Info(PATH_WAD_INFO)

    ####################################### MTL
    print(f">>> Parsing {WAD_NAME}.mtl")
    mat_names, mat_colors, colors = extract_colors(DIR_DA_WAD / f"{WAD_NAME}.mtl")


    ####################################### OBJ
    print(f">>> Parsing {WAD_NAME}.obj")
    with open(DIR_DA_WAD / f"{WAD_NAME}.obj", 'r') as file: obj = file.read()


    repr_indexes = [g.start() for g in re.finditer("\ng.*", obj)]
    n_repr = str(len(repr_indexes))
    print(f"...>>> Identified {n_repr} representation(s)")
    repr_indexes.append(len(obj))


    box_str = obj[:repr_indexes[0]]
    print(f"...>>> Extracting box vertices...")
    box = extract_vertices(box_str)

    xmin, xmax = np.unique(box[:,0])
    ymin, ymax = np.unique(box[:,1])
    zmin, zmax = np.unique(box[:,2])

    x_scale = info.xbox / (xmax - xmin)
    y_scale = info.ybox / (ymax - ymin)
    z_scale = info.zbox / (zmax - zmin)

    info.update(
        n_repr = int(n_repr),
        x_scale = x_scale,
        y_scale = y_scale,
        z_scale = z_scale,
        x_offset = xmin * x_scale,
        y_offset = ymin * y_scale,
        z_offset = zmin * z_scale,
    )


    for i,repr_start in enumerate(repr_indexes[:-1]):
        repr_end = repr_indexes[i + 1]
        repr = obj[repr_start : repr_end]

        print(f"...>>> Extracting vertices for repr {i}...")
        n_vertices = extract_vertices(repr, f"{WAD_NAME}-mesh_v{i}.npy")

        print(f"...>>> Extracting faces for repr {i}...")
        faces_str = prepare_faces(repr)

        if faces_str:
            extract_faces(faces_str, n_vertices, colors, f"{WAD_NAME}-mesh_f{i}.csv")
        else:
            print(f"...>>> No faces found for repr {i}, skipping...")


# ------------------------------------------------------------------------------
def visualize():
    info = Info(PATH_WAD_INFO)
    fig = go.Figure()

    plot_density(fig, info)
    plot_protein(fig, info)

    fig.update_layout(margin = dict(l = 0, r = 0, b = 0, t = 0))
    fig.update_scenes(xaxis_visible = False, yaxis_visible = False, zaxis_visible = False, bgcolor = "rgb(0,0,50)")
    fig.show()

    print("""
    +========= COLOR SCHEMA =========+
    | GREEN:    hydrophilic residues |
    | YELLOW:   hydrophobic residues |
    | RED:      active site          |
    +================================+
    """)


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = docs_wad["main"], formatter_class = RawDescriptionHelpFormatter)
    parser.add_argument("-0", "--rmsd", action = "store_true", help = docs_wad["rmsd"])
    parser.add_argument("-1", "--vmd", action = "store_true", help = docs_wad["vmd"])
    parser.add_argument("-2", "--calc", action = "store_true", help = docs_wad["calc"])
    parser.add_argument("-3", "--meshes", action = "store_true", help = docs_wad["meshes"])
    parser.add_argument("-4", "--plot", action = "store_true", help = docs_wad["plot"])
    args = parser.parse_args()

    if not (args.rmsd | args.vmd | args.calc | args.meshes | args.plot):
        print(">>> Refer to the documentation (wad.py -h) for how to follow correctly the WAD pipeline. Closing...")

    if args.rmsd: choose_frames()
    if args.vmd: prepare_vmd()
    if args.calc: calc_density()
    if args.meshes: obj_converter()
    if args.plot: visualize()


# //////////////////////////////////////////////////////////////////////////////
