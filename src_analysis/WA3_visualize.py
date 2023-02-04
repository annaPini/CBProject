from parameters import *
import pandas as pd
import plotly.graph_objects as go

# //////////////////////////////////////////////////////////////////////////////
def load_vertices_faces(info):
    vertices = []; faces = []
    for i in range(info.n_repr):
        vertices_dir = WA_DIR / f"{WA_NAME}-mesh_v{i}.npy"
        faces_dir =    WA_DIR / f"{WA_NAME}-mesh_f{i}.csv"

        vertices.append(np.load(vertices_dir))
        faces.append(pd.read_csv(faces_dir))

        vertices[i][:,0] *= info.x_scale
        vertices[i][:,1] *= info.y_scale
        vertices[i][:,2] *= info.z_scale

        vertices[i][:,0] -= info.x_offset
        vertices[i][:,1] -= info.y_offset
        vertices[i][:,2] -= info.z_offset

    return vertices, faces

def plot_layer(vertices, faces = None):
    kwargs = dict(
        x = vertices[:,0],
        y = vertices[:,1],
        z = vertices[:,2],
    )

    if faces is not None:
        kwargs['i'] = faces.i
        kwargs['j'] = faces.j
        kwargs['k'] = faces.k
        kwargs["facecolor"] = faces.color

    fig.add_trace(
        go.Mesh3d(**kwargs)
    )

def plot_density(info):
    density = np.load(WA_DIR / f"{WA_NAME}-wet_density.npy")

    X, Y, Z = np.mgrid[
        0 : info.xbox : complex(0, WA_DIVISIONS),
        0 : info.ybox : complex(0, WA_DIVISIONS),
        0 : info.zbox : complex(0, WA_DIVISIONS),
    ]

    fig.add_trace(
        go.Volume(
            x = X.flatten(),
            y = Y.flatten(),
            z = Z.flatten(),
            value = -density.flatten(),
            isomin = -DENSITY_TRESHOLD_UPPER,
            isomax = -DENSITY_TRESHOLD_LOWER,
            opacity = 0.05, # needs to be small to see through all surfaces
            opacityscale = "uniform",
            surface_count = 50, # needs to be a large number for good volume rendering
            colorscale = ["blue", "white", "magenta"],
        )
    )

def plot_protein(info):
    for v,f in zip(*load_vertices_faces(info)):
        plot_layer(v, f)

# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    info = Info(WA_INFO_PATH)
    fig = go.Figure()

    plot_density(info)
    plot_protein(info)

    fig.update_layout(margin = dict(l = 0, r = 0, b = 0, t = 0))
    fig.update_scenes(xaxis_visible = False, yaxis_visible = False, zaxis_visible = False, bgcolor = "rgb(0,0,50)")
    fig.show()

# //////////////////////////////////////////////////////////////////////////////

# GREEN:    hydrophilic residues
# YELLOW:   hydrophobic residues
# RED:      active site
