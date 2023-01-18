from parameters import *
import pandas as pd
import re

# //////////////////////////////////////////////////////////////////////////////
####################################### MATERIALS
def extract_colors(mtl_path):
    with open(mtl_path, 'r') as file: mtl = file.read()
    print("...>>> Extracting colors...")

    mat_names_raw = re.findall("\nnewmtl.*", mtl)
    mat_colors_raw = re.findall("\nKd.*", mtl)

    mat_names = map(lambda s: s[8:], mat_names_raw)
    mat_colors = map(lambda s: f"rgb({ s[4:].replace(' ', ',') })", mat_colors_raw)
    colors = dict(zip(mat_names, mat_colors))

    return mat_names, mat_colors, colors

####################################### VERTICES
def extract_vertices(repr, out_name = ''):
    vertices_raw = re.findall("\nv\s.*", repr)
    vertices = np.array(list(map(lambda s: s[3:].split(' '), vertices_raw)), dtype = float)

    if out_name:
        print(f"...>>> Saving vertices to '{out_name}'...")
        np.save(WA_DIR / out_name, vertices, allow_pickle = False)
        return vertices.shape[0]

    return vertices

####################################### FACES
def prepare_faces(repr):
    faces_start = re.search("\nusemtl.*", repr)
    if faces_start is None: return
    return repr[faces_start.start():]

def extract_faces(faces_str, n_vertices, colors, out_name):
    re_materials = re.finditer("usemtl.*", faces_str)
    materials = map(lambda m: (m.start(), m.group()), re_materials)


    matches = np.array(list(re.finditer("[\s\-]\d+[\/]", faces_str)))
    matches = np.reshape(matches, (matches.size // 3, 3))

    faces = pd.DataFrame(
        columns = ['i', 'j', 'k'],
        data = np.vectorize( lambda m: n_vertices + int( m.group()[:-1] ) )(matches)
    )

    raw_start = np.vectorize( lambda m: m.start() )(matches[:,0])
    color = np.full(raw_start.size, "rgb(_.--,_.--,_.--)")

    for id,m in materials:
        color[raw_start > id] = colors[m[7:]]
    faces["color"] = color

    print(f"...>>> Saving faces to '{out_name}'...")
    faces.to_csv(WA_DIR / out_name, index = False)


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    info = Info(WA_INFO_PATH)

    ####################################### MTL
    print(f">>> Parsing {WA_NAME}.mtl")
    mat_names, mat_colors, colors = extract_colors(WA_DIR / f"{WA_NAME}.mtl")


    ####################################### OBJ
    print(f">>> Parsing {WA_NAME}.obj")
    with open(WA_DIR / f"{WA_NAME}.obj", 'r') as file: obj = file.read()


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
        n_vertices = extract_vertices(repr, f"{WA_NAME}-mesh_v{i}.npy")

        print(f"...>>> Extracting faces for repr {i}...")
        faces_str = prepare_faces(repr)

        if faces_str:
            extract_faces(faces_str, n_vertices, colors, f"{WA_NAME}-mesh_f{i}.csv")
        else:
            print(f"...>>> No faces found for repr {i}, skipping...")

# //////////////////////////////////////////////////////////////////////////////
