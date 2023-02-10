from parameters import *
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.pca import PCA, cosine_content

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import plotly.express as px


################################################################################

PATH_GRO = DIR_DA_TRAJECTORIES / f"{CURRENT_RUN}.gro"
PATH_XTC = DIR_DA_TRAJECTORIES / f"{CURRENT_RUN}.xtc"

PATH_PCA_CUMVAR = DIR_DA_SPECIFIC / f"{CURRENT_RUN}-pca_cumvar.npy"
PATH_PCA_PCOMPS = DIR_DA_SPECIFIC / f"{CURRENT_RUN}-pca_pcomponents.npy"
PATH_PCA_SPACE  = DIR_DA_SPECIFIC / f"{CURRENT_RUN}-pca_space.npy"


traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))



################################################################################

if PATH_PCA_CUMVAR.exists():
    cumvar = np.load(PATH_PCA_CUMVAR)
    space  = np.load(PATH_PCA_SPACE)
    pcomps = np.load(PATH_PCA_PCOMPS)

else:
    print(">>> Running PCA...")
    rec_pca = PCA(traj, select = "name CA")
    rec_pca.run()
    print(">>> Done.")

    cumvar = rec_pca.cumulated_variance
    pcomps = rec_pca.p_components
    space = rec_pca.transform(traj.select_atoms("name CA"), 3)

    np.save(PATH_PCA_CUMVAR, cumvar, allow_pickle = False)
    np.save(PATH_PCA_PCOMPS, pcomps, allow_pickle = False)
    np.save(PATH_PCA_SPACE, space, allow_pickle = False)

################################################################################
# fg = px.line(x = np.arange(cumvar.shape[0]),
#         y = cumvar,
#        labels = {"x" : "components", "y" : "cumulated variance"},
#        range_x = [0, 100])
#
# n_pcs = np.where(cumvar > 0.95)[0][0]
# print("NPCS:", n_pcs)
# fg.show()



fig = plt.figure(figsize=(6,6))
ax = Axes3D(fig)


sc = ax.scatter(space[:,0], space[:,1], space[:,2],
                c = np.arange(space.shape[0]), s = 40, marker = "o", alpha=1)
# legend
plt.legend(*sc.legend_elements(), bbox_to_anchor=(1.05, 1), loc=2)



pca_data = pd.DataFrame(space,columns=["first_comp","second_comp","third_comp"])

# df = pca_data
last_frames = 500
df = pca_data.tail(last_frames)


fig = px.scatter_3d(df, x="first_comp", y="second_comp", z="third_comp",
              color=df.index.values, width=900,height=800)
fig.show()


fg0 = px.line(x=np.arange(pcomps[0].shape[0]),y=pcomps[0])
fg0.show()

fg1 = px.line(x=np.arange(pcomps[0].shape[0]),y=pcomps[1])
fg1.show()

fg2 = px.line(x=np.arange(pcomps[0].shape[0]),y=pcomps[2])
fg2.show()
