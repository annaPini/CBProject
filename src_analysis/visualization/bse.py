from .calculations import get_cov_matrix
from matplotlib import pyplot as plt
import seaborn as sns

# //////////////////////////////////////////////////////////////////////////////
def vis_BSE(bse_naive, bse_alte):
    plt.plot(bse_rmsd_all)
    plt.plot(bse_rmsd_all_B)

    cov = get_cov_matrix(size = rmsd.size, length = 50)
    sns.heatmap(cov)

    plt.show()

# //////////////////////////////////////////////////////////////////////////////
