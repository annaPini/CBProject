from parameters import *
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import multivariate_normal

#Naive version of BSE: the mean value is not stable since we are leaving out the last block.
def BSE_naive(arr):
    v = arr
    n = len(arr)
    l = []
    #for i in [2**j for j in range(1,10)]:
    for i in range(1,n//3):
        d = n//i
        w = v[:d*i].reshape(d,i)
        m = w.mean(axis=1)
        l.append(m.var()/(m.shape[0]-1))
    l = np.asarray(l)
    l = np.sqrt(l)
    return l

# Returns a vector s(m) with the std-dev computed according to the Block Analysis method. The average is stable, needs more testing.
def BSE_alternate(arr):
    v = arr
    n = len(arr)
    a = []
    s = []
    for i in range(1,n//3):
        d = n//i
        r = n-d*i
        w = v[:d*i].reshape(d,i)
        m = w.mean(axis=1)
        avg = m.mean()
        weights = np.full(d,i)
        if d*i != n:
            m_d = v[d*i:].mean()
            avg = (d*i*avg+r*m_d)/(d*i+r)
            weights=np.append(weights,r)
            m=np.append(m,m_d)
        m2 = (m - avg)**2
        weights = weights/weights.sum()
        m2 *= weights
        var = m2.sum()/(m2.shape[0]-1)
        s.append(var)
        a.append(avg)
    s = np.asarray(s)
    s = np.sqrt(s)
    a = np.asarray(a)
     #return a,s # a is the average, check that it is constant.
    return s

def calc_rmsd1D (rmsd_dir):
    rmsd = np.load(rmsd_dir)
    return rmsd[0] # riga 0, tutte le colonne
    #traj = mda.Universe(str(gro_dir),str(traj_dir))
    #rmsd = RMSD(traj,select='name CA', groupselections=['backbone', 'resid 20:50'])
    # rmsd.run()


#Create a covariance matrix with a given correlation length
def get_cov_matrix(size,length):
    x = np.arange(size)
    cov = np.exp(-(1/length)*(x-np.atleast_2d(x).T)**2)
    return cov

# Expected error on mean obtained from the covariance matrix of the random process generator
def get_err_on_avg(cov_mat):
    cov = cov_mat
    return np.sqrt(cov.mean())



if __name__ == "__main__":
    #for run in RUNS:

        #PATH_GRO     = DIR_DA_TRAJECTORIES / f"{run}.gro"
        #PATH_XTC     = DIR_DA_TRAJECTORIES / f"{run}.xtc"
        #RMSD_OUT  = DIR_DA_BLOCK_ANALYSIS / f"{run}-rmsd.npy"

    PATH_GRO_MT2_0     = DIR_DA_TRAJECTORIES / f"mt2_rep0.gro"
    PATH_XTC_MT2_0     = DIR_DA_TRAJECTORIES / f"mt2_rep0.xtc"
    RMSD_MT2_0  = DIR_DA_GENERAL / f"mt2_rep0-rmsd.npy"

    rmsd = calc_rmsd1D (RMSD_MT2_0)

    bse_rmsd_all = BSE_naive(rmsd)
    bse_rmsd_all_B = BSE_alternate(rmsd)

    # plt.plot(bse_rmsd_all)
    # plt.plot(bse_rmsd_all_B)

    cov = get_cov_matrix(size = rmsd.size,length=50)
    sns.heatmap(cov)

    plt.show()
