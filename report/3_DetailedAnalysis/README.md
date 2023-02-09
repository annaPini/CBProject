# Detailed Analysis
Analysis of **mt2_rep0** trajectory (ie protein that "opens up" at frame ~5500)

- Choosing the trajectory of interest and explaining why
    - Clustering to better justify the conformation families as observed in the RMSD map
    - PCA to compare with clustering

- Picking the frames of interest with RMSD
    - GIF of protein opening up (vs wild type analogous) [observation 1]
    - Compare cmap when protein is open vs when it's closed

- Observation of active site getting closed [observation 2]
    - VMD gifs of 3 different visualizations of the same region near the active site (closes up because of the mutation) (comparison with wt)
    - Example of using WAD to compare frames when it's open vs when it's closed (frame 0 vs frame 5000)
        - WAD vs analogous methods

- Pyinteraph to try to explain either of these observations [1,2]



## Pyinteraph
We have observed that during the trajectory the lower part of the protein enhances and at the same time the upper part moves too. We wondered if there are some correlations between these fluctuations, we wanted to explore more in detail what kind of interaction occur between the residues involved.
We adressed this problem by means of Pyinteraph.

## WAD
