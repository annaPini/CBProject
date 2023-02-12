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
One of the most intriguing part of the system under study was the water surrounding the protein. In the first part of this section we present a brief introduction on water models and in particular on the $TIP3P$ model, that was employed in our simulation. The second part is devoted to a simulation of the water distribution around the protein. We observe how water molecules arrange according to the hydrophobic or hydrophilic nature of the residues. Finally we estimate the accessible surface with the GROMACS tool SASA. 

Modelling water is a challenging task since its microscopic structure at the liquid state it is not known yet and many different models are available. The choice of the water model is determined by
- the size of the system
- the simulation period
- the physical properties that have to be simulated with the most accurancy. 

The most suitable model for our system was the $TIP3P$, This water model is widely used because it is not too computationally demanding. Even to simulate a simple protein in water, a large number of water molecule is required, and their positions and velocity have to be updated at every time step. It is then necessary to have a quite simple water model. More complex models are available, but it is possible to employ them only on small systems, with less than $1000$ water molecules, and for short simulation periods $(\simeq1ps)$. These are not the lenght and time scales of biophysical systems. 

This model, as many other frequently used, ignores the proton hopping from one molecule to another, which favours the cluster formation. It describes the water as a planar molecule, with three interaction sites, an Oxygen atom and two Hidrogens. The atoms interact due to the Coulomb interaction and are modelled as point like particles. The long range electrostatic interaction is truncated after a cutoff distance. The Oxygen atom also interacts with other molecules by means of the Lennard-Jones interaction.


  
