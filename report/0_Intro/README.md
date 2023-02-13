# Introduction
This documentation serves as a report for the project of Computational Biology course 2022-2023 (University of Trento). The repository contains the source code and initial input files necessary to reproduce all the steps detailed by the report, as well as instructions of use inside the documentation.

The focus of this project is to analyse the dynamic of a main protease protein. The software we used to simulate the trajectory is GROMACS, which integrates the equations of motion by means of the Verlet algorithm. Once the trajectories were obtained we analysed them, focusing on general aspects of the trajectories, as well as frames of interest, the active sites and the distribution of water around the protein.


## Physical system
The Sars-Cov-2 virus is organised in two polyproteins which contain shorter proteins. We focus on one of the two main protease proteins, $M^{pro}$. They play a fundamental role in the replication process of the virus, by cutting the polyproteins in shorter chains. The M^pro^ protein is one of the most important drug targets, due to the fact that it is absent in human body and that it is directly involved in the virus replication process ([source](https://www.biorxiv.org/content/10.1101/2022.06.28.497978v2)).

From the crystallographic structure the $M^{pro}$ results to be divided in two subunits (monomers). Each subunit consists of the same sequence of 306 aminoacids. Each one has an active site, known as the catalytic dyad, which consists in two residues: **H41-C145**. Structurally, it comprehends three domains: ```Domain I: residues 8-101```, ```Domain II: residues 102-184```, ```Domain III: residues 201-306```. There is also a long loop, residues 185-200, which connects domains II and III.


## Molecular dynamics
As mentioned before, all the trajectories are obtained by means of Molecular dynamics simulations. It is then worth to briefly explain this method. In the system under study there are not dissipative forces, so Hamiltonian formalism can be employed. The equations of motion are obtained by taking the time derivative of the generalized coordinates and momenta
$$\dot{q}_{\alpha} = \frac{\partial \mathcal{H}}{\partial p_{\alpha}}  \dot{p}_{\alpha} =- \frac{\partial \mathcal{H}}{\partial q_{\alpha}} .$$
These are first order differential equations and the algorithms we are going to introduce are meant to numerically solve them. The positions and the velocity of each atom have to be found as a function of time and by keeping the total energy constant. This means the total derivative of the Hamiltonian is
$$\frac{d \mathcal{H}}{dt} = 0 .$$
Our system is also ergodic, meaning that during its evolution it visits almost all the accessible microstates. If we apply the algorithm to integrate the equations of motions many times, we will always find different results, but, due to the ergodicity, the error does not diverge. Ergodicity also allows to exactly evaluate the macroscopic parameters by means of an average on time. For our purposes, a discrete average is necessary
$$ A = \braket{a} = \frac{1}{M}\sum_{n=1}^{M} a(x_{n \Delta t}),$$
where $\Delta t$ is the time step and $ a(x_{n \Delta t})$ is the macroscopic parameter evaluated at each time step on each element of the ensemble which is visited by the system.

### Basic components of a Molecular dynamics simulation
- The Model, that is the force field
- Calculation of forces and energies
- The algorithm to integrate the equations of motion

### Verlet algorithm (1967)
The Taylor expansion to the second order of the position of the atom at $(t+ \Delta t)$ is
$$\vec{r}_{i} (t+ \Delta t) = \vec{r}_{i}(t) + \vec{\dot{r}}_{i} \Delta t + \frac{1}{2}\vec{\ddot{r}}_{i}(\Delta t)^{2}.$$
This is what can be obtained by integrating the Hamilton equations of motion, in a discrete time domain.
If we recall the second law of Newton $\vec{F_{i}} = m_{i} \vec{\ddot{r}}_{i}$ and $\vec{\dot{r}}_{i}  = \vec{v}_{i}$, the equation above can be written as
$$\vec{r}_{i} (t+ \Delta t) = \vec{r}_{i}(t) + \vec{v}_{i} \Delta t + \frac{1}{2}\frac{\vec{F}_{i}}{m_{i}}(\Delta t)^{2}.$$
By integrating the equation of motion backward in time:
$$ \vec{r}_{i} (t- \Delta t) = \vec{r}_{i}(t) - \vec{v}_{i} \Delta t + \frac{1}{2}\frac{\vec{F}_{i}}{m_{i}}(\Delta t)^{2}$$
and if we sum the two equations above we obtain:
$$ \vec{r}_{i} (t+ \Delta t) = 2\vec{r}_{i}(t) - \vec{r}_{i} (t- \Delta t)+ \frac{1}{2}\frac{\vec{F}_{i}}{m_{i}}(\Delta t)^{2}.$$
It can be observed that in the equation above only even terms appear. Hence, the first correction is at the fourth order and not at the third one as it would have been without summing the two equations.

The velocity can be obtained by taking the difference between $\vec{r}_{i}(t+ \Delta t)$ and $\vec{r}_{i} (t- \Delta t)$:
$$\vec{v}_{i}(t) = \frac{\vec{r}_{i} (t+\Delta t) - \vec{r}_{i} (t- \Delta t) }{2 \Delta t}$$
The velocity results to be dependent on the temperature since it is related to the kinetic energy.

### Velocity Verlet algorithm (1982)
The main difference between this algorithm and the previous one is that the positions and the velocities are updated at the same time.
The starting point is the Taylor expansion of the position truncated to the second order:
$$ \vec{r}_{i} (t+ \Delta t) = \vec{r}_{i}(t) + \vec{\dot{r}}_{i}(t) \Delta t + \frac{(\Delta t)^{2}}{2m_{i}}\vec{F}_{i}(t).$$
Due to the time reversibility, which is guaranteed by the Hamiltonian formalism, it is possible to write
$$ \vec{r}_{i} (t) = \vec{r}_{i}(t+\Delta t) + \vec{\dot{r}}_{i}(t+\Delta t) \Delta t + \frac{(\Delta t)^{2}}{2m_{i}}\vec{F}_{i}(t+\ Delta t).$$
If $\vec{r}_{i} (t+ \Delta t) $ in (\eqref{vv_1}) is substituted in (\eqref{vv_2}) the following relation is obtained:
$$\vec{v}_{i} (t+ \Delta t) = \vec{v}_{i}(t) + \frac{(\Delta t)^{2}}{2m_{i}}[\vec{F}_{i}(t+\ Delta t)+\vec{F}_{i}(t+\ Delta t)].$$
To run a Molecular dynamics simulation some parameters are needed.
- Time step $\Delta t$
- Initial positions, which are provided by crystallographic experiments and collected in the PDB file
- Initial velocities, which are randomly taken from a Boltzmann distribution at temperature $T$. The system is characterized by a constant energy $E$, so it can be described using the microcanonical ensemble. In this ensemble the temperature is not constant and will change during the evolution of the system.

### Constraints
The algorithms we presented above do not take into account constraints. If reaction forces are present the formalism changes.
In our system, constrains are represented by chemical bonds. Such bonds oscillate with a high frequency; to resolve these fast oscillations a time step $\Delta t < 1 fs$ is needed. For our purposes a $\Delta t \simeq 1 fs$ is better because it allows us to run longer simulations. The chemical bonds will be seen as a steady bound by the Verlet algorithm, and they will constitute a constraint.
The Lagrange equations of motion obtained by means of the variational method if the constraint can be written in a differential form
$$\sum_{\alpha = 1}^{3N}a_{k \alpha} dq_{\alpha}=0.$$
This condition is verified for holonomic time independent constraints, not always for non-holonomic constraints. The fixed chemical bonds belong to the first category, so the generalized Euler-Lagrange equations are:
$$\frac{d}{dt} (\frac{\partial \mathcal{L}}{\partial \dot{q}_{\alpha}}) - \frac{\partial \mathcal{L}}{\partial q_{\alpha}} = - \sum_{k=1}^{N_{c}} \lambda_{k} a_{k \alpha}. $$
Here the $\lambda_{k}$ are the Lagrange multiplicators, a set of $N_{c}$ unknown numbers we have to find to apply the Verlet algorithm.
To find them another condition is necessary, namely the [TODO]

## Project structure
The crystallographic structures are available for many mutations of the $M^{pro}$, but the dynamics of the protein has not been explored yet. The aim of this project is to explore the dynamics of the $M^{pro}$, in particular of the [8DFN mutation](https://www.rcsb.org/structure/8DFN).

### Part 1: Molecular Dynamics Simulation
To explore the role of the mutation in the dynamics of the protein, the wild type structure is compared to the mutated one. In this case, the protein from the [7SI9 complex](https://www.rcsb.org/structure/7SI9) was taken as the wild type reference, after discarding the additional ligand present in the PDB file.

Given that the protein is composed of two monomers, we aimed to investigate how the quaternary interactions of the protein affects the dynamics of the subunits. For achieving this, each full protein evolution is compared to the evolution of its respective single subunit. As a result, a single simulation run consists of 4 trajectories.

During the evolution, the system does not explore all the accessible microstates. Hence, starting from the same initial configurations, different trajectories are obtained. Because of this, we ran two repetitions for each trajectory and compared them.

In summary, by following the workflow of this project 8 different trajectories are produced. To avoid confusion, a short nomenclature was employed to name the different trajectories, which is described as follows:

| Trajectory name | protein | subunits | repetition |
| :-----: | :-----: | :-----: | :-----: |
| **mt1_rep0** | mutated-type | 1 | first  |
| **mt1_rep1** | mutated-type | 1 | second |
| **mt2_rep0** | mutated-type | 2 | first  |
| **mt2_rep1** | mutated-type | 2 | second |
| **wt1_rep0** | wild-type    | 1 | first  |
| **wt1_rep1** | wild-type    | 1 | second |
| **wt2_rep0** | wild-type    | 2 | first  |
| **wt2_rep1** | wild-type    | 2 | second |


### Part 2: General Analysis
In the general analysis we focus on the equilibrium configurations and on the conformational changes of all the trajectories. The methods shown are **Root Mean Square Deviation (RMSD)**, **Root Mean Square Fluctuation (RMSF)**, **Radius of Gyration (RGYR)** and **Contact Map (CMAP)**.

### Part 3: Detailed Analysis
The following part of the analysis focuses mainly on one repetition of the mutated type protein, **mt2_rep0**.
The results from the methods introduced in the previous section are used as a basis to choose frames of interest from the trajectory, which are further studied. The residues from and near the active site are also explored in more detail.

For this section, the usage of RMSD values is expanded by introducing clustering methods, namely **Hierarchical Clustering** and **Principal Components Analysis (PCA)**. Additionally, an enhanced way of visualizing CMAP values is presented, denominated here **Dynamic Contact Map Differences (DCMD)**. The **Block Standard Error Analysis (BSE)** and the **Pyinteraph** methodology are also introduced during the analysis.

<!-- When we observed its evolution we observed how the chain around the active site opens, and at the same time in the upper part the chain closes on itself -->

### Part 4: Solvent Analysis
The last part of the analysis focuses on the behaviour of water during the simulation. We aimed to calculate the density of water molecules, in relation with the hydrophobic and hydrophilic residues of the protein. For achieving this, a new methodology named **Water Average Density (WAD)** is proposed. The analysis is further complemented by means of the **Solvent Accessible Surface Area (SASA)** command provided by GROMACS.
