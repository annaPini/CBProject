# Introduction
## Physical system
The Sars-Cov-2 is organised in two polyproteins which contain shorter proteins, among which two main protease proteins, M^pro^. They play a fundamental role in the replication of the virus, by cutting the polyproteins in shorter chains. Due to its function and because it is absent in the human body, M^pro^ is one of the most important drug targets. 
The protein is composed by 306 amino acids and comprehend three domains:
- Domain I: residues 8-101
- Domain II: residues 102-184
- Domain III: residues 201-306
There is also a long loop, residues 185-200, which connects domains II and III. 
From the crystallographic structure the M^pro^ results to be divided in two subunits, inactive monomers. Each subunit has an active site, known as the catalytic dyad, composed of the residues H41-C145. 
## Molecular dynamics
As mentioned before, all the simulations and the analysis are carried out by means of Molecular dynamics simulation. It is then worth to briefly introduce this method. In the system that we are considering there are not dissipative forces, then Hamiltonian formalism can be used. The equations of motion are obtained by taking the time derivative of the generalized coordinates and momenta
$$\dot{q}_{\alpha} = \frac{\partial \mathcal{H}}{\partial p_{\alpha}}  \dot{p}_{\alpha} =- \frac{\partial \mathcal{H}}{\partial q_{\alpha}} $$
These are first order differential equations and the algorithms we are going to introduce are meant to numerically solve them. The positions and the velocity of each atom have to be found as a function of time and by keeping the total energy constant. This means the total derivative of the Hamiltonian is zero: (derivate lungo I moti)
$$\frac{d \mathcal{H}}{dt} = 0 $$
Our system is also ergodic, meaning that during its evolution the system visits almost all the microstates which belong to the ensemble. If we apply the algorithm to integrate the equations of motions many times, we will always find different results. Due to ergodicity we know that the error does not diverge. This property is quite relevant because it allows to exactly evaluate the value of macroscopic parameters by means of an average on time. For our purposes, a discrete average is necessary, since the it is not possible for a computer to use a continuous variable
$$ A = \braket{a} = \frac{1}{M}\sum_{n=1}^{M} a(x_{n \Delta t}) $$
Where $\Delta t$ is the time step and $ a(x_{n \Delta t})$ is the macroscopic parameter evaluated at each time step on each element of the ensamble which is visited by the system.
\su## Basic components of a Molecular dynamics simulation
\begin{itemize}
\item The Model, that is the force field
\item Calculation of forces and energies
\item The algorithm to integrate the equations of motion
\end{itemize}
### Verlet algorithm (1967)
The Verlet algorithm considers the Taylor expansion to the second order of the position of the atom at $(t+ \Delta t)$
$$\vec{r}_{i} (t+ \Delta t) = \vec{r}_{i}(t) + \vec{\dot{r}}_{i} \Delta t + \frac{1}{2}\vec{\ddot{r}}_{i}(\Delta t)^{2}.$$
This is what can be obtained by integrating the Hamilton equations of motion, in a discrete time domain.
If we recall the second law of Newton $\vec{F_{i}} = m_{i} \vec{\ddot{r}}_{i}$ and $\vec{\dot{r}}_{i}  = \vec{v}_{i}$, the equation above can be written as follow
\vec{r}_{i} (t+ \Delta t) = \vec{r}_{i}(t) + \vec{v}_{i} \Delta t + \frac{1}{2}\frac{\vec{F}_{i}}{m_{i}}(\Delta t)^{2}.$$
By integrating backward in time we obtain
$$ \vec{r}_{i} (t- \Delta t) = \vec{r}_{i}(t) - \vec{v}_{i} \Delta t + \frac{1}{2}\frac{\vec{F}_{i}}{m_{i}}(\Delta t)^{2}.$$
If we sum the two equations above we obtain:
$$ \vec{r}_{i} (t+ \Delta t) = 2\vec{r}_{i}(t) - \vec{r}_{i} (t- \Delta t)+ \frac{1}{2}\frac{\vec{F}_{i}}{m_{i}}(\Delta t)^{2}.$$
It can be observed that in the equation above only even terms appear. This means that the first correction is at the fourth order and not at the third one as it would have been without summing the two equations. 
The velocity can be obtained by taking the difference between \eqref{-\Delta t} and \eqref{\Delta t}:
$$\vec{v}_{i}(t) = \frac{\vec{r}_{i} (t+\Delta t) - \vec{r}_{i} (t- \Delta t) }{2 \Delta t}$$ 
It is worth noticing than the velocity results to be dependent on the temperature since it is related to the kinetic energy.
\su## Velocity Verlet algorithm (1982)
This is another version of the Verlet algorithm described above. The main difference is that in this case the positions and the velocities are updated at the same time. 
The starting point is the Taylor expansion of the position truncated to the second order
$$ \vec{r}_{i} (t+ \Delta t) = \vec{r}_{i}(t) + \vec{\dot{r}}_{i}(t) \Delta t + \frac{(\Delta t)^{2}}{2m_{i}}}\vec{F}_{i}(t).$$
Due to the time reversibility, which is guaranteed because these equations are derived from Hamilton equations of motion, it is possible to write
$$ \vec{r}_{i} (t) = \vec{r}_{i}(t+\Delta t) + \vec{\dot{r}}_{i}(t+\Delta t) \Delta t + \frac{(\Delta t)^{2}}{2m_{i}}\vec{F}_{i}(t+\ Delta t).$$
If $\vec{r}_{i} (t+ \Delta t) $ in (\eqref{vv_1}) is substituted in (\eqref{vv_2}) the following relation is obtained:
$$\vec{v}_{i} (t+ \Delta t) = \vec{v}_{i}(t) + \frac{(\Delta t)^{2}}{2m_{i}}[\vec{F}_{i}(t+\ Delta t)+\vec{F}_{i}(t+\ Delta t)].$$
To run a Molecular dynamics simulation some parameters are needed. 
\begin{itemize}
\item Time step $\Delta t$
\item Initial positions, which are provided by crystallographic experiments and collected in the PDB file
\item Initial velocities, which are randomly taken from a Boltzmann distribution at temperature T. It is important noticing that the system we consider is characterized by a constant energy $E$, then it can be described using the microcanonical ensemble. In this ensemble the temperature is not constant and will change during the evolution of the system. 
\end{itemize}
## Constraints
The algorithms we presented until now are derived from Hamilton equations of motion in case no constraints are present. If reaction forces are present the formalism becomes a little more complicated. 
In our system the constrains are represented by the chemical bonds. Such bonds oscillate with an high frequency, then to resolve these fast oscillations a time step $\Delta t < 1 fs$ is needed. For our purposes a $\Delta t \simeq 1 fs$ is better because it allows to run longer simulations. The chemical bounds will be seen as a steady bound by the Verlet algorithm, and they will constitute a constraint. 
The Lagrange equations of motion obtained by means of the variational method if the constraint can be written in a differential form
$$\sum_{\alpha = 1}^{3N}a_{k \alpha} dq_{\alpha}=0.$$
This condition is verified for holonomic time independent constraints, not always for non-holonomic constraints. The fixed chemical bonds belong to the first category, then the generalized Euler-Lagrange equations are: 
$$\frac{d}{dt} (\frac{\partial \mathcal{L}}{\partial \dot{q}_{\alpha}}) - \frac{\partial \mathcal{L}}{\partial q_{\alpha}} = - \sum_{k=1}^{N_{c}} \lambda_{k} a_{k \alpha}. $$
Here the $\lambda_{k}$ are the Lagrange multiplicators, a set of $N_{c}$ unknown numbers we have to find in order to apply the Verlet algorithm. 
To find them another condition is necessary, namely the 
### Project structure
The crystallographic structures are available for many mutations of the M^{pro}, but the dynamics of the protein has not been explored so far. The aim of this project is to explore the dynamics of the main protease, M^{pro}, in particular the 8DFN mutation. 
To explore the role of the mutation in the dynamics of the protein, the wild structure is compared to the mutated one. Given that the protein is composed of two dimer, we tried to investigate how the dynamic of one subunit depends on the dynamics of the other, the full protein evolution is compared to the evolution of a single subunit. As a result, four trajectories are analyzed. 
The 
We focused on the active site, and we explored its evolution and fluctuations. 
We also focused on water, and we tried to understand why the TIP3P water model and the  force field were used. We also investigate how the water molecule arrange around the protein in relation with the hydrophilic and hydrophobic regions of the protein itself. The dynamic of water molecules was also briefly explored. We qualitatively investigated how water goes inside the protein. 
