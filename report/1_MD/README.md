# Molecular dynamics simulation
We were provided the PDB files for the 8DFN mutation main protease protein of the Sars-Cov2. The first step of the MD simulation were performed, namely
- the creation of the topology file
- the construction of the simulation box 
- the ions were added

At this point the system is complete but the protein is still in a high energy configuration. Hence it is necessary to find the configuration associated to the minimum value of the energy and then let it evolve. 

We set some of the parameters necessary to perform the simulation:
- the temperature at $310K$ 
- the lenght of the time step $0.2ns$
- the number of steps, hence the simulation period at $300 ns$

The entire MD simulation was performed on the cluster of the university. To choose the right number of cores to use we did some benchmarks. We runned the GROMACS command with a walltime of $10$ and $30$ minutes. We also repeated the benchmarks two times for each walltime. We noticed that the efficiency obtained with the same number of cores and walltime was not constant. This happend because the command is runned on different cores every time. 

The first step we performed was the **energy minimisation**. The protein explored many configurations to find the one associated to the minimum value of the molecular energy. 

An equivalent method through which this result can be achieved is to find the volume and the pressure values which correspond to the minimum of the molecular energy of the protein. 

The first step consists in the **NVT equilibration** the volume is fixed, and the pressure slightly changes, until the minimum energy configuration is found. The value of the pressure which corresponds to the minimum energy configuration will be employed in the following step. 

The second step is the **NPT equilibration**. The volume slightly changes, until the configuration which corresponds to the minimum energy is found. 

At this point we can proceed with the MD unrestrained simulation. 

The last step is the **post processing**. The periodic boundary condition and the diffusive dynamics are removed. We also had to center water. During our first attempt the water left the box while the protein was evolving. We also had to prevent the protein from exit from the box. In this case its evolution could not be represented in the correct way. 