# Computational Biophysics Project
This documentation serves as a report for the project of Computational Biophysics course 2022-2023 (University of Trento). The repository contains the source code and initial input files necessary to reproduce all the steps detailed by the report, as well as instructions of use inside the documentation.

The focus of this project is to analyse the dynamic of a main protease protein. The software we used to simulate the trajectory is [GROMACS](https://www.gromacs.org/), which integrates the equations of motion by means of the Verlet algorithm. Once the trajectories were obtained we analysed them, focusing on general aspects of the trajectories, as well as frames of interest, the active sites and the distribution of water around the protein. The analysis is mostly performed via Python 3 scripts. The library [MDAnalysis](https://www.mdanalysis.org/) facilitates dealing with GRO and XTC files, as well as performing the different analysis methodologies.

# Index
- [Introduction](00-Introduction)
- [Part 1: Molecular Dynamics](01-Molecular-Dynamics)
- [Part 2: General Analysis](02-General-Analysis)
- [Part 3: Detailed Analysis](03-Detailed-Analysis)
- [Part 4: Solvent Analysis](04-Solvent-Analysis)
- [Conclusions](05-Conclusions)
