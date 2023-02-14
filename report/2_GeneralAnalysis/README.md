The source code for the general, detailed and solvent analysis is found in the folder ```src_analysis```. It contains scripts for both calculation and visualization of the results. Numerical results are stored in a folder created while executing ```src_md/3_POST.sh```, called ```data_analysis```.

# RMSD
From the trajectories visualized in VMD it can be seen how the protein fluctuates during its evolution in water. The first frame, which is the PDB file taken from crystallographic experiments, represents the protein in a structure far from equilibrium.
During its evolution, the system should move towards a minimum of the free energy, and reach more stable configurations.
We wonder if the system reaches an equilibrium configuration, and if there are some transitions from a metastable configuration to another.
To explore this behaviour we can calculate the **root mean square deviation (RMSD)** for each trajectory.  
The RMSD between two frames is defined as

$$ RMSD(\Delta t) = \sqrt{\frac{1}{N}\sum_{i=1}^{N}(x_{i}(\Delta t)-x_{i}(0))^{2}}$$

where
- $N$ is the number of atoms in the protein
- $\Delta t$ is the time step
- $x_{i}(\Delta t)$ is the position of atom $i$ at time $\Delta t$

It is important to remember that this quantity is only an indicator of the stability of the system. If a plateau is reached it means that the protein can be at equilibrium, but it important to notice that this is not always the case.
Furthermore, before calculating the RMSD it is necessary to perform an alignment between the structures. In our case, this was performed in the post-processing described before.

## RMSD calculation
Calculating a RMSD map of ```6000 x 6000``` frames is computationally expensive. Firstly, one must consider that large XTC file should be opened as a stream (loading just the frame needed for a given operation) instead of loading all frames in memory. This means that, when calculating the RMSD matrix, the XTC "reading window" should shift from start to finish multiple times to load the required frames, which is highly demanding.

Considering this, the function ```calculate_general.extract_ca_coords``` was defined. It iterates over the XTC once, assigning the coordinates of the backbone carbon alpha atoms for each frame into a numpy array, then saving the whole array in ```data_analysis/general```. This data structure is much lighter and can be loaded into memory for any method that requires using the backbone positions of the protein.

Furthermore, as the RMSD matrix is symmetric, the function ```calculate_general.calc_rmsd``` calculates only half of the matrix, and then assigns the other half in a specular manner. This structure is also stored in ```data_analysis/general```, to avoid calculating the RMSD each time a different visualization is performed. This idea of pre-calculation is also applied for other subsequent analysis methodologies.

## RMSD visualization - 1D
A frame from the trajectory can be taken as a reference point for the RMSD with other frames, as a means of general comparison between conformations and their evolution. By taking the respective row from the RMSD matrix, a fixed 1D representation can be obtained. The classes ```RMSD_visualize.RMSD_1D``` and ```RMSD_visualize.RMSD_1D_Compare``` were defined with this purpose in mind (the latter allows for plotting two RMSD lines simultaneously).

In the following plots, the RMSD is visualized with respect to the initial frame, to estimate if the proteins reach an equilibrium during the evolution of the system (i.e. a plateau on RMSD values with respect to the original structure). Furthermore, each plot contains both simulation repetitions, which allows to directly compare them. In general, the noise is quite high, and it is not straightforward to identify the stable or equilibrium regions on the plots.

![RMSD1D-mt1.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD1D-mt1.png)
The first system, ```mt1```, presents large differences in their RMSD values. An important disturbance seems to occur around frame 900 for $rep0$ (but not for $rep1$) and a change in the trend of configurations can be noticed. Always in $rep0$ we observe an equilibrium region approximately between frame 2100 and frame 3800. However, before frame 2000 an extremely noisy region occurs, and after frame 4000 it seems that the configuration of the protein changes. Meanwhile, in $rep1$ the RMSD fluctuates around a constant value and is much less noisy than the $rep0$ one. The only relevant change occurs between frame 2000 and frame 3000, which might be a slight conformational change.

All classes inherited by ```RMSD_visualize.RMSD_1D``` also allow to "slide" the reference frame in an interactive manner. This allows for a better understanding of the whole RMSD matrix, as exemplified bellow.

![RMSD1D-mt1.gif](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD1D-mt1.gif)

The $rep1$ values slightly oscillate around a constant value, while the values for the $rep0$ show a more dynamic trend. It can be observed that when the reference frame is approximately between frame 2100 and frame 5000 the two RMSD tend to coincide in the central region. In general, the stability of this system is doubtful, so observing its wild-type analogous ```wt1``` perhaps provides more information.

![RMSD1D-wt1.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD1D-wt1.png)

The RMSD of $rep0$ is quite low and it seems to oscillate around a constant value. It seems that no major changes of configuration occur. On the other side, the RMSD of $rep1$ undergoes several changes. At frame 800 it changes abruptly, perhaps consequence of a similar disturbance as suffered by mt1_rep0. Until frame 4000 it seems that the protein stays in the same equilibrium configuration, only to then present extreme changes on the RMSD values.

Overall, the system with only one subunit seems to be highly unstable, sometimes apparently reaching equilibrium quickly while other times undergoing extreme changes because of slight perturbations. This is most likely related to the fact that this system does not operate in nature, as the protein is expected to be constituted by two monomers. Perhaps this is so because having two monomers provides structural stability to the polypeptides, and could explain why this is the form found in nature. Effectively, when checking the evolution of the dimeric systems (```wt2```, ```mt2```), a higher stability is observed:

![RMSD1D-wt2.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD1D-wt2.png)

![RMSD1D-mt2.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD1D-mt2.png)

It is noticeable how stable the ```wt2``` system is with respect to the other trajectories. However, this should not be a surprise, as this system is the wild-type dimeric form prevalent in nature. The ```mt2``` system is also quite stable, with seemingly only one relevant conformational change in both runs. This change occurs between frame 4000 and 6000, and, as it comes from the mutated system and its explanation is not evident, will be further analysed in [Detailed Analysis](03-Detailed-Analysis).

## RMSD visualization - 2D
It is also possible to build a RMSD matrix and represent it using an heatmap. The RMSD graphs represented above correspond to the first row in the respective matrix.

The matrix are symmetric, and the values on the diagonal are $0$.

![RMSD2D-mt1.jpg](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD2D-mt1.jpg)

In $rep0$ map it seems a change in configuration occurs around frame 800. Other changes inconfiguration occurs but the contrast is lower. In $rep1$ map the contrast is much lower, meaning that the protein does not undergoes major changes in its configurations. It is interesting to observe that two changes in the protein structure occur around frame 2500. In the 1D graph we could not grasp this information, because the RMSD are noisy.

![RMSD2D-mt2.jpg](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD2D-mt2.jpg)

In $rep0$ only a change in the configuration occur around frame 5500. In $rep1$ several changes occurs, and the value of RMSD is in general higher than the ones in $rep0$.

![RMSD2D-wt1.jpg](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD2D-wt1.jpg)

Here we observe two deeply different behaviour in the RMSD. In $rep0$ the values are substantially constant and it seems that no important changes in the configuration occur. In $rep1$ a clear change in the configuration occurs around frame 4000.

![RMSD2D-wt2.jpg](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSD2D-wt2.jpg)

The RMSD values in both repetitions seems to be similar. In $rep1$ a change in conformation seems to occur around frame 3000.

# RMSF
The protein is composed of many atoms which move during the evolution of the system. To understand which parts of the protein are more mobile it is useful to calculate the root mean square fluctuations for each atom.
The RMSF related to atom $i$ is defined as:
$$ RMSF_{i} = \sqrt{<\Delta \vec{r_{i}}>^{2}} $$
where $<\Delta \vec{r_{i}}>$ is the mean deviation of the atom from its equilibrium position.

In the graphs below it can be seen that some atoms have a larger RMSF, meaning that they are more mobile during the evolution of the protein.

![RMSF-mt1.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSF-mt1.png)

The RMSF values for $rep0$ are in general larger than the values calculated for $rep1$. This means that the atoms oscillate more around their equilibrium position in $rep0$ than in $rep1$.
The RMSF graphs are compatible with the RMSD ones. As regards $rep0$, we can observe non constant regions, which could correspond to changes in the conformation of the protein. In these regions, the atoms are expected to move more, than we expect an higher value of the RMSF. On the other hand, the RMSD values calculated for $rep1$ seem to oscillate around a constant value. This means that probably conformational changes do not occur, then the atoms are expected to slightly move around their equilibrium position.

![RMSF-mt1_rep0-mt2_rep0.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSF-mt1_rep0-mt2_rep0.png)

Here we present the RMSF of the whole protein with the one computed for only a subunit. In this way we can explore how the fluctuations of atoms in a given subunit are correlated to the fluctuations in the other subunit.
We observe that for the protein with two subunits the RMSF is lower than the one calculate for one subunit alone. This means that the system with two subunits is more stable than the system with only one subunit.


![RMSF-mt2.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSF-mt2.png)

![RMSF-mt2_rep0-wt2_rep0.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSF-mt2_rep0-wt2_rep0.png)

![RMSF-wt1.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSF-wt1.png)

![RMSF-wt1_rep0-wt2_rep0.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSF-wt1_rep0-wt2_rep0.png)

![RMSF-wt2.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RMSF-wt2.png)



# Radius of Gyration
From VMD it can be observed that even though the protein fluctuates, its shape remains constant. To check if the protein elongates during its evolution it is useful to calculate the radius of gyration:
$$ RGYR(\Delta t) = \sqrt{\frac{\sum_{i=1}^{N}m_{i}(r_{i}(\Delta t)-R_{cm})^{2}}{\sum_{i=1}^{N} m_{i}}}$$
where
- $m_{i}$ is the mass of atom $i$
- $r_{i}$ is the position of atom $i$ at time $\Delta t$
- $R_{cm}$ is the center of mass

In the following graphs we present the radii of gyration for the four system. The noise is large, then few assuptions can be done as regard the compatibility between the two sets of data and the configurational changes of the protein. The RMSF computed for $rep0$ amd $rep1$ seem to agree more for $mt2$ and $wt2$, then for the evolution for the whole protein. This could mean that a single subunit is less stable than two subunits linked togheter.

![RGYR-mt1.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RGYR-mt1.png)

The radius of gyration, RGYR, for $rep0$ presents a peak around frame 900, and two local maxima around frame 2000 and frame 4000. A non constant RGYR could mean that the protein is changing shape, then it is compatible with a change of the configuration. The RGYR for $rep1$ oscillates around $220nm$ and it slightly increases between frame 1000 and frame 3800. In this interval it is possible that occurs a conformational change.

![RGYR-mt2.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RGYR-mt2.png)

The values of the RGYR are extremely noisy, then it is quite difficult to understand if a change of conformation occurs. In $rep1$ a peak is observed around frame 4800, and in $rep0$ the same peak occurs around frame 5500. We think that the protein assumes the similar configurations around frame 4800 and 5500.

![RGYR-wt1.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RGYR-wt1.png)

The values of the RGYR quite agree approximately until frame 1800 and starting from frame 5500 to the end. Starting from frame 2000 the RGYR of $rep1$ increases and reaches a peak around frame 4000.

![RGYR-wt2.png](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/RGYR-wt2.png)

The RGYR of the two repetitions quite agree, it seems they both evolve with a sinusoidal trend.

# Contact Map
The protein is constituted by a folded chain composed of atoms linked together. Non subsequent atoms can be near enough to interact even if they are not adjacent. The contact maps shows the links between atoms only for a single frame, and a more detailed analysis of interactions in the protein is performed in the following part of this work.
The contact map is defined in 2D, then it is not possible to overlap the graphs as we did for RMSF, RMSD and RGYR.
Here we show as an example a contact maps calculated at the initial frame of the trajectory. The initial frame is the same for the two trajectories, namely rep0 and rep1, then the contact map results to be the same.

![CMAP-mt1mt2.jpg](https://github.com/annaPini/CBProject/blob/main/report/2_GeneralAnalysis/CMAP-mt1mt2.jpg)

The most dark parts represent the regions in which there is more contact between the residues, while the brightest parts correspond to the regions in which there is less contact. In the contact map on the right a discontinuos line in the indicates the region in which the first subunit ends and the secon one begins. It can be observed that the map is symmetric, and for the whole protein the pattern is repeated two times along the main diagonal.
