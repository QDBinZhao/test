![LOGO](https://github.com/QDBinZhao/test/blob/main/logo.png)
# ABR
ABR is a program suite for studying nuclear quantum dynamics in polyatomic molecular systems.
# ABR lite
The development of the ABR program is still ongoing, and the current release is a lite version that is tailored for studying the non-adiabatic quenching of OH(A) by collisions with H<sub>2</sub>,<br>
 OH(A) + H<sub>2</sub> &#8594; H + H<sub>2</sub>O (reactive quenching)<br>
                     &#8594; OH(X)+H<sub>2</sub> (non-reactive quenching)<br>
                     &#8594; OH(A)+H<sub>2</sub> ((in)elastic scattering)

This lite version can perform following calculations:
1. Select specific initial ro-vibrational states of the OH and H<sub>2</sub> reactants.
2. Calculate the fractional probabilities of the three channels.
3. Resolve the final state information in the non-reactive quenching channel and the (in)elastic scattering channel.

---

# Features
The ABR program suite is developed to study the ABR type of molecular systems.
Here, A and B represent two idividual atoms, while the R group represents loosely for a radical of a complete molecule (not necessarily a free radical), such as a hydroxyl group in water or a methyl group in methane.  

The coordinate system used in the description of the ABR system consists of three parts:  
1. two Jacobi vectors connecting the A, B, and the center of mass (COM) of the R group;
2. the internal coordinates of the R group;
3. the rotational angles that rotate to the R-fixed frame.

## The body-fixed frame that is determined by two Jacobi vectors
![coordinates](https://github.com/QDBinZhao/test/blob/main/coordinates.png)  
The two Jacobi vectors, ___r___<sub>i</sub> and ___R___<sub>i</sub>, i=1,2,3, shown in the above figure, determine the body-fixed (BF) frame of the system. The z axis lies along the ___R___<sub>i</sub> vector, and the ___r___<sub>i</sub> vector lies in the x-z plane of the BF frame and points in the positive direction of the x axis.

## The internal coordinates of the R group
For __three-atom systems__, the R group is an atom and has no internal nuclear degree of freedom (DOF).  
For __four-atom systems__, the R group is a two-atom linear radical/molecule
and has one internal nuclear DOF.  
For __six-atom systems__, the R group is a four-atom radical
and has six internal nuclear DOFs.
The description of the internal motion of the R group can use part or all of the six DOFs.  
For other polyatomic systems consisting of more than four atoms, similar part or full internal coordinates can be used to describe the internal motion of the R group.

## The rotational angels that rotate to the R-fixed frame
For __three-atom systems__, the single atom of the R group does not form a frame so that no rotational motion of the R group is present.  
For __four-atom systems__, description of the rotation to the linear R-fixed frame requires two Euler angles and the third one describing the rotation around the linear axis is absent.  
For __six-atom systems__, three Euler angles are required to describe the rotation to a non-linear R-fixed frame. If the R group is linear, only two Euler angles are required.  
Other polyatomic systems consisting of more than four atoms follow the six-atom systems.

---

# Field of research
The ABR program is used to solve the time-dependent and time-independent [Schr&ouml;dinger equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) that governs the wave functions of the quantum-mechanic nuclear motion.
The study relies on the [Born-Oppenheimer (BO) approximation](https://en.wikipedia.org/wiki/Born%E2%80%93Oppenheimer_approximation) that assumes the separability of the nuclear and electronic motion (due to the large mass difference between nuclei and electron).
The ABR program performs calculations of nuclear quantum dynamics on pre-existing analytic [potential energy surfaces (PESs)](https://en.wikipedia.org/wiki/Potential_energy_surface).
## What can't ABR program do
- The ABR program can't perform calculations for electronic motion.
- The ABR program can't fit _ab initio_ potential energies to an analytic form.

## What can ABR program do  
- The ABR program can describe time-dependent nuclear motion on a single PES or on coupled PESs.
- The ABR program can calculate eigenstates of the molecular system.

---

# What types of calculations can ABR program perform?
The ABR program can perform following types of jobs for a few molecular systems.

## Types of molecular systems
The _kcoord_ controller selects the type of molecular systems:
- kcoord = 1, 2, and 3 for three-atom molecular systems;
- kcoord = 13 and 22 for four-atom molecular systems;
- kcoord = 93 for three-atom molecular systems on a surface;
- kcoord = 15 and 24 for six-atom molecular systems;
- kcoord = 95 for five-atom molecular systems on a surface.

## Types of jobs
The _idtype_ controller selects the type of jobs:
- idtype = 0 for a single channel calculation from a single initial state
   - calculate total reaction probability using a flux operator placed right after the transition state;
   - calculate the total (in)elastic scattering probability using a flux operator placed in the reactant asymptotic region;
   - calculate the final state-resolved scattering probabilities and _S_-matrix elements in the (in)elastic channel.
- idtype = 3 and 31 for calculating state-to-state reaction probabilities and _S_-matrix elements using the product coordinate based (PCB) method within the initial state-selected wave packet (ISSWP) framework.
   - idtype = 3 for propagating the initial (focusing) wave packet close to the interaction region and store the overall wave function on disk for the second step;
   - idtype = 31 for continuous propagation of the wave packet in the product coordinate system after a wave function representation change from the reactant coordinates to the product coordinates. Final product state-resolved reaction probabilities and _S_-matrix elements can be obtained on a dividing surface in the product asymptotic region.
- idtype = 4 and 41 for calculating state-to-state reaction probabilities and _S_-matrix elements using the reactant coordinate based (RCB) method within the initial state-selected wave packet (ISSWP) framework.
   - idtype = 4 for propagating the initial wave packet all the way to the product asymptotic region and interpolating the wave function on the grids of a dividing surface in the product asymptotic region.
   - idtype = 41 for resolving final product state information.
- idtype = 5 and 51 for calculating state-to-state reaction probabilities and _S_-matrix elements using the reactant product decoupling (RPD) method within the initial state-selected wave packet (ISSWP) framework.
   - idtype = 5 for calculations that are the same as idtype = 0 except that the wave functions removed by the absorbing potential in the product channel is stored on the disk for the continuous calculations in the second step;
   - idtype = 51 for continuous propagation of the wave packet in the product coordinate system after wave function representation changes from the reactant coordinates to the product coordinates. Final product state-resolved reaction probabilities and _S_-matrix elements can be obtained on a dividing surface in the product asymptotic region.
- idtype = 6, 61, and 62 for the calculations in the quantum transition state (QTS) framework
   - idtype = 6 for the calculation of thermal flux eigenstates;
   - idtype = 61 for the propagation of thermal flux eigenstates in the reactant coordinates. Initial state-resolved total reaction probabilities of the reactant arrangement channel can be obtained on a dividing surface in the reactant asymptotic region.
   - idtype = 62 for the propagation of thermal flux eigenstates in the product coordinates after a wave function representation change from the reactant coordinates to the product coordinates. Initial state-resolved total reaction probabilities of the product arrangement channel can be obtained on a dividing surface in the product asymptotic region.
   - The whole _S_-matrix at a selected number of energy grids can be obtained combining the state-resolved information in the calculations with idtype = 61 and 62.
- idtype = 7 for the calculation of eigenstates of the whole system.

---

# Example systems
The ABR program has been employed to study the following reaction systems.
## Non-adiabatic quantum reaction dynamics on coupled potential energy surfaces
1. Non-adiabatic transitions in the F(<sup>2</sup>P) + CHD<sub>3</sub> &#8594; HF + CD<sub>3</sub> reaction on coupled PESs including vibronic and spin-orbit couplings.
   - Bin Zhao and Uwe Manthe, “Non-adiabatic transitions in the reaction of fluorine with methane.” [_J. Chem. Phys._ __152__, 231102 (2020)](https://aip.scitation.org/doi/10.1063/5.0013852).
2. Non-adiabatic quenching of OH(A) by collisions with H<sub>2</sub>.
   - Bin Zhao, Shanyu Han, Christopher Malbon, Uwe Manthe, David Yarkony and Hua Guo, “Full-dimensional quantum stereodynamics of the nonadiabatic quenching of OH(A<sup>2</sup>&#931;<sup>+</sup>) by H<sub>2</sub>”, accepted for publication in _Nature Chem._ (2021)

## State-to-state study within the Quantan transition state (QTS) framework
1. Integral and differential cross sections in the H + D<sub>2</sub>  &#8594; HD + D reaction.
   - Bin Zhao, Zhigang Sun, and Hua Guo, “Calculation of state-to-state differential and integral cross sections for atom-diatom reactions with transition-state wave packets”, [_J. Chem. Phys._ __140__, 234110 (2014)](https://aip.scitation.org/doi/10.1063/1.4883615).
2. State-to-state reaction probabilities in the H<sub>2</sub> + OH &#8652; H + H<sub>2</sub>O reaction:
   - Bin Zhao, Zhigang Sun, and Hua Guo, “State-to-State Mode Specificity: Energy Sequestration and Flow Gated by Transition State”, [_J. Am. Chem. Soc._ __137__, 15964 (2015)](https://pubs.acs.org/doi/10.1021/jacs.5b11404).
   - Bin Zhao, Zhigang Sun, and Hua Guo, “Calculation of the state-to-state S-matrix for tetra-atomic reactions with transition-state wave packets: H<sub>2</sub>/D<sub>2</sub> + OH &#8594; H/D + H<sub>2</sub>O/HOD”, [_J. Chem. Phys._ __141__, 154112 (2014)](https://aip.scitation.org/doi/10.1063/1.4898100).
   - Bin Zhao, Zhigang Sun, and Hua Guo, “State-to-state Dynamics of the Cl + H<sub>2</sub>O  &#8594; HCl + OH Reaction: Energy Flow into Reaction Coordinate and Transition-state Control of Product Energy Disposal”, [_J. Chem. Phys._ __142__, 241101 (2015)](https://aip.scitation.org/doi/10.1063/1.4922650).
3. State-to-state reaction probabilities in the F + H<sub>2</sub>O &#8594; HF + OH reaction
   - Bin Zhao and Hua Guo, “Modulations of Transition-State Control of State-to-State Dynamics in the F + H<sub>2</sub>O &#8594; HF + OH Reaction”, [_J. Phys. Chem. Lett._ __6__, 676 (2015)](https://pubs.acs.org/doi/10.1021/acs.jpclett.5b00071).
4. State-to-state reaction probabilities in the Cl + H<sub>2</sub>O &#8594; HCl + OH reaction.
   - Bin Zhao, Zhigang Sun, and Hua Guo, “State-to-state Dynamics of the Cl + H<sub>2</sub>O  &#8594; HCl + OH Reaction: Energy Flow into Reaction Coordinate and Transition-state Control of Product Energy Disposal”, [_J. Chem. Phys._ __142__, 241101 (2015)](https://aip.scitation.org/doi/10.1063/1.4922650).
5. State-to-state reaction probabilities and differential cross sections of the H + CH<sub>4</sub> &#8652; H<sub>2</sub> + CH<sub>3</sub> reaction.
   - Bin Zhao and Uwe Manthe, “Eight-Dimensional Wave Packet Dynamics within the Quantum Transition-State Framework: State-to-State Reactive Scattering for H<sub>2</sub> + CH<sub>3</sub> &#8652; H + CH<sub>4</sub>”, [_J. Phys. Chem. A_ __124__, 9400 (2020)](https://pubs.acs.org/doi/10.1021/acs.jpca.0c08461).
   - Differential cross sections of the H + CH<sub>4</sub> &#8594; H<sub>2</sub> + CH<sub>3</sub> reaction with different vibrational and rotational excitations in the CH<sub>4</sub> reactant, _in preparation_.
   - Differential cross sections of the H<sub>2</sub> + CH<sub>3</sub> &#8594; H + CH<sub>4</sub> reaction with different vibrational and rotational excitations in the H<sub>2</sub> and CH<sub>3</sub> reactants, _in preparation_.
6. Initial state-selected and state-to-state reaction probabilities of the H + CHD<sub>3</sub> &#8594; H<sub>2</sub> + CD<sub>3</sub> reaction.
   - Bin Zhao, “The symmetric C-D stretching spectator mode in the H + CHD<sub>3</sub> &#8594; H<sub>2</sub> + CD<sub>3</sub> reaction and its effect on dynamical modeling”, submitted to _Phys. Chem. Chem. Phys._
   - Bin Zhao and Uwe Manthe, "State-to-state reaction probabilities of the H + CHD<sub>3</sub> &#8594; H<sub>2</sub> + CD<sub>3</sub> reaction with multiple ro-vibrational excitations in the CHD<sub>3</sub> reactant.", in preparation
7. Sticking probability of the dissociative chemisorption of CH<sub>4</sub> on Ni(111) transition metal surface.
   - Bin Zhao, Uwe Manthe, and Hua Guo, "Effects of initial vibrational excitation and rotational orientation on the dissociative chemisorption of CH<sub>4</sub> on Ni(111) transition metal surface", in preparation

## State-to-state study using the reactant coordinate based (RCB) method in the initial state-selected wave packet (ISSWP) framework
1. State-to-state reaction probabilities of the abstraction and exchange channels in the H + H<sub>2</sub>O reaction
   - Bin Zhao, Zhigang Sun, and Hua Guo, “A reactant-coordinate-based wave packet method for full-dimensional state-to-state quantum dynamics of tetra-atomic reactions: Application to both the abstraction and exchange channels in the H + H<sub>2</sub>O reaction”, [_J. Chem. Phys._ __144__, 064104 (2016)](https://aip.scitation.org/doi/10.1063/1.4941671).
2. State-to-state reaction probabilities of the two product channels in the HD + OH reaction
   - Bin Zhao, Zhigang Sun, and Hua Guo, “State-to-state mode selectivity in the HD + OH reaction: Perspectives from two product channels”, [_J. Chem. Phys._ __144__, 214303 (2016)](https://aip.scitation.org/doi/10.1063/1.4952764).
3. Differential cross sections of the H<sub>2</sub> + OH reaction
   - Bin Zhao, Zhigang Sun, and Hua Guo, “A reactant-coordinate-based approach to state-to-state differential cross sections for tetratomic reactions”,  [_J. Chem. Phys._ __145__, 184106 (2016)](https://aip.scitation.org/doi/10.1063/1.4966966).
4. Differential cross sections of the H + HOD reaction
   - Bin Zhao, Zhigang Sun, and Hua Guo, “State-to-state mode specificity in H + DOH(vOH=1) → HD + OH(v2=0) reaction: vibrational non-adiabaticity or local-mode excitation?”, [_Phys. Chem. Chem. Phys._ __20__, 191 (2018)](https://pubs.rsc.org/en/content/articlelanding/2018/cp/c7cp07199j#!divAbstract).
   - Bin Zhao, Uwe Manthe, and Hua Guo, “(2018 PCCP Hot Article) Fermi Resonance Controlled Product Branching in the H+HOD reaction”, [_Phys. Chem. Chem. Phys._ __20__, 17029 (2018)](https://pubs.rsc.org/en/content/articlelanding/2018/cp/c8cp02279h#!divAbstract).

## State-to-state study using the reactant-product decoupling (RPD) method in the initial state-selected wave packet (ISSWP) framework
1. Differential cross section of the D<sub>2</sub> + OH &#8594; D + DOH reaction
   - Bin Zhao, Zhigang Sun, and Hua Guo, “State-to-state differential cross sections for D<sub>2</sub> + OH &#8594; D + DOH reaction: Influence of vibrational excitation of OH reactant”, [_J. Chem. Phys._ __145__, 134308 (2016)](https://aip.scitation.org/doi/10.1063/1.4964322).
2. Differential cross section of the H + CD<sub>4</sub> &#8594; HD + CD<sub>3</sub> reaction
   - Bin Zhao, “Product state-pair correlated differential cross sections of the H + CD<sub>4</sub> &#8594; HD + CD<sub>3</sub> reaction at a collision energy of 0.72 eV”, in preparation

---

# Program language, parallelism, and external library dependence
The present version of the ABR program is written in Fortran.

Hybrid [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface) and [Open Multi-Processing (OpenMP)](https://en.wikipedia.org/wiki/OpenMP) are used to facilitate efficient distributed-memory and shared-memory multiprocessing programming for different computer systems and different molecular systems.  

The ABR program uses following packapges:
- the [ARPACK package](https://www.caam.rice.edu/software/ARPACK/) for efficient calculations of selected eigenstates.
- the [Linear Algebra Package (LAPACK)](http://www.netlib.org/lapack/) and [Basic Linear Algebra Subprograms (BLAS)](http://www.netlib.org/blas/) for efficient linear algebra operations.

---

# Parallel efficiency
## Scaling with the number of computing nodes
The current ABR program relies on the efficient transformation between the FBR and DVR, which allows an efficient application of operators on the wave function. The current MPI implemenation for distributed-memory programming requires a complete redistribution of the wave function among the computation nodes at each propagation step. The overhead of MPI communication is minimized by overlapping computation with communication so that the ABR program show a favorable quasi-linear scaling with respect to the number of computing nodes.

Some test calculations were performed in 2019 for the F + CHD<sub>3</sub> &#8594; HF + CD<sub>3</sub> reaction on the OCuLUS cluster at [the Paderborn Center for Parallel Computing (PC²)](https://pc2.uni-paderborn.de/about-pc2). The details of the nodes are listed in Table 1.

Table 1. Details of the compute nodes and the mode of parallelization.

| Location |  Paderborn               |
|--------------|-----------------|
|CPUs per node | 16        |
| CPU type     | Intel Xeon E5-2670, 2.6GHz |
|main memory per node	| 64 GB |
|interconnect|	QDR InfiniBand PCIe3, 40 Gbit/s Mellanox|
|accelerators used (such as GPUs)	| no |
|number of MPI processes per node	|16 |
|number of OpenMP threads per MPI process 	1|

The test calculations were performed using all the 16 cores on the computing nodes. The scaling is listed in Table 2 as the elapsed runtime (wall-clock time) T with respect to the number of nodes N. The resource usage evaluated as the product of runtime and number of nodes are also listed. In tests 1, 3, and 4, the runtime was systematically reduced when more computation nodes were used. The test calculations generally required similar amount of computation resources at a given problem size. It should be noted that the ABR program requires extensive MPI communication during the calculation. Latency of MPI communication is affected by connections among the computing nodes. Parallel efficiency depends on whether the nodes are connected to a common InfiniBand switch or not. The nodes in tests 1, 3 and 4 were randomly assigned by OpenCCS resource manager and they were not connected to a common InfiniBand switch. In test 2, the nodes were required to be connected to a common switch. Comparing tests 1 and 2 that both use 10 computations nodes, it is clear that test 2 required a short runtime due to less MPI latency in the communication. It should also be noted that the chunk size of data in each MPI communication affects the efficiency of the overlapped communication and computation. A optimized chunk size depends on the molecular system size and underlying hardware.

Table 2. Scaling of the elapsed runtime (wall-clock time) T with respect to the number of nodes N.  

|	|N	|T(N)	|N*T(N)	|comment|
|--------------|-----------------|--------------|-----------------|----|
|1	|10	|331	|3310	|Connected to different switches|
|2	|10	|289	|2890	|Connected to common switch|
|3	|20	|227	|4540	|Connected to different switches|
|4	|30	|151	|4530	|Connected to different switches|



## Massive parallel computation with a large amount of RAM
In the study of the non-adiabatic effects in the F(<sup>2</sup>P) + CHD<sub>3</sub> &#8594; HF + CD<sub>3</sub> reaction, a large amount of RAM was used to represent the wave functions on six coupled potential energy surfaces and to store the energy resolved scattering wave functions on a ultra-fine energy grid to resolve the sharp resonance features. This study used 50 computing nodes with 2 Intel® Xeon® Gold 6148 Processor each and 10 TB RAM in total.

## Parallel file system
The ABR program stores values of the potential energy operator on the disk so that they can be frequently used in the follow-up calculations. For polyatomic molecular systems, the direct-product type grids are large, and a parallel file system is helpful for efficient writing/reading of the potential energy data.

---

# Program installation and execution
## Package dependence
The compilation and installation of the ABR program require a Fortran compiler (such as the [intel Fortran compiler](https://www.open-mpi.org/), the [GNU Fortran compiler](https://gcc.gnu.org/fortran/), etc.), a MPI implementation (such as [MPICH](https://www.mpich.org/), [OpenMPI](https://www.open-mpi.org/), etc.), [GNU make](https://www.gnu.org/software/make/), the [ARPACK package](https://www.caam.rice.edu/software/ARPACK/), and the [Linear Algebra Package (LAPACK)](http://www.netlib.org/lapack/) and [Basic Linear Algebra Subprograms (BLAS)](http://www.netlib.org/blas/).

In the src folder exists an example makefile that has been used to compile the ABR program on two computing clusters: the [Center for Advanced Research Computing](https://carc.unm.edu/) at the University of New Mexico, New Mexico, USA and the [Paderborn Center for Parallel Computing (PC<sup>2</sup>)](https://pc2.uni-paderborn.de/), North Rhine-Westphalia, Germany. Following compiler and packages have been used:
- MPI implementation: OpenMPI/3.1.1
- Fortran compiler: Intel Fortran compiler in the Intel compiler and libraries Composer
- LAPACK and BLAS: Intel MKL library
- ARPACK: stand alone ARPACK library

Other installation options have not been fully tested. Please contact us in case of installation problem.

## Installation procedure
- Install Fortran compiler and MPI implementation. Update the FC variable in the make file according to your specific selection.
- Install LAPACK and BLAS packages. If intel MKL is used, update the MKL_PATH variable and update the library linking flag LDFLAGS variable in the make file.
- Install ARPACK package. Link the ARPACK package.
- Supply subroutines for evaluating potential energy on the grid. Update PES and PESPATH for the name and path of the PES subroutine. The LINKPES variable in the make file is for additional linking option of the PES subroutine.
- Finally, the 'make' command does the trick to compile the program.

## Execution procedure
- First of all, if the PES subroutine requires external data, please make sure the program can access these data in the execution folder of the program.
- A job submission script exists in the example folder. For example, four steps are required in the calculations of the non-adiabatic quenching of OH(A) by H<sub>2</sub>:
   1. Calculate all the potential energy values and the diabatic-adiabatic transformation matrices on all the discrete grids and store them on the disk.
   2. Calculate the potential energy values and the diabatic-adiabatic transformation matrices on an asymptotic analyzing plane.
   3. Optimizing the potential energy values and decompose them to different computing nodes.
   4. Perform the real-time propagation of the wave packet.
- Analyze the DATA that are generated and written on the disk.

---

# Code contribution
ABR program suite is developed and implemented by Bin Zhao during his post-doc with Prof. Hua Guo in the University of New Mexico, Albuquerque (A) and with Prof. Uwe Manthe in University Bielefeld, Bielefeld (B).
It is based on Prof. Minghui Yang's code used for the calculations of
total reaction probabilities of
the F + CHD<sub>3</sub> &#8594; HF + CD<sub>3</sub> reaction
for vanishing total system angular momentum (J<sub>tot</sub>=0).

---

# License
ABR program is an open source program that uses the GNU GPL v3 license. Please find the included license file for details.

---

# Contact
The ABR program is maintained by Bin Zhao.
For report of bug or supports on the use of the ABR program, please contact Bin Zhao.

---

# Acknowledgement
Bin Zhao sincerely acknowledges the indispensable discussions and supports
from Prof. Uwe Manthe, Prof. Hua Guo, and Prof. Zhigang Sun during the development of the ABR program.

