Quantum Monte Carlo Homework 

# Project Overview
The objective of the project is to develop a Quantum Monte Carlo (QMC) code capable of calculating the ground state energies of various molecular systems, including single and multi-electron configurations. This program will implement the Variational Monte Carlo (VMC) and Pure Diffusion Monte Carlo (PDMC) algorithms, allowing accurate estimation of the ground state energies. The key feature of the program is its adaptability to different molecular systems, such as the hydrogen atom, the helium atom, diatomic molecules such as H2 and H2+, and triatomic molecules such as H3+, with flexible input parameters including nuclear coordinates, wave function parameters, and simulation parameters. By reading input data from a configuration file, the program will seamlessly switch between VMC and PDMC modes and adjust calculations according to the specified molecular system. 

# Installation
1. Clone this repository: git clone https://github.com/Rory144/Code_Quantum_Monte_Carlo_by_Rosell_Martin-.git 
2. Navigate to the project directory: cd Code_Quantum_Monte_Carlo_by_Rosell_Martin

# Usage

This project presents different codes: 

1. data.txt ->  This file contains all the input parameters required to specify before executing the program, called by the main program.  

2. xyz_atomic.txt -> This file contains the nucleus coordinates and the atomic numbers required to specify that will be called by the main program. 

3. xyz_angstroms.txt -> This file is the same file that the previous one but the coordinates are shown in angstroms. 

4. conversion_units.py -> This program allows to convert angstroms coordinates to atomic unit, using as input file the xyz_angstroms and generating an output file named xyz_atomic.txt 
 
5. main_program.py -> The main program is in charge of reading the information in file.txt files (previous) and importing them to the other programs with the help of the coordinates.py and parameters.py modules. Also depending on the reading of the type of calculation, it executes the Variational or the Pure diffusion program of Monte Carlo. 

6. Local_energy.py -> This program compute the values of the potential, atomic orbitals, wavefunction, kinetic energy and finally, the local energy that is neccessary to the following programs: 

7. Variational_Monte_Carlo.py -> This program have all the calculations to compute the energy of the Variational Monte Carlo method.

8. Pure_Difussion_Monte_Carlo.py -> This program have all the calculations to compute the energy of the Pure Diffussion Monte Carlo method. 

# Execution 

# License
This project is licensed under the MIT License - see the LICENSE.md file for details.

