Quantum Monte Carlo Homework 

# Project Overview
The objective of the project is to develop a Quantum Monte Carlo (QMC) code capable of calculating the ground state energies of various molecular systems, including single and multi-electron configurations. This program will implement the Variational Monte Carlo (VMC) and Pure Diffusion Monte Carlo (PDMC) algorithms, allowing accurate estimation of the ground state energies. The key feature of the program is its adaptability to different molecular systems, such as the hydrogen atom, the helium atom, diatomic molecules such as H2 and H2+, and triatomic molecules such as H3+, with flexible input parameters including nuclear coordinates, wave function parameters, and simulation parameters. By reading input data from a configuration file, the program will seamlessly switch between VMC and PDMC modes and adjust calculations according to the specified molecular system. 

# Installation
1. Clone this repository: git clone https://github.com/Rory144/Code_Quantum_Monte_Carlo_by_Rosell_Martin-.git 
2. Navigate to the project directory: cd Code_Quantum_Monte_Carlo_by_Rosell_Martin

# Usage

This project presents different codes: 

1. data.txt ->  This file contains all the input parameters required to specify before executing the program, called by the main_program.  

2. Geomtry directory -> This directory contains the files showing the charge of each atom, the number of electrons and the geometry of each atom or molecule being analyzed. For execution the files must be in the folder where the main program and the read file are located. 

3. read_files.py -> This program is in charge of reading the geometry files and the necessary parameters of the execution. It also converts the units of the xyz file that are in angstroms to atomic units. 

4. Local_energy.py -> This program calculates the necessary functions to compute the local energy such as kinetic energy, potential, wave function and atomic orbitals. 

5. method.py -> This program is in charge of computing the Variational Monte Carlo(VMC) and Pure diffusion Monte Carlo(PDMC) methods.

6. main_program.py -> This program is in charge of importing the three programs mentioned above, importing from Local_energy.py the functions needed to compute the methods, and from read_file.py the functions needed to read the information. In addition, it distinguishes whether the calculation type is VMC or PDMC and executes one or the other. Finally, it writes the results to an external file.

# Execution 
To run this program you must execute the program main_program.py. To do so, type the following command: 
python3 main_program.py 

Note: Remember that both the data.txt file and the file.xyz files must be in the same directory as the programs required for execution. 
# License
This project is licensed under the MIT License - see the LICENSE.md file for details.

