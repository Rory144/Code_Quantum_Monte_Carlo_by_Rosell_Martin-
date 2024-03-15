# ------------------------------------------------------------------------
#                         Main_program.py 
#-------------------------------------------------------------------------
#This program involving reading data from two different files:
#(xyz_atomic.txt and data.txt) 
#and performing calculations based on the parameters extracted from data.txt
# ------------------------------------------------------------------------

# Import neccesary modules and libraries
import parameters
import coordinates 
import numpy as np 
import pandas as pd
from MC import * 

# Function to read the coordinates from xyz_atomic.txt file
def read_coordinates(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

        # Split the first line into two parts and convert them to floats
        coordinates.Z = [float(coord) for coord in lines[1].split()]
        coordinates.ne = [float(coord) for coord in lines[2].split()]
        
        # Read the atoms & xyz coordinates
        coordinates.nxyz = pd.read_csv(filename, skiprows=3, index_col=False,
                       delimiter=r"\s+", names=['atom', 'x', 'y', 'z'])

        # Change coordinates from angstrom to atomic units
        conversion_units = 1.8897259886
        coordinates.nxyz[['x','y','z']] *= conversion_units
        
        # Print parameters for verification
        print("="*50)
        print(" Data of xyz file:")
        print("="*50)
        print("Number atomic, (Z):", coordinates.Z)
        print("Number of electrons, (ne):", coordinates.ne )
        print("Nucleus coordinates, (R):", coordinates.nxyz )

if __name__ == "__main__":
    filename = "H3+.xyz"
    read_coordinates(filename)

# Function to read parameters from data.txt file
def read_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

        parameters.calculation_type = lines[0].split('#')[0].strip()
        parameters.a = float(lines[1].split('#')[0].strip())
        parameters.N_MC = int(lines[2].split('#')[0].strip())
        parameters.dt = float(lines[4].split('#')[0].strip())
        parameters.tau = float(lines[5].split('#')[0].strip())
        parameters.Eref = float(lines[6].split('#')[0].strip())
        

        # Print parameters 
        print("="*50)
        print("Parameters:")
        print("="*50)
        print("Calculation Type:", parameters.calculation_type)
        print("Amplitude, (a):", parameters.a)
        print("Number of MC steps, (N_MC):", parameters.N_MC)
        print("Value of dt:", parameters.dt)
        print("Value of tau:", parameters.tau)
        print("Value of reference energy:", parameters.Eref)

if __name__ == "__main__":
    filename = "data.txt"
    read_data(filename)

    if parameters.calculation_type == "VMC":
        from Variational_Monte_Carlo import VMC
        VMC(parameters.a, parameters.N_MC)
    elif parameters.calculation_type == "PDMC":
        from Pure_Difussion_Monte_Carlo import PDMC
        PDMC(parameters.a, parameters.N_MC, parameters.dt, parameters.tau, parameters.Eref)
    else:
        print("Invalid calculation type.")
