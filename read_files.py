# ------------------------------------------------------------------------
#                         Main_program.py 
#-------------------------------------------------------------------------
#This program involving reading data from two different files:
#(xyz_atomic.txt and data.txt) 
#and performing calculations based on the parameters extracted from data.txt
# ------------------------------------------------------------------------

# Import neccesary modules and libraries
import numpy as np 
import pandas as pd

def read_coordinates(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

        Z = int(lines[1].strip())
        ne = int(lines[2].strip())
#        R_angs = lines[3].strip().split()
 #       R = [1.8897259886 * float(coor) for coor in R_angs]
        
        R = pd.read_csv(filename, skiprows=3, index_col=False, 
                       delimiter=r"\s+", names=['x', 'y', 'z'])

        # Change coordinates from angstrom to atomic units
        conversion = 1.8897259886
        R[['x','y','z']] *= conversion
        
        # Print parameters for verification
        print("="*50)
        print(" Data of xyz file:")
        print("="*50)
        print("Atomic number, (Z):", Z)
        print("Number of electrons, (ne):", ne )
        print("Nucleus coordinates, (R):", R )

    return Z, ne, R

if __name__ == "__main__":
    filename = "H3+.xyz"
    read_coordinates(filename)

def read_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

        calculation_type = lines[0].split('#')[0].strip()
        
        if calculation_type == "VMC":
            calculation_type == 'VMC'
        
        elif calculation_type == "PDMC":
            calculation_type = 'PDMC'
         
        else:
            print("Invalid calculation type.")
            exit() 

        a = float(lines[1].split('#')[0].strip())
        N_MC = int(lines[2].split('#')[0].strip())
        dt = float(lines[3].split('#')[0].strip())
        tau = float(lines[4].split('#')[0].strip())
        Eref = float(lines[5].split('#')[0].strip())
        

        # Print parameters 
        print("="*50)
        print("Parameters:")
        print("="*50)
        print("Calculation Type:", calculation_type)
        print("Amplitude, (a):", a)
        print("Number of MC steps, (N_MC):", N_MC)
        print("Value of dt:", dt)
        print("Value of tau:", tau)
        print("Value of reference energy:", Eref)
        print("="*50)
    
    return calculation_type, a, N_MC, dt, tau, Eref
if __name__ == "__main__":
    filename = "data.txt"
    read_data(filename)

