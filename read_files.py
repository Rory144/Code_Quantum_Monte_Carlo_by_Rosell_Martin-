# Import neccesary modules and libraries
import numpy as np 
import pandas as pd

# Function to read coordinates from an xyz file
def read_coordinates(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        
        # Extract atomic number (Z) from the second line        
        Z = [int(coord) for coord in lines[1].split()]
        # Extract number of electrons (ne) from the third line
        ne = int(lines[2].strip())
        # Read the coordinates from the file skipping the first 3 lines and using whitespace as delimiter
        R = pd.read_csv(filename, skiprows = 3, index_col=False, 
                       delimiter=r"\s+", names=['x', 'y', 'z'])

        # Change coordinates from angstrom to atomic units
        conversion = 1.8897259886
        R[['x','y','z']] *= conversion
        
        # Print parameters
        print("="*50)
        print(" Data of xyz file:")
        print("="*50)
        print("Atomic number, (Z):", Z)
        print("Number of electrons, (ne):", ne )
        print("Nucleus coordinates, (R):", R )

    return Z, ne, R

# Main function for testing the read_coordinates function
if __name__ == "__main__":
    filename = "H2+.xyz"
    read_coordinates(filename)

# Function to read data from a text file
def read_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        
        # Extract calculation type from the first line
        calculation_type = lines[0].split('#')[0].strip()
        
        # Check if the calculation type is valid and what type is 
        if calculation_type == "VMC":
            calculation_type == 'VMC'
        
        elif calculation_type == "PDMC":
            calculation_type = 'PDMC'
         
        else:
            print("Invalid calculation type.")
            exit() 
        
        # Extract the rest of parameters from subsequent lines
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

# Main function for testing the read_data function
if __name__ == "__main__":
    filename = "data.txt"
    read_data(filename)

