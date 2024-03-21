# Import necessary modules and libraries
import sys
from read_files import read_coordinates, read_data
from Local_energy import psi, e_loc
from method import * 

#Redirect standard output to output file
sys.stdout = open("H2+_VMC.out", "w")

# Read coordinates from external file
filename_xyz = "H2+.xyz"
Z, ne, R_nucleus = read_coordinates(filename_xyz)

R = []
for index, row in R_nucleus.iterrows():
  R.extend([row['x'], row['y'], row['z']])

# Read data from external file
filename_data = "data.txt"
calculation_type, a, N_MC, dt, tau, Eref = read_data(filename_data) 

# Calculate VMC and PDMC methods
if calculation_type == "VMC":
    print("The result obtained with the Variational Monte Carlo (VMC) method are: ")

    X0 = [VMC(a, N_MC, dt, ne, Z, R) for i in range(30)]

    # Energy
    X = [ x for (x, _) in X0 ]
    E, deltaE = ave_error(X)

    print(f"The energy and the error are: {E} +/- {deltaE}")
    print("."*50)

    # Acceptance rate
    X = [ x for (_, x) in X0 ]
    A, deltaA = ave_error(X)
    print(f"The aceptance rate and the error are: {A} +/- {deltaA}")

elif calculation_type == "PDMC":
    print(" The results obtained with the Pure dufussion Monte Carlo (PDMC) method")
    print("."*50)

    X0 = [ PDMC(a, N_MC, dt, tau, Eref, ne, R, Z) for i in range(30)]
    
    #Energy
    X = [ x for (x, _) in X0 ]
    E, deltaE = ave_error(X)

    print(f"The energy and the error are: {E} +/- {deltaE}")
    print("."*50)
    
    # Acceptance rate
    X = [ x for (_, x) in X0 ]
    A, deltaA = ave_error(X)
    print(f"The aceptance rate and the error are: {A} +/- {deltaA}")
    print("="*50)

else:
    print("Invalid calculation type.")

# Close the output file: 
sys.stdout.close()

