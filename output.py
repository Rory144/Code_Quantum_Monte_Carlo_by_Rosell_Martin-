# Import necessary modules and libraries
from read_files import read_coordinates, read_data
from Local_energy import psi, e_loc
from MC import * 

# Read coordinates from external file
filename_xyz = "H3+.xyz"
Z, ne, R = read_coordinates(filename_xyz)

# Read data from external file
filename_data = "data.txt"
calculation_type, a, N_MC, dt, tau, Eref = read_data(filename_data) 

##----------------------------------------------------

X0 = [VMC(a, ne, Z, R, N_MC, dt) for i in range(30)]

# Energy
X = [ x for (x, _) in X0 ]
E, deltaE = ave_error(X)

# Acceptance rate
A, deltaA = ave_error(X)

print("="*50)
print("Generalized Metropolis Sampling")

#-------------------------------------------------------------------------

X0 = [ PDMC(a, ne, Z, R, N_MC, dt, tau, Eref) for i in range(30)]

print("="*50)
print("Pure dufussion Monte Carlo")
print("."*50)

#Energy
X = [ x for (x, _) in X0 ]
E, deltaE = ave_error(X)

print(f"The energy and the error are: {E} +/- {deltaE}")
print("."*50)

# Acceptance rate
A, deltaA = ave_error(X)
print(f"The aceptance rate and the error are: {A} +/- {deltaA}")


# Calculate VMC and PDMC methods
#if calculation_type == "VMC":
#    result_VMC = calculate_methods.calculate_VMC(potential, a, N_MC, dt, tau, Eref)
#    print("VMC Result:", result_VMC)
#
#elif calculation_type == "PDMC":
#    result_PDMC = calculate_methods.calculate_PDMC(potential, a, N_MC, dt, tau, Eref)
#    print("PDMC Result:", result_PDMC)
#
#else:
#    print("Invalid calculation type.")

