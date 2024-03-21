#import neccesary libraries and modules
from math import sqrt
import numpy as np 
from Local_energy import psi, phi, e_loc

# Function to calculate the average and error of an array
def ave_error(arr):

    # Determine the length of the array
    M = len(arr)
    assert(M>0)

    # If there's only one element in the array, return the element and error as zero
    if M == 1:
        average = arr[0]
        error   = 0.

    else:
        # Calculate the average of the array
        average = sum(arr)/M
        # Calculate the variance of the array
        variance = 1./(M-1) * sum( [ (x - average)**2 for x in arr ] )
        # Calculate the error from the variance
        error = sqrt(variance/M)

    # Return the calculated average and error as a tuple
    return (average, error)

# Function to compute the drift vector 
def drift(a, r, R):
    
     # Determine the number of electrons and nucleus 
     N = int(len(r)/3)
     M = int(len(R)/3) 

     # Create an array to store drift values
     drift = np.zeros(len(r))

     # Loop through each electron
     for i in range(N): 
         ii = 3 * i
         r_i = r[ii:ii+3]
        
         #Initialize each component of numerator and then the denominator
         numerator_x = 0. 
         numerator_y = 0. 
         numerator_z = 0. 
         denominator = 0. 
            
         # Loop through each nucleus
         for A in range(M): 
            AA = 3 * A
            R_A = R[AA:AA+3]

            # Calculate distance between electron i and nucleus A
            r_iA = np.sqrt((r_i[0] - R_A[0])**2 +
               (r_i[1] - R_A[1])**2 +
               (r_i[2] - R_A[2])**2)

            # Avoid division by zero
            if r_iA == 0:
                    return -float("inf") # Return negative infinity if distance is zero
            else:
                Cx_iA = (r_i[0]-R_A[0])/r_iA 

                Cy_iA = (r_i[1]-R_A[1])/r_iA 

                Cz_iA = (r_i[2]-R_A[2])/r_iA
                
                numerator_x += Cx_iA * phi(a, r_i, R_A) 
                numerator_y += Cy_iA * phi(a, r_i, R_A) 
                numerator_z += Cz_iA * phi(a, r_i, R_A) 
                
                denominator += phi(a, r_i, R_A)

         # Calculate drift components for electron (i)
         drift[ii] = -a * numerator_x/denominator 
         drift[ii+1] = -a * numerator_y/denominator 
         drift[ii+2] = -a * numerator_z/denominator
    
     return drift 

# Function to perform the Variational Monte Carlo (VMC) simulation
def VMC(a, N_MC, dt, ne, Z, R): 
    
    #Initialize the total energy and acceptance 
    energy = 0. 
    accept = 0.
    
    # Generate random initial positions for the particles
    ro = np.random.normal(loc=0., scale=1., size = 3 * ne)
    # Compute initial drift
    do = drift(a, ro, R)
    # Compute square of the initial drift
    d2o = np.dot(do, do)
    # Evaluate wavefunction at initial position
    psio = psi(a, ro, R)
    
    # Loop through Monte Carlo steps
    for istep in range(N_MC):
        
        chi = np.random.normal(loc=0., scale=1., size=3*ne)

        # Evaluate the local energy at r_old
        energy += e_loc(a, ro, R, Z)

        # Compute new position
        rnew = ro + dt * do + np.sqrt(dt) * chi

        # Compute new drift
        dnew = drift(a, rnew, R)
        d2new = np.dot(dnew, dnew)

        # Evaluate the wavefunction at new position
        psinew = psi(a, rnew, R)
        
        # Compute the Metropolis algorithm 
        exponent = np.dot((rnew - ro), (dnew + do)) + 1./2. * dt * (d2new - d2o)

        # Compute the acceptance ratio 
        ratio_A = min(1., ((psinew / psio)**2 * np.exp(-exponent)))

        # Draw an uniform random number
        v = np.random.uniform()

        # If v <= ratio ->  accept the move: set r_{n+1} = rnew
        if v <= ratio_A:
            accept += 1
            r_n1 = rnew
            d_n1 = dnew
            d2_n1 = d2new
            psi_n1 = psinew

        # else, reject the move: set r_{n+1} = ro
        else:
            r_n1 = ro
            d_n1 = do
            d2_n1 = d2o
            psi_n1 = psio
        
        # Update position, drift, wavefunction, and square of drift
        ro = r_n1
        do = d_n1
        d2o = d2_n1
        psio = psi_n1

    return energy/N_MC, accept/N_MC

# Function to perform the Pure Difussion Monte Carlo (PDMC) simulation
def PDMC(a, N_MC, dt, tau, Eref, ne, R, Z):
    
    #Initialize total energy, the acceptance and the normalization factor
    energy = 0. 
    accept = 0
    norma = 0.

    # Start with W(r0) = 1, tau_current = 0
    w = 1.
    tau_current = 0.

    # Generate random position and a  drift vector around ro
    ro = np.random.normal(loc=0., scale=1., size=3*ne)
    do = drift(a, ro, R)
    # Compute square of initial drift
    d2o = np.dot(do, do)
    # Evaluate wavefunction at initial position
    psio = psi(a, ro, R)

    # Loop through Monte Carlo steps
    for istep in range(N_MC):
        
        # Evaluate the local energy at rn
        local_e = e_loc(a, ro, R, Z)

        # Compute the contribution to the weight & update it
        w *= np.exp(-dt * (local_e - Eref))
        # Add up the weight for the normalization factor and the energy 
        norma += w
        energy += w * local_e

        # Update tau
        tau_current += dt

        # Reset when the long projection time has been reached
        if tau_current > tau:
            w = 1.0
            tau_current = 0. 

        chi = np.random.normal(loc=0., scale=1., size=3*ne)

        # Compute new position
        rnew = ro + dt * do + np.sqrt(dt) * chi

        # Compute new drift and the psi at the new position 
        dnew = drift(a, rnew, R)
        d2new = np.dot(dnew, dnew)
        psinew = psi(a, rnew, R)
        
        # Metropolis algorithm 
        first_part = np.dot((dnew+do), (rnew-ro))
        second_part = 1./2. * (d2new - d2o) *dt 
        exponent = first_part + second_part 

        # Compute the ratio
        ratio_A = min(1., ((psinew / psio)**2 * np.exp(-exponent)))

        # Draw an uniform random number
        v = np.random.uniform()

        # If v <= ratio ->  accept the move: set r_{n+1} = rnew
        if v <= ratio_A:
            accept += 1
            r_n1 = rnew
            d_n1 = dnew
            d2_n1 = d2new
            psi_n1 = psinew

        # else, reject the move: set r_{n+1} = ro
        else:
            r_n1 = ro
            d_n1 = do
            d2_n1 = d2o
            psi_n1 = psio
        
        # Update position, drift, wavefunction, and square of drift
        ro = r_n1
        do = d_n1
        d2o = d2_n1
        psio = psi_n1

    return energy/norma, accept/N_MC
