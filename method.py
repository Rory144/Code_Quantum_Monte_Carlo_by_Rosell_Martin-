
# Importing necessary libraries and modules
from math import sqrt
import numpy as np 
from Local_energy import phi, psi, e_loc

def ave_error(arr):

    M = len(arr)
    assert(M > 0)

    # "M should be greater than 0 because is in the denominator"
    if M == 1:
        average = arr[0]
        error   = 0.0
    else:
        average = sum(arr)/M
        variance = 1./(M-1) * sum( [(x - average) **2 for x in arr ] )
         # To correct for bias
        error = sqrt(variance/M)

    return (average, error)

# Write a function to compute the dvec vector:
def drift(a, r, R):
    # Number of electrons and nuclei
    M = int(len(str(R))/3)
    N = int(len(r)/3)

    drift = np.zeros(len(r))

    for i in range(N):
        ii = 3 * i
        r_i = r[ii:ii + 3]
        numerator = 0.
        denominator = 0.

        for A in range(M):
            AA = 3 * A
            R_A = R[AA:AA + 3]

            r_iA2 = ((r_i[0] - R_A[0])**2 +
                     (r_i[1] - R_A[1])**2 +
                     (r_i[2] - R_A[2])**2)

            # Evitar división por cero
            if r_iA2 != 0:
                r_iA = np.sqrt(r_iA2)
                C_iA = np.sqrt((r_i[0] - R_A[0])**2 +
                               (r_i[1] - R_A[1])**2 +
                               (r_i[2] - R_A[2])**2)

                numerator += C_iA * phi(a, r_i, R_A)
                denominator += phi(a, r_i, R_A)

        # Multiplicar después de acumular
        drift[ii:ii+3] = -a * (numerator / denominator)

    return drift

def VMC(a, ne, R, Z, N_MC, dt):
    
    e_loc  = 0.
    accep = 0

    r_old   = np.random.normal(loc = 0., scale = 1.0, size = 3 * ne)
    d_old   = drift(a, r_old, R)
    d2_old  = np.dot(d_old, d_old)
    psi_old = psi(a, r_old, R)

    for istep in range(N_MC):
        chi = np.random.normal(loc=0., scale=1.0, size= 3 * ne)
        ## The Gaussian random number with zero mean and variance dt is -> chi

        #Evaluation of the local energy at r_old.
        e_loc += e_loc(a, r_old, R, Z)

        # Compute a new position:
        r_new = r_old + dt * d_old + np.sqrt(dt) * chi

        # Evaluate psi and dvec vector at the new position
        d_new   = drift(a, r_new, R)
        d2_new  = np.dot(d_new, d_new)
        psi_new = psi(a, r_new, R)

        # Compute the ratio for Metropolis
        first_part  = np.dot((d_new + d_old), (r_new - r_old))
        second_part = 0.5 * dt * (d2_new - d2_old)
        exponent = first_part + second_part 

        ratio_A = min(1., (psi_new/psi_old)**2 * np-exp(-exponent)) 

        # Draw a uniform random number, v,  between 0 to 1:
        v = np.random.uniform() 

        # If u is less or equal than the ratio_A -> accept the move: set r{n+1} = r_new

        if v <= ratio_A :
            accep += 1

            r_n1  = r_new
            d_n1   = d_new
            d2_n1  = d2_new
            psi_n1 = psi_new

        # else, reject the move: set r{n+1} = r_old 
        else: 
            r_n1  = r_old
            d_n1   = d_old
            d2_n1  = d2_old
            psi_n1 = psi_old
      
        # Update the position   
        r_old = r_n1
        d_old = d_n1
        d2_old = d2_n1
        psi_old = psi_n1
    
    return e_loc/N_MC, accep/N_MC

def PDMC(a, ne, Z, R, N_MC, dt, tau, Eref):

    e_loc = 0.
    accep = 0
    norma = 0.

    #Starting with w(r0) in 1 and tau in 0
    w           = 1.
    tau_current = 0.

    # Evaluation the local energy at r_old
    r_old   = np.random.normal(loc = 0., scale = 1.0, size = 3 * ne)
    d_old   = drift(a, r_old, R)
    d2_old  = np. dot(d_old, d_old) 
    psi_old = psi(a, r_old, R)

    for istep in range(N_MC):
        #Evaluation of the local energy at r_old 
        e_loc = e_loc(a, r_old, R, Z)
        
        # Compute the contribution to the weight and update it at r_old
        w *= np.exp(-dt*(e_loc - Eref))

        # Accumulate the weight for the normalization factor and the energy
        norma += w
        energy += w * e_loc

        # Update tau
        tau_current += dt

        # Reset when the long projection time is reached
        if tau_current > tau:
            w = 1.
            tau_current = 0.

        #Compute chi as a random number
        chi = np.random.normal(loc = 0., scale = 1.0, size = 3 * ne)

        #Compute the new position
        r_new = r_old + dt * d_old + np.sqrt(dt) * chi

        # Evaluate the dvec vector and the psi at the new position 
        d_new = drift(a, r_new, R)
        d2_new = np.dot(d_new, d_new)
        psi_new = psi(a, r_new, R)

        # Compute the ratio
        first_part  = np.dot((d_new + d_old), (r_new - r_old))
        second_part = 0.5 * dt * (d2_new - d2_old)
        exponent = first_part + second_part 

        ratio_A = min(1., (psi_new/psi_old)**2 * np-exp(-exponent))  

        #Draw a unifrom random number, v 
        v = np.random.uniform()  
        
        # If u is less or equal than the ratio_A -> accept the move: set r{n+1} = r_new
        if v <= ratio_A:
            accep += 1
            r_n1    = r_new
            d_n1  = d_new
            d2_n1  = d2_new
            psi_n1 = psi_new

        # else, reject the move: set r{n+1} = r_old 
        else:         
            r_n1  = r_old
            d_n1   = d_old
            d2_n1  = d2_old
            psi_n1 = psi_old
        
        # Update the position   
        r_old = r_n1
        d_old = d_n1
        d2_old = d2_n1
        psi_old = psi_n1

    return energy/norma, accep/N_MC

#
