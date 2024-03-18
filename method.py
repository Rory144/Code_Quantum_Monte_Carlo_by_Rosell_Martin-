from math import sqrt
import numpy as np 
from Local_energy import psi, phi, e_loc

def ave_error(arr):
    M = len(arr)
    assert(M>0)

    if M == 1:
        average = arr[0]
        error   = 0.

    else:
        average = sum(arr)/M
        variance = 1./(M-1) * sum( [ (x - average)**2 for x in arr ] )
        error = sqrt(variance/M)

    return (average, error)

def drift(a, r, R):

     N = int(len(r)/3)
     M = int(len(R)/3) 

     drift = np.zeros(len(r))

     for i in range(N): 
         ii = 3 * i
         r_i = r[ii:ii+3]
         numerator = 0. 
         denominator = 0. 

         for A in range(M): 
            AA = 3 * A
            R_A = R[AA:AA+3]
            
            r_iA = np.sqrt((r_i[0] - R_A[0])**2 +
               (r_i[1] - R_A[1])**2 +
               (r_i[2] - R_A[2])**2)
            print('h',r_iA)
            if r_iA == 0:
                    return -float("inf")
            else:

                C_iA = ((r_i[0]-R_A[0])/r_iA +
                       (r_i[1]-R_A[1])/r_iA +
                       (r_i[2]-R_A[2])/r_iA 
                        )
            numerator = C_iA * phi(a, r_i, R_A) 
            denominator = phi(a, r_i, R_A)

            drift = -a * numerator/denominator 
    
            return drift 

def VMC(a, N_MC, dt, ne, Z, R): 
    
    energy = 0. 
    accept = 0.

    ro = np.random.normal(loc=0., scale=1., size = 3 * ne)
    do = drift(a, ro, R)
    d2o = np.dot(do, do)
    psio = psi(a, ro, R)
    
    for istep in range(N_MC):
        
        chi = np.random.normal(loc=0., scale=1., size=3*ne)

        # Evaluate the local energy at r_n
        energy += e_loc(a, ro, R, Z)

        # Compute new position
        rnew = ro + dt * do + np.sqrt(dt) * chi

        # Compute new drift
        dnew = drift(a, rnew, R)
        d2new = np.dot(dnew, dnew)

        # Evaluate Psi at the new position
        psinew = psi(a, rnew, R)
        
        exponent = np.dot((rnew - ro), (dnew + do)) + 1./2. * dt * (d2new - d2o)

        # Compute the ratio
        ratio_A = min(1., ((psinew / psio)**2 * np.exp(-exponent)))

        # Draw an uniform random number
        v = np.random.uniform()

        # If v <= ratio ->  accept the move: set r_{n+1} = rnew
        if v <= ratio_A:
            accept += 1
            r_n1 = rnew
            d_n1 = dpnew
            d2_n1 = d2new
            psi_n1 = psinew

        # else, reject the move: set r_{n+1} = ro
        else:
            r_n1 = ro
            d_n1 = do
            d2_n1 = d2o
            psi_n1 = psio

        # Update the position
        ro = r_n1
        do = d_n1
        d2o = d2_n1
        psio = psi_n1

        return e_local/N_MC, accept/N_MC

def PDMC(a, N_MC, dt, tau, Eref, ne, R, Z):
    
    energy = 0. 
    accept = 0
    norma = 0.

    # Start with W(r0) = 1, tau_current = 0
    w = 1.
    tau_current = 0.

    # Generate random position and a  drift vector around ro
    ro = np.random.normal(loc=0., scale=1., size=3*ne)
    do = drift(a, ro, R)
    psio = psi(a, ro, R)
    
    for istep in range(N_MC):
        # Evaluate the local energy at rn
        e_loc = e_loc(a, ro, R, Z)

        # Compute the contribution to the weight & update it
        w *= np.exp(-dt * (e_loc - Eref))

        # Add up the weight for the normalization factorand the energy 
        norma += w
        energy += w * e_loc

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
        psinew = psi(a, rnew, R)
       
        exponent = np.dot((rnew - ro), (dnew + do)) + 1./2. * dt * (d2new - d2o)

        # Compute the ratio
        ratio_A = min(1., ((psinew / psio)**2 * np.exp(-exponent)))

        # Draw an uniform random number
        v = np.random.uniform()

        # If v <= ratio ->  accept the move: set r_{n+1} = rnew
        if v <= ratio_A:
            accept += 1
            r_n1 = rnew
            d_n1 = dpnew
            d2_n1 = d2new
            psi_n1 = psinew

        # else, reject the move: set r_{n+1} = ro
        else:
            r_n1 = ro
            d_n1 = do
            d2_n1 = d2o
            psi_n1 = psio

        # Update the position
        ro = r_n1
        do = d_n1
        d2o = d2_n1
        psio = psi_n1

        return e_local/N_MC, accept/N_MC
