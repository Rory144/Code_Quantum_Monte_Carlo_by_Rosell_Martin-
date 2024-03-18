# Importing necessary libraries and modules 
import numpy as np


# Calculation of atomic orbital (phi) function
def phi(a, r, R):

    phi = 0.
    # Calculate the distance between electron and nucleus
    r_iA2 = (r[0] - R[0])**2 + (r[1] - R[1])**2 + (r[2] - R[2])**2
    if r_iA2 < 0:
        return 0 
    else: 
        r_iA = np.sqrt(r_iA2) 
    # Calculate phi 
    phi += (a**3/np.pi)**0.5 * np.exp(-a * r_iA)

    return phi 

# Calculation of wavefunction (psi) function
def psi(a, r, R):
    
    N = int(len(r)/3)  
    M = int(len(R)/3)  

    #Inizialize psi and the sum  
    psi = 1. 
    
    # Loop over electron indices
    for i in range(N):
        r_i = r[3*i : 3*i + 3 ]
        ii = 3 * i
        sum = 0. 
        # Loop over nuclei indices
        for A in range(M):
            R_A = R[3*A : 3*A + 3] 
            sum += phi(a, r_i, R_A) 
        psi *= sum
    return psi

# Calculation of the electron-electron potential (Vee). 
def Vee(r):

    # Number of electrons and nuclei
    N =int(len(r)/3)
    
    # Inizialize the potential
    Vee = 0.0

    # Outer loop for electron indices (i)
    for i in range(0, N):
        ii = 3 * i
        # Inner loop for electron indices (j)
        for j in range(i + 1, N):
            jj = 3 * j

            # Calculate the distance between electrons
            r_ij = np.sqrt((r[jj] - r[ii]) ** 2 + (r[jj + 1] - r[ii + 1]) ** 2 + (r[jj + 2] - r[ii + 2]) ** 2)
            
            # Avoid division by zero
            if r_ij == 0.:
                print("The potential at r=0 diverges")
                return float("inf")
            else:
                Vee += 1.0 / r_ij

    return Vee
    
# Calculation of the nuclei_nuclei potential (Vnn)
def Vnn(Z, R): 
    
    #Number of nucleus
    M =int(len(R)/3)
    
    #Initialize the potential
    Vnn = 0.0
    
    # For one nucleus, the potential is null, so 
    if M == 1:
        return 0.
    # For more than 1 nucleus: 
    else:
    
        #Outer loop for nuclei indices (A)
        for A in range(M + 1):
            AA = 3 * A
            R_A = R[AA:AA+3]
            # Inner loop for nuclei indices (B)
            for B in range(A + 1, M):
                BB = 3 * B
                R_B = R[BB:BB+3]
                
                # Calculate the distance between nuclei
                R_AB = np.sqrt((R_A[0] - R_B[0])**2 + (R_A[1] - R_B[1])**2 + (R_A[2] - R_B[2])**2)       
         
                # Avoid division by zero
                if R_AB == 0:
                    print("potential at r=0 diverges")
                    return float("inf")
                else: 
                    Vnn += Z[A]* Z[B]/ R_AB
        return Vnn

# Calculation of the nucleus-electron potential (Ven)
def Ven(Z, r, R): 
   
    #Number of electrons and nucleus
    N =int(len(r)/3)
    M =int(len(R)/3)
    
    if len(Z) != m:
        raise ValueError('The number of Z values must be the same as the number of nucleus')

    # Initialize the potential
    Ven = 0.0

    #Outer loop for nuclei indices (A)
    for i in range(N):
        ii = 3 * i
        r_i = r[ii:ii+3]

        #Inner loop for nuclei indices (B)
        for A in range(M):
            AA = 3 * A
            R_A = R[AA:AA+3]
                
            # Calculate the distance between nuclei
            r_iA = np.sqrt((r_i[1] - R_A[1])**2 + (r_i[2] - R_A[2])**2 + (r_i[3] - R_A[3])**2)               
            # Avoid division by zero
            if r_iA == 0:
               print("potential at r=0 diverges")
               return -float("inf")
            else: 
                Ven += float(Z[A])/r_iA

        return -Ven

# Calculation of the total potential, (V_total)
def V_total(Z, r, R): 
    
    #Number of electrons and nucleus
    N =int(len(r)/3)
    M =int(len(R)/3)
    
    # Initialize the potential
    V_total = 0.0
    
    # Sum up the contributions from electron-electron, nuclei-nuclei, and electron-nuclei potentials
    V_total += Vee(r) + Vnn(Z, R) + Ven(Z, r, R)

    return V_total 

# Calculation of kinetic energy function
def kinetic(a, r, R):
    
    #Number of electrons and nucleus
    N =int(len(r)/3)
    M =int(len(R)/3)
    
    sum = 0. 

    # Loop over electron indices
    for i in range (N):
        ii = 3*i
        r_i = r[ii:ii+3]

        numerator = 0. 
        denominator = 0. 

        # Loop over nuclei indices    
        for A in range(M):
            AA = 3*A
            R_A = R[AA:AA+3]

            # Calculate the distance between electron and nucleus
            r_iA = np.sqrt((r[0] - R[0])**2 + (r[1] - R[1])**2 + (r[2] - R[2])**2)
            
            if r_iA == 0.:
                return float("inf")
            else:   
            
                D_iA = (3. / np.abs(r_iA)) - (1. + a * np.abs(r_iA)) / np.abs(r_iA)
                numerator *= - a * D_iA * phi(a, r_i, R_A)

        sum += (numerator / denominator) 
    kinetic = -0.5 * sum 

    return kinetic 

# Calculation of local energy function (e_loc)
def e_loc(a, r, R, Z):
    
    #Inizialize the local energy 
    e_loc = 0

    # Sum up kinetic energy and total potential energy
    e_loc += kinetic(a, r, R) + V_total(Z, r, R)
    
    return e_loc

