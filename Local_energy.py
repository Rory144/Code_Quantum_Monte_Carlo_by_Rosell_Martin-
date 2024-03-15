# Importing necessary libraries and modules 
import numpy as np

# Definition of parameters imported from other modules
#a = parameters.a   # Amplittude of the wavefunction 'a'
#Z = coordinates.Z  # Array containing atomic numbers of nuclei
#R = coordinates.R  # Array containing nuclei coordinates

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

# Print the atomic orbital
#print(f'The atomic orbital is:', phi(a, r, R))
#print("="*50) 

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

# Print the wavefunction
#print(f'The wavefunction is:', psi(a, r, R))
#print("="*50)

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
    
# Print the electron-electron potential
#print("=" * 50)
#print(f"The electron-electron potential is:", Vee(r) )

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
            # Inner loop for nuclei indices (B)
            for B in range(A + 1, M):
                BB = 3 * B
                
                # Calculate the distance between nuclei
                R_AB = np.sqrt((R[BB] - R[AA])**2 + (R[BB+1] - R[AA+1])**2 + (R[BB+2] - R[AA+2])**2)       
         
                # Avoid division by zero
                if R_AB == 0:
                    print("potential at r=0 diverges")
                    return float("inf")
                else: 
                    Vnn += Z[A]* Z[B]/ R_AB
        return Vnn

# Print the nuclei-nuclei potential
#print("-"*50) 
#print("The nuclei_nuclei potential is:", Vnn(Z, R))

# Calculation of the nucleus-electron potential (Ven)
def Ven(Z, r, R): 
   
    #Number of electrons and nucleus
    N =int(len(N)/3)
    M =int(len(R)/3)

    # Initialize the potential
    Ven = 0.0

    # Outer loop for electron indices (i)
    for i in range (N):
        ii = 3*i
        # Inner loop for nuclei indices (A)
        for A in range(M):
            AA = 3*A

            # Calculate the distance between electron and nucleus
            r_iA = np.sqrt((r[ii] - R[AA])**2 + (r[ii + 1] - R[AA + 1])**2 + (r[ii + 2] - R[AA + 2])**2)

            # Avoid division by zero
            if r_iA == 0:
               print("Potential at r=0 diverges")
               return -float("inf")
            else:
                Ven +=  Z[A] / r_iA
    return -Ven

#Print the nuclei-nuclei potential
#print("-"*50) 
#print("The electron_nuclei potential is:", Ven(Z, r, R))
#print("="*50) 

# Calculation of the total potential, (V_total)
def V_total(Z, r, R): 
    
    #Number of electrons and nucleus
    N =int(len(N)/3)
    M =int(len(R)/3)
    
    # Initialize the potential
    V_total = 0.0
    
    # Sum up the contributions from electron-electron, nuclei-nuclei, and electron-nuclei potentials
    V_total += Vee(r) + Vnn(Z, R) + Ven(Z, r, R)

    return V_total 

# Print the total potential
#print("The total potential is:", V_total(Z, r, R))
#print("="*50) 

# Calculation of kinetic energy function
def kinetic(a, r, R):
    
    #Number of electrons and nucleus
    N =int(len(N)/3)
    M =int(len(R)/3)
    
    sum = 0. 

    # Loop over electron indices
    for i in range (N):
        ii = 3*i
        
        numerator = 0. 
        denominator = 0. 

        # Loop over nuclei indices    
        for A in range(M):
            AA = 3*A
            
            # Calculate the distance between electron and nucleus
            r_iA = np.sqrt((r[ii] - R[AA])**2 + (r[ii + 1] - R[AA + 1])**2 + (r[ii + 2] - R[AA + 2])**2)
            
            if rij == 0.:
                return float("inf")
            else:   
            
                D_iA = (3. / np.abs(r_iA)) - (1. + a * np.abs(r_iA)) / np.abs(r_iA)
                numerator *= - a * D_ij * phi(a, r_i, R_A)

        suma += (numerator / denominator) 
    kinetic = -0.5 * suma 

    return kinetic 

# Print the kinetic energy
#print(f'The kinetic energy is:', kinetic(a, r, R))
#print("="*50) 

# Calculation of local energy function (e_loc)
def e_loc(a, r, R, Z):
    
    #Inizialize the local energy 
    e_loc = 0

    # Sum up kinetic energy and total potential energy
    e_loc += kinetic(a, r, R) + V_total(Z, r, R)
    
    return e_loc

# Print the local energy
#print(f'The local energy is:', e_loc(a, Z, r, R))
#print("="*50) 
