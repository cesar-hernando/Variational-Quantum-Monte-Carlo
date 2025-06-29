"""
In this file, we define functions to perform the MCMC sampling and energy calculation.
"""

import numpy as np
import matplotlib.pyplot as plt

def trial_wavefunction(x, alpha, hamiltonian_label):
    """
    This function defines the trial wavefunction, which is to be optimized.
    
    Parameters
    ----------
    alpha: np.ndarray
        array of variational parameters
    x: np.ndarray 
        contains degrees of freedom 
    hamiltonian_label: str
        Name of system that is being minimized
    Returns
    -------
    psi: np.imag
        wavefunction evaluated at input positions
    """

    # For each Hamiltonian, we define a different trial wave function, which is taken from [Thijssen]

    if hamiltonian_label == 'QHO_1D':
        psi = np.exp(-alpha*x**2).astype(np.complex128)
    elif hamiltonian_label == 'Hydrogen':
        psi = np.exp(-alpha*np.linalg.norm(x)).astype(np.complex128)

    elif hamiltonian_label == "Helium":
        x1 = x[0]
        x2 = x[1]

        r1=np.linalg.norm(x1)
        r2=np.linalg.norm(x2)
        r12=np.linalg.norm(x1-x2)
        
        psi = np.exp(-2*r1-2*r2+r12/(2*(1+alpha*r12))).astype(np.complex128)

    return psi    

def probability_density(psi):
    """
    Parameters
    ----------
    psi: np.ndarray
        array of wavefunction evaluated at input positions
    
    Returns
    -------
    density: np.ndarray
        array of wavefunction evaluated at input positions
    """
    
    density = psi*np.conjugate(psi)
    return density

def metropolis(alpha, num_samples, burn_in_period, d, thinning_factor, hamiltonian_label):
    """
    This function performs the MCMC algorithm to sample from the probability distribution dictated by psi.
    
    Parameters
    ----------
    alpha: np.ndarray
        array of variational parameters
    num_samples: int
        Number of desired samples from distriubution
    burn_in_period: int
        Number of initial samples to be discarded from MCMC method  
    d: float
        Maximum displacement for trial move
    hamiltonian_label: str
        Name of system that is being minimized
        
    Returns
    -------
    positions_sample: np.ndarray
        position samples according to probability distribution
    """
    
    if hamiltonian_label == 'QHO_1D':
        positions_samples = np.zeros((1, num_samples))       # Array storing the samples of the positions following the desired probability distribution
        position_initial = np.array([0])                     # Initial state (position in the x axis)
        
    elif hamiltonian_label == 'Hydrogen':
        positions_samples = np.zeros((3, num_samples))       # Array storing the samples of the positions following the desired probability distribution
        position_initial = np.array([1,0,0])                 # Initial state (cannot initialize in origin as the hamiltonian explodes there)

    elif hamiltonian_label == "Helium":
        positions_samples = np.zeros((2, 3, num_samples))    #Array storing the samples of the positions following the desired probability distribution
        position_initial = np.array([[1,0,0], [0,1,0]])      #Initial state 
    
    counter = 0                            # Stores the number of Metropolis moves (including burn-in)
    counter_samples = 0                    # Stores the number of saved samples (after burn in, used as index)
    counter_accepted = 0                   # Stores the number of accepted trial moves (after burn_in_period)
    counter_rejected = 0                   # Stores the number of rejected (after burn_in_period)

    pos = position_initial                                                                 # Stores the previous position
    prob_prev = probability_density(trial_wavefunction(pos, alpha, hamiltonian_label))     # Stores the previous probability density

    while counter_samples < num_samples:
        # Do trial move
        pos_trial = pos + d*np.random.uniform(-1, 1, size=pos.shape) 
        prob_trial = probability_density(trial_wavefunction(pos_trial, alpha, hamiltonian_label))

        if prob_trial > prob_prev or np.random.uniform(0,1) < prob_trial/prob_prev: # Condition for accepting trial move
            pos = pos_trial         # Update position
            prob_prev = prob_trial  # Update probability density
            if counter>=burn_in_period:
                counter_accepted += 1
            
        else:
            if counter>=burn_in_period:
                counter_rejected += 1
        
        if counter >= burn_in_period:   # After burn in period, we add the sampled positions to the output list
            if hamiltonian_label=="Helium":
                    positions_samples[:, :,counter-burn_in_period] = pos  
            else:
                    positions_samples[:,counter-burn_in_period] = pos  

            counter_samples += 1  

        counter += 1

    # We return the list with the Markov samples and the acceptance ratio, thinned out to account for correlations (data blocking)
    
    if hamiltonian_label=="Helium":
        positions_samples_thinned = positions_samples[:,:,::thinning_factor]
    else:
        positions_samples_thinned = positions_samples[:,::thinning_factor] 
    
    return positions_samples_thinned, counter_accepted/(counter_accepted + counter_rejected) 

    
def local_energy(positions_sample, alpha, hamiltonian_label):
    """
    This function computes the local energy at a specified position and for a given alpha.
    Additionally, it computes a custom expectation value, which is used for the optimization.
    
    Parameters
    ----------
    positions_sample: np.ndarray
        contains position sampled from MCMC samples
    alpha: np.ndarray
        array of variational parameters
    hamiltonian_label: str
        Name of system that is being minimized
        
    Returns
    -------
    E_loc: np.array
        array of local energies
    """

    # compute vectorized H|psi> analytically  

    if hamiltonian_label == 'QHO_1D':
        E_loc = alpha+positions_sample**2*(0.5-2*alpha**2)
        exp_val= (positions_sample**2)*E_loc # This is E_loc * (d ln(\psi)/d\alpha)
        
    elif hamiltonian_label == 'Hydrogen':
        # According to Thijssen, chap 12

        R = np.linalg.norm(positions_sample)        # compute distance to the origin of all the sampled positions
        E_loc = -1/R - 0.5*alpha*(alpha - 2/R)
        exp_val = R * E_loc                         # computes custom expectation value <R E_loc>

    elif hamiltonian_label == "Helium":

        # Here we compute E_loc and exp_val = E_loc * (d ln(\psi)/d\alpha), which is used for the gradient descent update [Thijssen]

        r1=np.linalg.norm(positions_sample[0], axis=0)
        r2=np.linalg.norm(positions_sample[1], axis=0)
        positions_12=positions_sample[0]-positions_sample[1]
        r12=np.linalg.norm(positions_12, axis=0)
        unit_vectors_12=positions_sample[0]/np.linalg.norm(positions_sample[0], axis=0)-positions_sample[1]/np.linalg.norm(positions_sample[1], axis=0)
        
        scalar_product=np.sum(positions_12*unit_vectors_12)

        E_loc= -4 + scalar_product/(r12*(1+alpha*r12)**2)-1/(r12*(1+alpha*r12)**3) - 1/(4*(1+alpha*r12)**4) + 1/r12
        exp_val = (r12**2/(2*(1+alpha*r12)**2))*E_loc 
        
    return E_loc, exp_val

def energy_function(positions_sample, alpha, hamiltonian_label):
    """
    This function computes the energy as a sum of local energies at the sampled positions.
    
    Parameters
    ----------
    positions_sample: np.ndarray
        contains degrees of freedom 
    alpha: np.ndarray
        array of variational parameters
    hamiltonian_label: str
        Name of system that is being minimized
        
    Returns
    -------
    energy_alpha: float
        position samples according to probability distribution
    exp_value: float 
        expectation value of custom_exp_value, used for the optimization
    """

    local_energies=[]
    local_energy_std = None
    custom_exp_values=[]
    
    # For a given alpha, we compute the local energy, the average energy and exp_value (term required for gradient descent)

    for position_sample in positions_sample:
        
        E_loc, exp_value = local_energy(position_sample, alpha, hamiltonian_label)
        local_energies.append(E_loc)
        custom_exp_values.append(exp_value)
    energy_alpha=np.mean(local_energies)
    local_energy_std = np.std(local_energies)
    exp_value=np.mean(custom_exp_values)
        
    return energy_alpha, local_energy_std, exp_value

def variational_method(alpha, energy_iter, convergence_crit, GD_stepsize, num_samples, burn_in_period, d, thinning_factor, hamiltonian_label):
    """
    This function implements the variational method, finding the optimal parameter alpha corresponding to the 
    ground state wave function.
    
    Parameters
    ----------
   
    alpha: np.ndarray
        array of variational parameters
    energy_iter: int
        Number of MCMC procedures per optimization step
    convergence_crit: float
        Criterion for achieving convergence
    GD_stepsize: float
        Stepsize for gradient descent update
    num_samples: int
        Number of samples per MCMC walker
    burn_in_period: int
        Period at the beginning of MCMC sampling which is to be discarded
    d: float
        MCMC trial step size
    thinning_factor: int
        Factor by which to thin out the resulting samples in order to mitigate correlations    
    hamiltonian_label: str
        Name of system that is being minimized
        
    Returns
    -------
    energies: np.ndarray
        Array of energies over different variational parameters
    energies_std: np.ndarray
        Array of energy standard deviations over different variational parameters
    alphas: np.ndarray
        Array of variational parameters found in variational algorithm
    """
    
    energies = []
    energies_std = []
    alphas = []
    counter = 0

    while True:
        # Initialise arrays of quantities to track for every iteration of the variational method
        energies_alpha = []
        energy_std_alpha = []
        exp_local_alpha = []
        exp_MCMC = []
        acceptance_ratios = []
        alphas.append(alpha)

        for iter in range(energy_iter):   
                # Obtain position samples from MCMC
                position_samples, acceptance_ratio = metropolis(alpha, num_samples, burn_in_period, d, thinning_factor, hamiltonian_label)
                
                # We distinguish three cases associated with three different physical systems: Hydrogen atom, 1D Quantum Harmonic Oscillator and Helium atom
                if hamiltonian_label == "Hydrogen":
                    
                    position_samples = np.stack((position_samples[0], position_samples[1], position_samples[2]), axis=1)
            
                    # Compute alpha-dependent energy, variance of local energies for fixed alpha and custom expectation value for optimisation
                    energy, local_energy_std, local_exp_val = energy_function(position_samples,alpha,hamiltonian_label)
                    exp_MCMC.append(np.mean(np.linalg.norm(position_samples,axis=1)))

                elif hamiltonian_label == "QHO_1D":

                    # Compute alpha-dependent energy, variance of local energies for fixed alpha and custom expectation value for optimisation
                    energy, local_energy_std, local_exp_val = energy_function(position_samples,alpha,hamiltonian_label)
                    exp_MCMC.append(np.mean(position_samples**2))

                elif hamiltonian_label=="Helium":

                    # Compute alpha-dependent energy, variance of local energies for fixed alpha and custom expectation value for optimisation
                    positions_sample= np.transpose(position_samples, (2, 0, 1))
                    energy, local_energy_std, local_exp_val = energy_function(positions_sample,alpha,hamiltonian_label)

                    MCMC_sublist=[]

                    for position_sample in positions_sample:
                        positions_12=position_sample[0]-position_sample[1]
                        r12=np.linalg.norm(positions_12, axis=0)
                        MCMC_sublist.append(r12**2/(2*(1+alpha*r12)**2))

                    exp_MCMC.append(np.mean(MCMC_sublist))

                # Store the energy of the current alpha, its standard deviation, the acceptance ratio of Metropolis
                # and the variable that stores the term local_exp_val = E_loc * (d ln(\psi)/d\alpha) required for gradient descent [Thijssen]
                energies_alpha.append(energy)
                energy_std_alpha.append(local_energy_std)
                acceptance_ratios.append(acceptance_ratio)
                exp_local_alpha.append(local_exp_val)
            
        print(f"Average acceptance ratio = {np.mean(acceptance_ratios)}")
        
        # Store the average energy and its standard deviation (over the energies obtained in each batch of size energy_iter)
        energies.append(np.mean(energies_alpha))
        energies_std.append(np.mean(energy_std_alpha))
        
        print(f"Energy for variational parameter alpha = {alpha} is {energies[-1]} with standard deviation {energies_std[-1]}")

        # Here we establish the convergence criterium to stop the simulation, which is different for Helium since the std does not go to zero
        # because the ground state is not exactly expressed via the trial wave function.
        if (hamiltonian_label=="Hydrogen" or hamiltonian_label=="QHO_1D") and (energies_std[-1] < convergence_crit):
            break 
        elif hamiltonian_label=="Helium":
            # For Helium the standard deviation won't converge as the wavefunctions are not exact solutions
            # If the energyy does not change (within the threshold given by convergence_crit) during 5 iterations, we consider that it has converged
            if counter==5:
                break     
            
            if len(energies)>3: # We need to add this as it only makes sense to calculate the relative difference once we have at least two stored energies
                relative_diff=abs((energies[-1]-energies[-2])/energies[-1]) # This represents how much the estimated energy changes for consecutive alphas (iterations)
                if relative_diff<convergence_crit:
                    counter+=1
                else:
                    counter=0
        

        # Update variational parameter alpha according to gradient descent update rule
        update_direction = 2*(np.mean(exp_MCMC)*np.mean(energies_alpha)-np.mean(exp_local_alpha))
            
        alpha = alpha - GD_stepsize * (update_direction)

    return energies, energies_std, alphas


def unit_conversion(value, std):
    """
    This function returns the energy in SI units.

    Parameters
    ----------
    value: float
        Energy in normalised units
    std: float
        Standard deviation in normalised units
    hamiltonian_label: str
        Name of system that is being minimized
        
    Returns
    -------
    value_SI: float
        Energy in SI units
    std_SI: float 
        Standard deviation in SI units
    """

    a0=5.29177211
    fs_cte=7.297352
    c=2.99792458
    h_bar=1.05457182
    eta=((fs_cte*h_bar*c)/a0)*(1e-18)
    
    value_SI = value*eta
    std_SI = std*eta

    return value_SI, std_SI 

def plotting(hamiltonian_label, alphas, energies, energies_std, ground_state_energy=None):
    """
    This function is used to plot the results obtained in the variational algorithm.
    
    Parameters
    ----------
   
    hamiltonian_label: str
        Name of system that is being minimized
    alphas: np.ndarray
        Array of variational parameters found in variational algorithm
    energies: np.ndarray
        Array of energies over different variational parameters
    energies_std: np.ndarray
        Array of energy standard deviations over different variational parameters
    ground_state_energy: float
        Ground state energy of the Hamiltonian
    Returns
    -------
    
    """
    
    # Define the labels, titles, and reference energy and alpha optimal values [Thijssen] for each Hamiltonian 
    if hamiltonian_label == 'QHO_1D':
        label_name = 'Energy'
        title_name = 'Quantum Harmonic Oscillator in 1D'
        ground_state_energy = 0.5
        variational_opt = 0.5
        
    elif hamiltonian_label == 'Hydrogen':
        label_name = 'Energy'
        title_name= 'Hydrogen atom'
        ground_state_energy = -0.5
        variational_opt = 1

    elif hamiltonian_label == 'Helium':
        label_name = 'Energy'
        title_name= 'Helium atom'
        variational_opt = alphas[energies.index(ground_state_energy)]

    energies_var = [e**2 for e in energies_std]

    # Figure 1: alpha-dependent energies
    plt.figure()
    plt.errorbar(alphas, energies, yerr=energies_var)
    plt.xlabel('Variational parameter (alpha)')
    plt.ylabel(label_name)
    plt.axhline(y=ground_state_energy, color='green', linestyle='--', label='Ground state energy')

    plt.title(f'Variational energy optimization of {title_name} using MCMC approach')
    plt.show()

    # Figure 2: Evolution of variational parameter with respect to the optimization iteration
    plt.figure()
    plt.plot(alphas, label='Variational parameter (alpha)')
    plt.xlabel('Iteration')
    plt.ylabel('alpha')
    plt.axhline(y=variational_opt, color='green', linestyle='--', label='Optimal parameter')

    if hamiltonian_label=="QHO_1D" or hamiltonian_label=="Hydrogen":
        plt.axhline(y=variational_opt, color='green', linestyle='--', label='Optimal parameter')
    elif hamiltonian_label=="Helium":
        plt.axhspan(variational_opt - 0.03, variational_opt + 0.03, color='grey', alpha=0.3)


    plt.title(f'Evolution of variational parameter')
    plt.show()
