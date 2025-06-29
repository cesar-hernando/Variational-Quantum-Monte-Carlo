"""
In this file the parameters used throughout the simulation are defined.
"""

# Choose hamiltonian_label from: 'QHO_1D', 'Hydrogen' and 'Helium'
hamiltonian_label = 'Helium'

# MCMC parameters
num_samples = 30000    # Number of total samples
burn_in_period = 4000  # Number of initial Markov steps that are not considered
thinning_factor = 2

'''
The values for d were empirically obtained by observing a ~0.5 acceptance ratio in the Metropolis algorithm.
We used the file test_metropolis.py to monitor the acceptance_ratio for different values of d

The values for alpha were chosen so to be far enough from the optimal alpha given by Thijssen to see convergence.

The variable convergence_crit establish the condition of convergence that stops the optimization loop. When the relative difference
of consecutive energies is less than convergence_crit for 5 consecutive loops, we consider that we have reached the ground state.
We distinguish three cases associated with three different physical systems: Hydrogen atom, 1D Quantum Harmonic Oscillator and Helium atom
'''

if hamiltonian_label=="QHO_1D":
    d=1.7 
    alpha=0.8
    convergence_crit = 1e-3
elif hamiltonian_label=="Hydrogen":
    d = 1.45 
    alpha=0.8
    convergence_crit = 1e-3
elif hamiltonian_label=="Helium":
    d=0.48 
    alpha = 0.12
    convergence_crit = 1e-3

# Number of batches of samples obtained for each alpha in the Metropolis algorithm in order to obtain the mean and standard deviation of the energy
energy_iter=20

# Learning rate used for gradient descent
GD_stepsize = 0.25