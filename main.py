"""
This is the main file coordinating the project execution.
"""

import numpy as np
import functions
import parameters

############### PARAMETERS ##############################

# Retrieve parameters for quantum system and optimisation
hamiltonian_label, alpha, num_samples = parameters.hamiltonian_label, parameters.alpha, parameters.num_samples
burn_in_period, d, thinning_factor = parameters.burn_in_period, parameters.d, parameters.thinning_factor
energy_iter, GD_stepsize, convergence_crit = parameters.energy_iter, parameters.GD_stepsize, parameters.convergence_crit

############### VARIATIONAL METHOD ##############################

energies, energies_std, alphas = functions.variational_method(alpha, energy_iter, convergence_crit, GD_stepsize, num_samples, burn_in_period, d, thinning_factor, hamiltonian_label)

############### PLOTS ##############################

if hamiltonian_label=="QHO_1D":
    print(f"\n Ground state parameter is {alpha} with energy {energies[-1]}")
    functions.plotting(hamiltonian_label, alphas, energies, energies_std)
elif hamiltonian_label=="Hydrogen":
    print(f"\n Ground state parameter is {alpha} with energy {energies[-1]}")
    energy_SI, std_SI=functions.unit_conversion(energies[-1], energies_std[-1])
    print(f"\n In SI units, the Ground energy is {energy_SI}({std_SI})")
    functions.plotting(hamiltonian_label, alphas, energies, energies_std)
elif hamiltonian_label=="Helium":
    print(f"\n Ground state parameter is {alpha} with energy {min(energies)}")
    ground_state_energy = min(energies)
    ground_state_std = energies_std[energies.index(ground_state_energy)]
    energy_SI, std_SI=functions.unit_conversion(ground_state_energy, ground_state_std)
    print(f"\n In SI units, the Ground energy is {energy_SI}({std_SI})")
    functions.plotting(hamiltonian_label, alphas, energies, energies_std, ground_state_energy)