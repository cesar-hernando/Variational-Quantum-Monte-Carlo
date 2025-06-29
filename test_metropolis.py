"""
(Temporal) File to test the implementation of metropolis algorithm
"""

import numpy as np
import matplotlib.pyplot as plt
import functions

hamiltonian_label = 'Hydrogen'

alpha = 1
num_samples = 10000
burn_in_period = 1000
d = 1.5
positions_samples, acceptance_ratio = functions.metropolis(alpha, num_samples, burn_in_period, d, hamiltonian_label)

print(f"Acceptance ratio = {acceptance_ratio}")

# Select what to plot dependidng of the Hamiltonian chosen
if hamiltonian_label == 'QHO_1D':
    samples = positions_samples[0,:]
    label_name = 'x (a.u)'
    title_name = 'Quantum Harmonic Oscillator in 1D'
elif hamiltonian_label == 'Hydrogen':
    samples = np.linalg.norm(positions_samples, axis=0)   # L2 norm (Euclidean distance to origin)
    label_name = 'r (a.u.)'
    title_name= 'Hydrogen atom'

# Plot the sampled position for each Markov step after the burn-in period
plt.figure()
plt.plot(range(num_samples), samples)
plt.xlabel('Markov Step')
plt.ylabel(label_name)
plt.title(f'Markov chain samples for {title_name}')
plt.show()

# Plot a histogram to see if the probability distribution matches the desired one
plt.figure()
plt.hist(samples, bins=100, density=True, alpha=0.7, color='blue')
plt.xlabel(label_name)
plt.ylabel('Probability Density')
plt.title(f'Histogram of Sampled Positions for {title_name}')
plt.show()
