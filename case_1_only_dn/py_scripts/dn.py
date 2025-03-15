import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def read_msd_file(filename):
    time_steps, msd_values = [], []
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('#') and line.strip():
                timestep, msd = map(float, line.split())
                time_steps.append(timestep)
                msd_values.append(msd)
    return np.array(time_steps), np.array(msd_values)

def calculate_diffusion_coefficient(time_steps, msd_values):
    slope, _, _, _, _ = linregress(time_steps, msd_values)
    return slope / 6  # D = slope / (2 * dimensions), where dimensions = 3

# List of MSD files and corresponding monomer numbers
msd_files = ['msd_12.dat', 'msd_16.dat', 'msd_24.dat', 'msd_32.dat']
monomer_numbers = [12, 16, 24, 32]

diffusion_coefficients = []

for file in msd_files:
    time_steps, msd_values = read_msd_file(file)
    D = calculate_diffusion_coefficient(time_steps, msd_values)
    diffusion_coefficients.append(D)

# Plot D vs n(monomer)
plt.figure(figsize=(10, 6))
plt.plot(monomer_numbers, diffusion_coefficients, 'bo-')
plt.xlabel('Number of Monomers (n)')
plt.ylabel('Diffusion Coefficient (D)')
plt.title('Diffusion Coefficient vs Number of Monomers')
plt.grid(True)

# Fit a power law: D = A * n^B
log_n = np.log(monomer_numbers)
log_D = np.log(diffusion_coefficients)
slope, intercept, _, _, _ = linregress(log_n, log_D)
A = np.exp(intercept)
B = slope

plt.plot(monomer_numbers, A * np.array(monomer_numbers)**B, 'r--', 
         label=f'Fit: D = {A:.2e} * n^{B:.2f}')
plt.legend()

plt.xscale('log')
plt.yscale('log')
plt.savefig('diffusion_coefficient.png')

plt.show()

print(f"Power law fit: D = {A:.2e} * n^{B:.2f}")
