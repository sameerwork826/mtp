import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import glob

# List of MSD files (ensure filenames match actual data)
msd_files = ["msd_12_10.dat","msd_20_10.dat","msd_24_10.dat"]
monomer_counts = [12,20,24]  # Corresponding N values

diffusion_constants = []  # Store computed D values

for file, N in zip(msd_files, monomer_counts):
    # Load MSD data, skipping comment lines
    data = np.loadtxt(file, comments='#')

    time = data[:, 0]  # Time values
    msd = data[:, 1]   # MSD values

    # Fit the linear region: Avoid t = 0 (ballistic regime)
    start, end = 5, len(time)  # Adjust as needed
    slope, intercept, _, _, _ = linregress(time[start:end], msd[start:end])

    D = slope / 6  # Since MSD = 6Dt
    diffusion_constants.append(D)

    # Plot MSD vs. time with fit
    plt.plot(time, msd, label=f"N={N}")
    plt.plot(time[start:end], slope * time[start:end] + intercept, '--', label=f"Fit N={N}")

plt.xlabel("Time")
plt.ylabel("MSD")
plt.legend()
plt.title("MSD vs. Time for Different Polymer Chain Lengths")
plt.savefig("msd_vs_time.png")
plt.show()

# Convert lists to numpy arrays for fitting
N_values = np.array(monomer_counts)
D_values = np.array(diffusion_constants)

# Fit D vs. N using power law D ~ N^(-gamma)
log_N = np.log(N_values)
log_D = np.log(D_values)
gamma, log_prefactor, _, _, _ = linregress(log_N, log_D)

# Plot D vs. N
plt.plot(N_values, D_values, 'o-', label="Computed D")
plt.plot(N_values, np.exp(log_prefactor) * N_values**gamma, '--', label=f"Fit: D ~ N^{gamma:.2f}")
plt.xlabel("Number of Monomers (N)")
plt.ylabel("Diffusion Constant (D)")
plt.legend()
plt.title("Diffusion Constant vs. Polymer Chain Length")
plt.xscale("log")
plt.yscale("log")
plt.savefig("diffusion_constant_vs_chain_length.png")
plt.show()

print(f"Fitted relation: D ~ N^{gamma:.2f}")
