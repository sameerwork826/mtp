
import matplotlib.pyplot as plt
import numpy as np

# Function to read data from file
def read_data(filename):
    time_steps, msd_values = [], []
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('#') and line.strip():
                timestep, msd = map(float, line.split())
                time_steps.append(timestep)
                msd_values.append(msd)
    return np.array(time_steps), np.array(msd_values)

# Read data
time_steps, msd_values = read_data('msd_20_10.dat')

# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot MSD in normal scale
ax1.plot(time_steps, msd_values, 'b-', label='MSD')
ax1.set_xlabel('Time Step')
ax1.set_ylabel('MSD')
ax1.set_title('MSD vs Time (Normal Scale)')
ax1.legend()
ax1.grid(True)

# Plot MSD in log-log scale
ax2.loglog(time_steps[1:], msd_values[1:], 'r-', label='MSD')
ax2.set_xlabel('Time Step')
ax2.set_ylabel('MSD')
ax2.set_title('MSD vs Time (Log-Log Scale)')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('msd_20_10.png')
plt.show()
