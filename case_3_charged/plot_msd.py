# import numpy as np
# import matplotlib.pyplot as plt
# import sys

# # ===== CHANGE THIS VARIABLE TO MATCH YOUR FILE PATH =====
# # Example: input_file = "msd_data.txt"
# input_file = "msd_all.dat"  # <-- Change this to your file name
# # ======================================================

# def read_msd_data(filename):
#     """Read MSD data from a file, skipping comment lines."""
#     timesteps = []
#     msd_values = []
    
#     with open(filename, 'r') as file:
#         for line in file:
#             # Skip comment lines
#             if line.startswith('#'):
#                 continue
            
#             # Parse data lines
#             try:
#                 parts = line.strip().split()
#                 if len(parts) == 2:
#                     timestep = float(parts[0])
#                     msd = float(parts[1])
#                     timesteps.append(timestep)
#                     msd_values.append(msd)
#             except ValueError:
#                 continue
    
#     return np.array(timesteps), np.array(msd_values)

# def plot_msd(timesteps, msd_values, output_prefix="msd_plot"):
#     """Create normal and log-log plots of MSD data."""
#     # Create a figure with two subplots
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
#     # Plot 1: Normal scale
#     ax1.plot(timesteps, msd_values, 'o-', color='blue')
#     ax1.set_xlabel('Time Step')
#     ax1.set_ylabel('MSD')
#     ax1.set_title('Mean Square Displacement vs Time')
#     ax1.grid(True)
    
#     # Plot 2: Log-log scale
#     # Filter out zero values for log-log plot
#     nonzero_indices = msd_values > 0
#     timesteps_nonzero = timesteps[nonzero_indices]
#     msd_nonzero = msd_values[nonzero_indices]
    
#     ax2.loglog(timesteps_nonzero, msd_nonzero, 'o-', color='red')
#     ax2.set_xlabel('Time Step (log scale)')
#     ax2.set_ylabel('MSD (log scale)')
#     ax2.set_title('Mean Square Displacement vs Time (Log-Log)')
#     ax2.grid(True, which="both", ls="-")
    
#     # Fit a power law to the data (MSD ~ t^α)
#     if len(timesteps_nonzero) > 1:
#         log_t = np.log10(timesteps_nonzero[1:])  # Skip first point if it's zero
#         log_msd = np.log10(msd_nonzero[1:])
        
#         # Linear fit in log-log space
#         coeffs = np.polyfit(log_t, log_msd, 1)
#         alpha = coeffs[0]  # Diffusion exponent
        
#         # Plot the fit line
#         fit_y = 10**(coeffs[1]) * timesteps_nonzero**(alpha)
#         ax2.plot(timesteps_nonzero, fit_y, '--', color='green', 
#                 label=f'Fit: MSD ~ t^{alpha:.3f}')
#         ax2.legend()
    
#     plt.tight_layout()
#     plt.savefig(f"{output_prefix}.png", dpi=300)
#     plt.savefig(f"{output_prefix}.pdf")
    
#     print(f"Plots saved as {output_prefix}.png and {output_prefix}.pdf")
    
#     # Show the plot
#     plt.show()

# # You can also customize the output filename prefix here
# output_prefix = "msd_plot"  # Default output name

# def main():
#     # Use the global variable if no command line arguments are provided
#     global input_file, output_prefix
    
#     if len(sys.argv) > 1:
#         input_file = sys.argv[1]
#     if len(sys.argv) > 2:
#         output_prefix = sys.argv[2]
    
#     try:
#         print(f"Reading data from: {input_file}")
#         timesteps, msd_values = read_msd_data(input_file)
#         if len(timesteps) == 0:
#             print("No valid data found in the file.")
#             return
            
#         plot_msd(timesteps, msd_values, output_prefix)
#     except Exception as e:
#         print(f"Error: {e}")

# if __name__ == "__main__":
#     main()

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
time_steps, msd_values = read_data(r'C:\Users\nande\Desktop\mtp\case_2\msd_data\msd_32_10.dat')

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
plt.savefig('msd_32_10.png')
plt.show()
