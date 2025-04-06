import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import re

def extract_system_size(filename):
    """Extract system size from filename (e.g., msd_all_12_1.dat -> 12)"""
    match = re.search(r'msd_all_(\d+)_', filename)
    if match:
        return int(match.group(1))
    return None

def extract_seed(filename):
    """Extract seed number from filename (e.g., msd_all_12_1.dat -> 1)"""
    match = re.search(r'msd_all_\d+_(\d+)\.dat', filename)
    if match:
        return int(match.group(1))
    return None

def load_and_average_msd(directory, system_size):
    """Load all MSD files for a given system size and average them"""
    pattern = f"msd_all_{system_size}_*.dat"
    files = glob.glob(os.path.join(directory, pattern))
    
    if not files:
        print(f"No files found for system size {system_size}")
        return None, None, None
    
    all_data = []
    
    # Load data from each file
    for file in files:
        try:
            # Use loadtxt with comments parameter to handle # lines
            data = np.loadtxt(file, comments='#')
            seed = extract_seed(file)
            print(f"Loaded {file} (seed {seed})")
            all_data.append(data)
        except Exception as e:
            print(f"Error loading {file}: {e}")
            continue
    
    if not all_data:
        print(f"No valid data loaded for system size {system_size}")
        return None, None, None
    
    # Check if all files have the same number of time points
    if not all(data.shape[0] == all_data[0].shape[0] for data in all_data):
        print(f"Warning: Not all files for system size {system_size} have the same number of time points")
        # Find the minimum length to truncate all arrays
        min_length = min(data.shape[0] for data in all_data)
        all_data = [data[:min_length] for data in all_data]
    
    # Stack data into a 3D array: [file_index, time_index, column_index]
    stacked_data = np.stack(all_data)
    
    # Extract time from the first column of the first file
    time = stacked_data[0, :, 0]
    
    # Extract MSD values from the second column of all files
    all_msd = stacked_data[:, :, 1]
    
    # Calculate average and standard deviation of MSD
    avg_msd = np.mean(all_msd, axis=0)
    std_msd = np.std(all_msd, axis=0)
    
    return time, avg_msd, std_msd

def plot_msd_normal(directory, system_sizes, plot_title="Mean Square Displacement"):
    """Plot average MSD for multiple system sizes in normal scale"""
    plt.figure(figsize=(10, 8))
    
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'pink']  # Colors for different system sizes
    markers = ['o', 's', '^', 'D', '<', ">"]  # Different markers for different system sizes
    
    for i, size in enumerate(system_sizes):
        time, avg_msd, std_msd = load_and_average_msd(directory, size)
        if time is not None:
            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]
            plt.plot(time, avg_msd, marker=marker, markersize=4, 
                     label=f"System Size {size}", color=color, linestyle='-', 
                     markevery=max(1, len(time)//20))  # Plot fewer markers for clarity
            plt.fill_between(time, avg_msd - std_msd, avg_msd + std_msd, alpha=0.2, color=color)
    
    plt.xlabel("Time", fontsize=14)
    plt.ylabel("MSD", fontsize=14)
    plt.title(f"{plot_title} (Normal Scale)", fontsize=16)
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig("msd_normal_scale.png", dpi=300)
    plt.show()

def plot_msd_log(directory, system_sizes, plot_title="Mean Square Displacement"):
    """Plot average MSD for multiple system sizes on log scale axes"""
    plt.figure(figsize=(10, 8))
    
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'pink']  # Colors for different system sizes
    markers = ['o', 's', '^', 'D', '<', '>']  # Different markers for different system sizes
    
    for i, size in enumerate(system_sizes):
        time, avg_msd, std_msd = load_and_average_msd(directory, size)
        if time is not None:
            # Filter out any zero or negative values for log scale
            valid_indices = (time > 0) & (avg_msd > 0)
            if np.any(valid_indices):
                valid_time = time[valid_indices]
                valid_msd = avg_msd[valid_indices]
                
                color = colors[i % len(colors)]
                marker = markers[i % len(markers)]
                
                plt.plot(valid_time, valid_msd, marker=marker, markersize=5, 
                         label=f"System Size {size}", color=color, linestyle='-',
                         markevery=max(1, len(valid_time)//15))  # Plot fewer markers for clarity
            else:
                print(f"Warning: No positive values for system size {size}, cannot plot on log scale")
    
    plt.xscale('log')  # Set x-axis to log scale
    plt.yscale('log')  # Set y-axis to log scale
    plt.xlabel("Time (log scale)", fontsize=14)
    plt.ylabel("MSD (log scale)", fontsize=14)
    plt.title(f"{plot_title} (Log Scale)", fontsize=16)
    plt.grid(True, which="both", ls="--")
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig("msd_log_scale.png", dpi=300)
    plt.show()

# Calculate diffusion coefficient
def calculate_diffusion(directory, system_sizes, time_range=None):
    """
    Calculate diffusion coefficient by fitting MSD = 2*D*t
    Optional time_range parameter to specify start and end times for fitting
    """
    plt.figure(figsize=(10, 8))
    
    diffusion_coeffs = []
    system_sizes_with_data = []
    
    for size in system_sizes:
        time, avg_msd, std_msd = load_and_average_msd(directory, size)
        if time is None:
            continue
            
        # If time range is specified, use only data within that range
        if time_range:
            start_time, end_time = time_range
            mask = (time >= start_time) & (time <= end_time)
            fit_time = time[mask]
            fit_msd = avg_msd[mask]
        else:
            fit_time = time
            fit_msd = avg_msd
            
        # Fit MSD = 2*D*t (linear fit through origin)
        # For 1D: MSD = 2*D*t
        coef = np.dot(fit_time, fit_msd) / np.dot(fit_time, fit_time)
        D = coef / 2  # Diffusion coefficient
        
        diffusion_coeffs.append(D)
        system_sizes_with_data.append(size)
        
        print(f"System Size {size}: D = {D:.6f}")
    
    # Plot diffusion coefficient vs system size
    plt.figure(figsize=(8, 6))
    plt.plot(system_sizes_with_data, diffusion_coeffs, 'o-', markersize=8)
    plt.xlabel("System Size", fontsize=14)
    plt.ylabel("Diffusion Coefficient (D)", fontsize=14)
    plt.title("Diffusion Coefficient vs System Size", fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("diffusion_vs_size.png", dpi=300)
    plt.show()
    
    return system_sizes_with_data, diffusion_coeffs

# Example usage
if __name__ == "__main__":
    # Set the directory where your MSD files are located
    data_directory = r"."  # Current directory, change if needed
    
    # Define the system sizes you want to analyze
    system_sizes = [12, 16, 20, 24, 28, 32]
    
    # Create separate normal and log plots
    plot_msd_normal(data_directory, system_sizes)
    plot_msd_log(data_directory, system_sizes)
    
    # Calculate diffusion coefficients (optional)
    # Uncomment to use - specify time range if needed
    # calculate_diffusion(data_directory, system_sizes, time_range=(10, 100))