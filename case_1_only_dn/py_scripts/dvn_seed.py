import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from scipy.optimize import curve_fit

def extract_system_size(filename):
    """Extract system size from filename (e.g., msd_12_1.dat -> 12)"""
    match = re.search(r'msd_(\d+)_', filename)
    if match:
        return int(match.group(1))
    return None

def extract_seed(filename):
    """Extract seed number from filename (e.g., msd_12_1.dat -> 1)"""
    match = re.search(r'msd_\d+_(\d+)\.dat', filename)
    if match:
        return int(match.group(1))
    return None

def load_and_average_msd(directory, system_size):
    """Load all MSD files for a given system size and average them"""
    pattern = f"msd_{system_size}_*.dat"
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

def linear_fit(t, D):
    """Linear fit function for MSD(t) = 2*d*D*t, assuming d=3 dimensions"""
    return 6 * D * t  # 2*3*D*t for 3D

def calculate_diffusion_coefficient(time, msd, std_msd=None):
    """Calculate diffusion coefficient from MSD data"""
    # Filter out zero time point and ensure positive values
    valid_indices = time > 0
    time_valid = time[valid_indices]
    msd_valid = msd[valid_indices]
    
    # Use only the linear regime (typically early times)
    # We'll use the first third of the data points to avoid subdiffusive regime
    cutoff = len(time_valid) // 3
    time_fit = time_valid[:cutoff]
    msd_fit = msd_valid[:cutoff]
    
    if std_msd is not None:
        std_msd_valid = std_msd[valid_indices]
        std_msd_fit = std_msd_valid[:cutoff]
        # Use standard deviation as weights for the fit
        sigma = std_msd_fit
        # Replace zeros with small values to avoid division by zero
        sigma[sigma < 1e-10] = 1e-10
        weights = 1.0 / sigma
    else:
        weights = None
    
    try:
        # Fit the MSD curve to obtain D
        popt, pcov = curve_fit(linear_fit, time_fit, msd_fit, p0=[0.1], sigma=weights)
        D = popt[0]
        D_error = np.sqrt(pcov[0, 0])
        return D, D_error
    except Exception as e:
        print(f"Error fitting diffusion coefficient: {e}")
        return None, None

def analyze_diffusion_vs_size(directory, system_sizes):
    """Analyze diffusion coefficient vs. system size"""
    diffusion_coeffs = []
    diffusion_errors = []
    
    for size in system_sizes:
        time, avg_msd, std_msd = load_and_average_msd(directory, size)
        if time is not None:
            D, D_err = calculate_diffusion_coefficient(time, avg_msd, std_msd)
            if D is not None:
                print(f"System size {size}: D = {D:.6e} ± {D_err:.6e}")
                diffusion_coeffs.append(D)
                diffusion_errors.append(D_err)
            else:
                print(f"Failed to calculate D for system size {size}")
                diffusion_coeffs.append(np.nan)
                diffusion_errors.append(np.nan)
        else:
            diffusion_coeffs.append(np.nan)
            diffusion_errors.append(np.nan)
    
    return np.array(system_sizes), np.array(diffusion_coeffs), np.array(diffusion_errors)

def power_law(x, a, b):
    """Power law function: y = a * x^b"""
    return a * np.power(x, b)

def fit_power_law(x, y, yerr=None):
    """Fit power law to data"""
    # Remove NaN values
    valid = ~np.isnan(y)
    x_valid = x[valid]
    y_valid = y[valid]
    
    # Take log of both x and y for linear fit
    log_x = np.log(x_valid)
    log_y = np.log(y_valid)
    
    if yerr is not None:
        yerr_valid = yerr[valid]
        # Propagate errors to log space
        log_yerr = yerr_valid / y_valid
        weights = 1.0 / log_yerr
        weights[~np.isfinite(weights)] = 1.0
    else:
        weights = None
    
    # Linear fit in log space
    coeffs = np.polyfit(log_x, log_y, 1, w=weights)
    slope = coeffs[0]
    intercept = coeffs[1]
    
    # Convert back to power law parameters
    a = np.exp(intercept)
    b = slope
    
    # Calculate error using bootstrap method
    if yerr is not None:
        n_bootstrap = 1000
        slopes = []
        for _ in range(n_bootstrap):
            # Generate bootstrap sample
            indices = np.random.randint(0, len(log_x), len(log_x))
            x_boot = log_x[indices]
            y_boot = log_y[indices]
            w_boot = weights[indices] if weights is not None else None
            
            # Fit bootstrap sample
            try:
                boot_coeffs = np.polyfit(x_boot, y_boot, 1, w=w_boot)
                slopes.append(boot_coeffs[0])
            except:
                continue
        
        # Calculate standard deviation of bootstrap slope estimates
        b_err = np.std(slopes) if slopes else 0.0
    else:
        # Approximate error from the fit
        residuals = log_y - (slope * log_x + intercept)
        b_err = np.sqrt(np.sum(residuals**2) / (len(log_x) - 2)) / np.sqrt(np.sum((log_x - np.mean(log_x))**2))
    
    return a, b, b_err

def plot_diffusion_vs_size(system_sizes, diffusion_coeffs, diffusion_errors):
    """Plot diffusion coefficient vs. system size"""
    plt.figure(figsize=(10, 8))
    
    # Filter out NaN values
    valid = ~np.isnan(diffusion_coeffs)
    sizes_valid = system_sizes[valid]
    diff_valid = diffusion_coeffs[valid]
    err_valid = diffusion_errors[valid]
    
    # Plot data points with error bars
    plt.errorbar(sizes_valid, diff_valid, yerr=err_valid, fmt='o', markersize=8, 
                 color='blue', ecolor='black', capsize=5, label='Data')
    
    # Fit power law
    a, b, b_err = fit_power_law(sizes_valid, diff_valid, err_valid)
    print(f"Power law fit: D = {a:.6e} * N^({b:.4f} ± {b_err:.4f})")
    
    # Plot fit line
    x_fit = np.linspace(min(sizes_valid)*0.9, max(sizes_valid)*1.1, 100)
    y_fit = power_law(x_fit, a, b)
    plt.plot(x_fit, y_fit, 'r-', label=f'Fit: D ∝ N^{b:.4f}±{b_err:.4f}')
    
    # Set log scales
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel('System Size (N)', fontsize=14)
    plt.ylabel('Diffusion Coefficient (D)', fontsize=14)
    plt.title('Diffusion Coefficient vs. System Size', fontsize=16)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    plt.savefig('diffusion_vs_size.png', dpi=300)
    plt.show()

# Example usage
if __name__ == "__main__":
    # Set the directory where your MSD files are located
    data_directory = "."  # Current directory, change if needed
    
    # Define the system sizes you want to analyze
    system_sizes = [12, 16, 24, 32]
    
    # Analyze diffusion coefficient vs. system size
    sizes, diffusion_coeffs, diffusion_errors = analyze_diffusion_vs_size(data_directory, system_sizes)
    
    # Plot diffusion coefficient vs. system size
    plot_diffusion_vs_size(sizes, diffusion_coeffs, diffusion_errors)