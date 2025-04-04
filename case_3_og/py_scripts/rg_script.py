import os
import numpy as np
import matplotlib.pyplot as plt
import glob

def read_rg_file(file_path):
    """Read radius of gyration data from a LAMMPS output file."""
    try:
        with open(file_path, 'r') as f:
            # Load the data file
            data = np.loadtxt(file_path)
            
            # Check the shape of data to handle it correctly
            if isinstance(data, np.ndarray):
                if data.ndim == 0:  # Single value
                    return float(data)
                elif data.ndim == 1:  # 1D array
                    return float(data[-1])  # Return the last value
                elif data.ndim == 2:  # 2D array
                    # Assuming last column contains Rg values, and we want the last row
                    return float(data[-1, -1])
            else:
                return float(data)
                
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def main():
    # System sizes
    system_sizes = [12, 16, 20, 24, 28, 32]
    
    # Dictionary to store average Rg values for each system size
    avg_rg = {}
    std_rg = {}  # Also store standard deviation for error bars
    
    # Current directory (change this to the directory containing your data files)
    data_dir = "."
    
    # Process each system size
    for size in system_sizes:
        rg_values = []
        
        # Process each seed
        for seed in range(1, 5):  # Seeds: 1, 2, 3, 4
            # Construct file path pattern
            file_pattern = f"rg_all_{size}_{seed}*.dat"
            matching_files = glob.glob(os.path.join(data_dir, file_pattern))
            
            if matching_files:
                for file_path in matching_files:
                    rg = read_rg_file(file_path)
                    if rg is not None:
                        # Format the output based on the magnitude
                        if abs(rg) < 0.001 and rg != 0:
                            print(f"File {file_path}: Rg = {rg:.6e}")
                        else:
                            print(f"File {file_path}: Rg = {rg:.6f}")
                        rg_values.append(rg)
            else:
                print(f"No files found for system size {size}, seed {seed}")
        
        # Calculate average Rg for this system size if we have data
        if rg_values:
            avg_rg[size] = np.mean(rg_values)
            std_rg[size] = np.std(rg_values)
            
            # Format output based on magnitude
            if abs(avg_rg[size]) < 0.001 and avg_rg[size] != 0:
                print(f"System size {size}: Average Rg = {avg_rg[size]:.6e} ± {std_rg[size]:.6e}")
            else:
                print(f"System size {size}: Average Rg = {avg_rg[size]:.6f} ± {std_rg[size]:.6f}")
        else:
            print(f"No data available for system size {size}")
    
    # Create plot if we have data
    if avg_rg:
        # Sort by system size
        sizes = sorted(avg_rg.keys())
        rg_values = [avg_rg[size] for size in sizes]
        rg_errors = [std_rg[size] for size in sizes]
        
        plt.figure(figsize=(10, 6))
        
        # Plot with error bars
        plt.errorbar(sizes, rg_values, yerr=rg_errors, fmt='o-', linewidth=2, 
                     markersize=8, capsize=5, label='Data')
        
        plt.xlabel('System Size (N)', fontsize=14)
        plt.ylabel('Radius of Gyration (Rg)', fontsize=14)
        plt.title('Average Radius of Gyration vs. System Size', fontsize=16)
        plt.grid(True, alpha=0.3)
        
        # Use scientific notation for y-axis if values are very small
        if max([abs(v) for v in rg_values if v is not None]) < 0.001:
            plt.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
        
        # Optional: Fit power law (Rg ~ N^ν)
        if len(sizes) > 1:
            try:
                # Filter out zero or negative values for log calculation
                valid_indices = [i for i, v in enumerate(rg_values) if v > 0]
                if valid_indices:
                    valid_sizes = [sizes[i] for i in valid_indices]
                    valid_rg = [rg_values[i] for i in valid_indices]
                    
                    # Convert to log scale for linear fitting
                    log_sizes = np.log(valid_sizes)
                    log_rg = np.log(valid_rg)
                    
                    # Linear fit
                    fit = np.polyfit(log_sizes, log_rg, 1)
                    exponent = fit[0]
                    prefactor = np.exp(fit[1])
                    
                    # Plot fitted line
                    fit_x = np.linspace(min(sizes), max(sizes), 100)
                    fit_y = prefactor * fit_x**exponent
                    plt.plot(fit_x, fit_y, '--r', label=f'Fit: Rg = {prefactor:.4e} × N^{exponent:.4f}')
                    plt.legend(fontsize=12)
                    
                    print(f"Power law fit: Rg = {prefactor:.4e} × N^{exponent:.4f}")
                else:
                    print("Cannot perform power law fit - zero or negative values detected")
            except Exception as e:
                print(f"Error fitting power law: {e}")
        
        plt.tight_layout()
        plt.savefig('radius_of_gyration_vs_system_size.png', dpi=300)
        plt.show()
        print("Plot saved as 'radius_of_gyration_vs_system_size.png'")
    else:
        print("No data available for plotting")

if __name__ == "__main__":
    main()