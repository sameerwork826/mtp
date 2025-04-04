# Molecular Dynamics Simulation of Polymer Systems

This repository contains molecular dynamics simulation data and analysis scripts for polymer systems, focusing on different cases of polymer behavior.

## Repository Structure

### Case 1: Diffusion coeficient relation calulation 
- Analysis of polymer systems with varying chain lengths (8, 12, 16, 24, 32)
- Mean Square Displacement (MSD) analysis
- Radius of Gyration (Rg) calculations
- Radial Distribution Function (RDF) data

### Case 2:Fixed Density for polymer 10 percent
- Analysis of polymer systems with varying chain lengths (8, 12, 16, 24, 32)
- Mean Square Displacement (MSD) analysis
- Radius of Gyration (Rg) calculations
- Radial Distribution Function (RDF) data

### Case 3: Charged Polymers
- Simulation of charged polymer systems
- Domain and cluster analysis
- RDF analysis for charged systems
- Coulombic interaction studies

## Tools and Technologies
- LAMMPS for molecular dynamics simulations
- Python scripts for data analysis and visualization
- Matplotlib for generating plots

## Analysis Capabilities
- Diffusion coefficient calculations
- Domain formation analysis
- Cluster size distribution
- Structural properties (Rg, RDF)

## Usage
Each case directory contains:
- LAMMPS input files
- Data files for polymer systems
- Python scripts for analysis
- Output data and visualization

To run the analysis scripts:
```bash
cd case_3_charged/py_scripts
python domains.py
```

## Requirements
- Python 3.12
- NumPy
- Matplotlib
- Ovito to run simulation dumps