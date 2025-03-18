import numpy as np

def write_polymer_data(filename, box_size, n_chains_a, n_chains_b, monomers_per_chain):
    """
    Generate LAMMPS data file for polymer chains with clear separation
    """
    # Calculate total number of atoms and bonds
    n_chains = n_chains_a + n_chains_b
    atoms_per_chain = monomers_per_chain
    total_atoms = n_chains * atoms_per_chain
    total_bonds = n_chains * (atoms_per_chain - 1)
    
    # Initialize arrays for coordinates
    positions = np.zeros((total_atoms, 3))
    molecule_ids = np.zeros(total_atoms, dtype=int)
    atom_types = np.zeros(total_atoms, dtype=int)
    
    # Refined spacing parameters
    monomer_spacing = 1.0  # Space between monomers in a chain
    chain_spacing_xy = 2.5  # Space between chains in xy plane
    
    # Calculate grid size for layout
    grid_size = int(np.ceil(np.sqrt(n_chains_a)))  # Grid for each type
    
    atom_idx = 0
    chain_idx = 0
    
    # Place type A chains (arranged in a grid on xy plane, alternating z levels)
    for i in range(n_chains_a):
        grid_x = i % grid_size
        grid_y = (i // grid_size) % grid_size
        grid_z = 0  # Type A starts at bottom
        
        # Calculate starting position
        start_x = grid_x * chain_spacing_xy
        start_y = grid_y * chain_spacing_xy
        start_z = grid_z * 4.0  # Vertical separation between types
        
        # Create linear chain
        for monomer in range(atoms_per_chain):
            positions[atom_idx] = [
                start_x,
                start_y,
                start_z + (monomer * monomer_spacing)
            ]
            molecule_ids[atom_idx] = chain_idx + 1
            atom_types[atom_idx] = 1  # Type A
            atom_idx += 1
        chain_idx += 1
    
    # Place type B chains (similar grid but shifted in z)
    for i in range(n_chains_b):
        grid_x = i % grid_size
        grid_y = (i // grid_size) % grid_size
        grid_z = 1  # Type B starts higher
        
        # Calculate starting position with offset
        start_x = grid_x * chain_spacing_xy + chain_spacing_xy/2  # Offset in x
        start_y = grid_y * chain_spacing_xy + chain_spacing_xy/2  # Offset in y
        start_z = grid_z * (atoms_per_chain * monomer_spacing + 4.0)  # Higher z level
        
        # Create linear chain
        for monomer in range(atoms_per_chain):
            positions[atom_idx] = [
                start_x,
                start_y,
                start_z + (monomer * monomer_spacing)
            ]
            molecule_ids[atom_idx] = chain_idx + 1
            atom_types[atom_idx] = 2  # Type B
            atom_idx += 1
        chain_idx += 1
    
    # Scale positions to fit within box size
    max_coords = np.max(positions, axis=0)
    scale_factors = np.array(box_size) / (max_coords + 2.0)  # Leave some margin
    positions = positions * scale_factors
    
    # Write data file
    with open(filename, 'w') as f:
        # Header
        f.write("LAMMPS data file for polymer chains\n\n")
        f.write(f"{total_atoms} atoms\n")
        f.write(f"{total_bonds} bonds\n\n")
        f.write("2 atom types\n")
        f.write("1 bond types\n\n")
        
        # Box size
        f.write(f"0.0 {box_size[0]} xlo xhi\n")
        f.write(f"0.0 {box_size[1]} ylo yhi\n")
        f.write(f"0.0 {box_size[2]} zlo zhi\n\n")
        
        # Masses
        f.write("Masses\n\n")
        f.write("1 1.0  # Type A\n")
        f.write("2 1.0  # Type B\n\n")
        
        # Atoms section
        f.write("Atoms # molecular\n\n")
        for i in range(total_atoms):
            f.write(f"{i+1} {molecule_ids[i]} {atom_types[i]} {positions[i][0]:.3f} {positions[i][1]:.3f} {positions[i][2]:.3f}\n")
        
        # Bonds section
        f.write("\nBonds\n\n")
        bond_id = 1
        for chain in range(n_chains):
            start_atom = chain * atoms_per_chain + 1
            for i in range(atoms_per_chain - 1):
                f.write(f"{bond_id} 1 {start_atom + i} {start_atom + i + 1}\n")
                bond_id += 1

# Parameters
box_size = [32.0, 32.0, 32.0]
n_chains_a =341  # Number of type A chains
n_chains_b =341# Number of type B chains
monomers_per_chain =12 # Monomers per chain

# Generate the data file
write_polymer_data("polymer_system_12_25.data", box_size, n_chains_a, n_chains_b, monomers_per_chain)