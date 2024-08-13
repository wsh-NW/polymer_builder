import numpy as np
# USE: import Constants as cst

# TODO: Replace constants in scripts and import Constants

# GENERAL
box_boundaries = [0, 30, 0, 30, 0, 50]
box_lengths = [30, 30, 50]

# POLYMER
polymer_attempts = 10 # TODO: likely redundant -> in future just use chain attempts
chain_attempts = 100
carbon_attempts = 20

# CRYSTAL
lattice_spacing = 3.859 # Pt-Pt and Pt-Ni spacing
crystal_atom_depth = 3 # z lattice units deep
crystal_height = crystal_atom_depth * lattice_spacing # depth of crystal in angstroms
fcc_basis = np.array([[0.0, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.5, 0.0, 0.5],
                      [0.0, 0.5, 0.5]]) * lattice_spacing
bcc_basis = np.array([[0.0, 0.0, 0.0],
                      [0.5, 0.5, 0.5]]) * lattice_spacing

