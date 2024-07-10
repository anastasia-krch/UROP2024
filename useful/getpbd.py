import h5py
from ase import Atoms
from ase.visualize import view

# Function to create ASE Atoms object from HDF5 data
def create_atoms(symbols, positions):
    symbols = [symbol.decode('utf-8') for symbol in symbols]
    return Atoms(symbols=symbols, positions=positions)

# Open the HDF5 file
with h5py.File('trajectory_data.h5', 'r') as f:
    steps = list(f.keys())
    all_atoms = []
    
    for step in steps:
        group = f[step]
        positions = group['positions'][:]
        symbols = group['symbols'][:]
        atoms = create_atoms(symbols, positions)
        all_atoms.append(atoms)

# Visualize the first frame using ASE viewer
view(all_atoms)
