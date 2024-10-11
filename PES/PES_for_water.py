from ase import Atoms
from ase.io import read, write
from mace.calculators import mace_mp
import numpy as np
from ase.optimize import BFGS
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation as R 

def rotate_water_molecule(water, angles):
    """Rotate the water molecule by given angles (in degrees)."""
    rotation = R.from_euler('xyz', angles, degrees=True)
    rotated_positions = rotation.apply(water.get_positions())
    water.set_positions(rotated_positions)
    return water

def objective_function(angles, atoms, fulleroid, water, x, y, z, calc):
    """Objective function to minimize: the potential energy after rotating the water molecule."""
    rotated_water = rotate_water_molecule(water.copy(), angles)
    rotated_water.translate([x, y, z])
    combined_atoms = rotated_water + fulleroid
    
    combined_atoms.set_calculator(calc)
    
    # Perform an energy calculation
    potential_energy = combined_atoms.get_potential_energy()
    
    return potential_energy

def calculate_energy_landscape(input_file, output_file, spacing=1.0, angle_increment=10.0):
    # Load the atomic structure (fullerene + water)
    atoms = read(input_file)

    # Assuming first three atoms are water
    water = atoms[:3]
    fulleroid = atoms[3:]

    # Determine the bounding box of the fullerene
    min_corner = np.min(fulleroid.get_positions(), axis=0)
    max_corner = np.max(fulleroid.get_positions(), axis=0)

    # Expand the box a bit beyond the fullerene structure
    margin = 3.0  
    min_corner -= margin
    max_corner += margin

    # box where the water molecule will be placed
    x_range = np.arange(min_corner[0], max_corner[0], spacing)
    y_range = np.arange(min_corner[1], max_corner[1], spacing)
    z_range = np.arange(min_corner[2], max_corner[2], spacing)

    # Set up a calculator for potential energy
    calc = mace_mp(model="medium", device="cuda", default_dtype="float32")
    atoms.set_calculator(calc)


    with open(output_file, 'w') as f:
        f.write("x y z rot_x rot_y rot_z PE\n")  # Header for the output file

        # Iterate over all grid points
        for x in x_range:
            for y in y_range:
                for z in z_range:
                    # Initial guess for the angles (0 degrees for all)
                    initial_angles = np.array([0.0, 0.0, 0.0])
                    result = minimize(objective_function, initial_angles, 
                                      args=(atoms, fulleroid, water, x, y, z, calc),
                                      method='BFGS')

                    # Extract the optimized angles and energy
                    best_rotation = result.x
                    lowest_energy = result.fun

                    f.write(f"{x} {y} {z} {lowest_energy}\n")

    print(f"Energy landscape data saved to {output_file}")

input_file = 'new_initial.xyz'  
output_file = 'energy_landscape.dat'
spacing = 0.5 

calculate_energy_landscape(input_file, output_file, spacing)
