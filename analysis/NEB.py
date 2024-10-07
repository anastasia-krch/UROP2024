import ase
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import read
from ase.mep import NEB
from mace.calculators import mace_off
from ase.optimize import BFGS
from ase import Atoms
from ase.io import read, write

# Read initial and final structures
MD_file_name_initial = 'initial.xyz'
MD_file_name_final = 'final.xyz'
atoms_initial = ase.io.read(MD_file_name_initial, format='xyz')
atoms_final = ase.io.read(MD_file_name_final, format='xyz')

# Set up calculators for initial and final images
calc = mace_mp(model="medium", device="cuda", default_dtype="float64")

atoms_initial.calc = calc
atoms_final.calc = calc

# Optimize the initial and final structures
qn = BFGS(atoms_initial, trajectory='initial.traj')
qn.run(fmax=0.05)
qn = BFGS(atoms_final, trajectory='final.traj')
qn.run(fmax=0.05)

initial = read('initial.traj')#we need to rewrite from xyz to traj
final = read('final.traj')

# Generate intermediate images
n_images = 15  # General number of images
images = [initial]
for i in range(1, n_images):
    image = initial.copy()
    image.calc = mace_off(model="medium", device="cuda", default_dtype="float64")
    images.append(image)
images.append(final)

# Create NEB object and interpolate to initialize the path
neb = NEB(images)
neb.interpolate()

# Refine the path by inserting additional images between specific pairs
refined_images = []

n_refined = 5  # Number of additional images to add between each pair

for i in range(2):  # Loop over the first 2 pairs
    refined_images.append(images[i])  # Add the original image

    # Create and insert n_refined images between images[i] and images[i+1]
    for j in range(1, n_refined + 1):
        new_image = Atoms(
            symbols=images[i].get_chemical_symbols(),
            positions=(1 - j / (n_refined + 1)) * images[i].get_positions() + (j / (n_refined + 1)) * images[i + 1].get_positions(),
            cell=images[i].get_cell(),
            pbc=images[i].get_pbc()
        )
        new_image.calc = calc_final

        # Optimize each newly created image
        qn = BFGS(new_image, trajectory='new_image.traj')
        qn.run(fmax=0.05)
        new_image = read('new_image.traj')
        new_image.calc = mace_off(model="medium", device="cuda", default_dtype="float64")
        refined_images.append(new_image)

# Add the remaining images (from image 2 onwards)
refined_images.extend(images[2:])

# Create a new NEB object with the refined set of images
neb_refined = NEB(refined_images)
neb_refined.interpolate()

# Run the NEB optimization
qn = BFGS(neb_refined, trajectory='neb.traj')
qn.run(fmax=0.05)
