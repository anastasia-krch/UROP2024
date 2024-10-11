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

initial = read('initial.traj')# need to rewrite from xyz to traj
final = read('final.traj')

# Generate intermediate images
n_images = 15  # General number of images - 1
images = [initial]
for i in range(1, n_images):
    image = initial.copy()
    image.calc = mace_off(model="medium", device="cuda", default_dtype="float64")
    images.append(image)
images.append(final)

# Create NEB object and interpolate to initialize the path
neb = NEB(images)
neb.interpolate()



# Run the NEB optimization
qn = BFGS(neb_refined, trajectory='neb.traj')
qn.run(fmax=0.05)
