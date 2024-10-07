import sys
from ase.io import read
from mace.calculators import mace_mp
from ase.constraints import FixCom, FixAtoms

from ase.calculators.socketio import SocketClient

# Define atoms object
atoms = read("init.xyz", 0)
# Set ASE calculator
calc = mace_mp(model="medium", device="cuda", default_dtype="float32")

atoms.set_calculator(calc)

# Create Client
client = SocketClient(unixsocket="driver-works-100")
client.run(atoms, use_stress=True)
