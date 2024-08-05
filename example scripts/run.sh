import os
import sys
from ase.io import read
from mace.calculators import mace_mp
from ase.calculators.socketio import SocketClient

# Get the driver name from the argument
driver = sys.argv[1]

# Define atoms object
atoms = read("init.xyz", 0)

# Set ASE calculator
calc = mace_mp(model="medium", device="cuda", default_dtype="float32")
atoms.calc = calc

# Create Client
client = SocketClient(unixsocket=f"{driver}")
client.run(atoms, use_stress=True)
