from ase.io import read
from mace.calculators import mace_mp
from ase.calculators.socketio import SocketClient

# Define atoms object
atoms = read("init.xyz", 0)

# Set ASE calculator #################
calcs = []
calc = mace_mp(model="medium", device="cuda", default_dtype="float32")
atoms.set_calculator(calc) #equivalent to atoms.calc = macemp in https://mace-docs.readthedocs.io/en/latest/guide/foundation_models.html#pretrained-mace-mp-0-models

# Create Client
host = "driver"
client = SocketClient(unixsocket=host)
client.run(atoms, use_stress=True)
